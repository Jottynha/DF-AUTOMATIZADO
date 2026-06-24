from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

from . import analysis as dock6_analysis
from . import docking as dock6_docking
from .prepare import prepare_dock6_inputs


@dataclass(frozen=True)
class Dock6LibraryResult:
    ligand_label: str
    ligand_file: str | None
    run_dir: str | None
    scored_mol2: str | None
    dock_out: str | None
    grid_score: float | None
    rmsd: float | None
    rank: int | None
    status: str
    error: str | None = None


def _ensure_obabel(obabel_exe: str = "obabel") -> str:
    if shutil.which(obabel_exe):
        return obabel_exe
    for candidate in ("obabel", "babel", "obconvert"):
        if shutil.which(candidate):
            return candidate
    raise FileNotFoundError("Open Babel não foi encontrado no PATH.")


def _run(command: list[str], log_file: Path | None = None) -> None:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        log_file.write_text(
            "\n".join(part for part in [completed.stdout.strip(), completed.stderr.strip()] if part),
            encoding="utf-8",
        )
    if completed.returncode != 0:
        raise RuntimeError(
            "\n".join(
                [
                    f"Comando falhou: {' '.join(command)}",
                    completed.stdout.strip(),
                    completed.stderr.strip(),
                ]
            )
        )


def _safe_label(value: str, fallback: str = "ligand") -> str:
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._-")
    return label or fallback


def _download_url_to_file(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "DF-AUTOMATIZADO/1.0"})
    with urllib.request.urlopen(req) as resp:
        out_path.write_bytes(resp.read())


def _download_pubchem_cid(cid: str, out_path: Path) -> None:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{urllib.parse.quote(cid)}/SDF"
    _download_url_to_file(url, out_path)


def _download_chembl_id(chembl_id: str, out_path: Path) -> None:
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{urllib.parse.quote(chembl_id)}.sdf"
    try:
        _download_url_to_file(url, out_path)
    except urllib.error.HTTPError:
        fallback = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{urllib.parse.quote(chembl_id)}"
        _download_url_to_file(fallback, out_path)


def _download_zinc_id(zinc_id: str, out_path: Path) -> None:
    attempts = [
        f"https://zinc15.docking.org/substances/{urllib.parse.quote(zinc_id)}/?format=sdf",
        f"https://zinc15.docking.org/substances/{urllib.parse.quote(zinc_id)}.sdf",
        f"https://zinc15.docking.org/substances/{urllib.parse.quote(zinc_id)}/download?format=sdf",
    ]
    last_error: Exception | None = None
    for url in attempts:
        try:
            _download_url_to_file(url, out_path)
            return
        except Exception as exc:  # pragma: no cover
            last_error = exc
    raise RuntimeError(f"Falha ao baixar {zinc_id} do ZINC: {last_error}")


def _download_library_from_spec(spec: str, out_ligands: Path) -> tuple[Path, list[str]]:
    parts = spec.split(":", 1)
    if len(parts) != 2:
        raise ValueError("Spec inválida. Use pubchem:, chembl: ou zinc:")

    source, ids = parts[0].lower(), parts[1]
    ids_list = [item.strip() for item in ids.split(",") if item.strip()]
    if not ids_list:
        raise ValueError("Nenhum identificador informado na biblioteca.")

    combined = out_ligands / f"downloaded_{source}.sdf"
    out_ligands.mkdir(parents=True, exist_ok=True)

    with combined.open("wb") as outfh:
        for identifier in ids_list:
            tmp = out_ligands / f"{source}_{_safe_label(identifier)}.sdf"
            if source == "pubchem":
                _download_pubchem_cid(identifier, tmp)
            elif source == "chembl":
                _download_chembl_id(identifier, tmp)
            elif source == "zinc":
                _download_zinc_id(identifier, tmp)
            else:
                raise ValueError(f"Fonte desconhecida: {source}")

            data = tmp.read_bytes()
            outfh.write(data)
            if not data.endswith(b"\n"):
                outfh.write(b"\n")
            if not data.rstrip().endswith(b"$$$$"):
                outfh.write(b"$$$$\n")

    return combined, ids_list


def _split_sdf_records(sdf_content: str) -> list[str]:
    return [rec.strip() for rec in re.split(r"\$\$\$\$\s*(?:\n|$)", sdf_content) if rec.strip()]


def _convert_to_mol2(
    input_file: Path,
    output_file: Path,
    obabel_exe: str = "obabel",
    extra_args: Iterable[str] = (),
) -> Path:
    obabel = _ensure_obabel(obabel_exe)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    command = [obabel, str(input_file), "-O", str(output_file)]
    command.extend(extra_args)
    _run(command, log_file=output_file.with_suffix(output_file.suffix + ".log"))
    if not output_file.exists() or output_file.stat().st_size == 0:
        raise RuntimeError(f"Open Babel não gerou MOL2 válido: {output_file}")
    return output_file


def _prepare_ligand_files(ligand_library: str, out_ligands: Path, obabel_exe: str) -> list[Path]:
    """
    Prepara uma biblioteca de ligantes para DOCK6.

    Regra importante:
    - Se o ligante já estiver em MOL2, ele não deve ser convertido novamente.
      Isso evita saídas como ligand_0001.mol2.mol2.
    - SDF/SMI são convertidos para MOL2.
    - Diretórios podem conter .mol2, .sdf e .smi.
    """
    _ensure_obabel(obabel_exe)

    library_path = Path(ligand_library)
    out_ligands.mkdir(parents=True, exist_ok=True)

    allowed_suffixes = {".mol2", ".sdf", ".smi", ".smiles"}

    if str(ligand_library).lower().startswith(("pubchem:", "chembl:", "zinc:")):
        downloaded, _ids = _download_library_from_spec(ligand_library, out_ligands)
        source_files = [downloaded]
    elif library_path.is_dir():
        source_files = sorted(
            file for file in library_path.iterdir()
            if file.is_file() and file.suffix.lower() in allowed_suffixes
        )
    elif library_path.is_file():
        source_files = [library_path]
    else:
        raise FileNotFoundError(f"Biblioteca de ligantes não encontrada: {library_path}")

    if not source_files:
        raise RuntimeError(f"Nenhum ligante válido encontrado em: {ligand_library}")

    prepared_files: list[Path] = []
    counter = 1

    for source_file in source_files:
        suffix = source_file.suffix.lower()

        if suffix == ".mol2":
            # MOL2 já é o formato esperado pelo DOCK6.
            destination = out_ligands / source_file.name

            # Se a origem já é exatamente o destino, apenas reutiliza.
            if source_file.resolve() != destination.resolve():
                destination.write_bytes(source_file.read_bytes())

            text = destination.read_text(encoding="utf-8", errors="ignore")
            if "@<TRIPOS>MOLECULE" not in text or "@<TRIPOS>ATOM" not in text:
                raise RuntimeError(f"MOL2 inválido: {destination}")

            prepared_files.append(destination)
            counter += 1
            continue

        if suffix in {".smi", ".smiles"}:
            lines = [
                line.strip()
                for line in source_file.read_text(encoding="utf-8", errors="ignore").splitlines()
                if line.strip() and not line.strip().startswith("#")
            ]

            for line in lines:
                smiles = line.split()[0]
                label = _safe_label(line.split()[1]) if len(line.split()) > 1 else f"ligand_{counter:04d}"
                tmp_smi = out_ligands / f"{label}.smi"
                tmp_smi.write_text(smiles + "\n", encoding="utf-8")
                out_mol2 = out_ligands / f"{label}.mol2"
                _convert_to_mol2(tmp_smi, out_mol2, obabel_exe=obabel_exe)
                prepared_files.append(out_mol2)
                counter += 1
            continue

        if suffix == ".sdf":
            content = source_file.read_text(encoding="utf-8", errors="ignore")
            records = _split_sdf_records(content)

            # Se o SDF tem múltiplos registros, separa para preservar nomes individuais.
            if len(records) > 1:
                for record in records:
                    label = f"ligand_{counter:04d}"
                    tmp_sdf = out_ligands / f"{label}.sdf"
                    tmp_sdf.write_text(record.rstrip() + "\n$$$$\n", encoding="utf-8")
                    out_mol2 = out_ligands / f"{label}.mol2"
                    _convert_to_mol2(tmp_sdf, out_mol2, obabel_exe=obabel_exe)
                    prepared_files.append(out_mol2)
                    counter += 1
            else:
                label = _safe_label(source_file.stem)
                tmp_sdf = out_ligands / f"{label}.sdf"
                if source_file.resolve() != tmp_sdf.resolve():
                    tmp_sdf.write_text(content, encoding="utf-8")
                out_mol2 = out_ligands / f"{label}.mol2"
                _convert_to_mol2(tmp_sdf, out_mol2, obabel_exe=obabel_exe)
                prepared_files.append(out_mol2)
                counter += 1
            continue

        raise ValueError(f"Formato de ligante não suportado: {source_file}")

    return prepared_files


def _search_pubchem_by_smiles(smiles: str, threshold: int = 90, max_records: int = 50) -> list[str]:
    """
    Busca CIDs similares no PubChem usando fastsimilarity_2d.

    Usa saída TXT em vez de JSON para evitar JSONDecodeError em respostas grandes.
    Também tenta usar MaxRecords no próprio endpoint para evitar baixar milhares
    de CIDs quando o usuário pediu poucos ligantes.
    """
    esc = urllib.parse.quote(smiles, safe="")

    urls: list[str] = []
    if max_records is not None and max_records > 0:
        urls.append(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
            f"fastsimilarity_2d/smiles/{esc}/cids/TXT?Threshold={threshold}&MaxRecords={max_records}"
        )

    urls.append(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        f"fastsimilarity_2d/smiles/{esc}/cids/TXT?Threshold={threshold}"
    )

    last_error: Exception | None = None

    for url in urls:
        req = urllib.request.Request(url, headers={"User-Agent": "DF-AUTOMATIZADO/1.0"})

        try:
            with urllib.request.urlopen(req) as resp:
                text_response = resp.read().decode("utf-8", errors="ignore")
        except urllib.error.HTTPError as exc:
            last_error = exc
            # Se MaxRecords não for aceito por alguma versão do endpoint,
            # tenta a URL sem MaxRecords.
            continue

        cids: list[str] = []
        for raw_line in text_response.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if not line.isdigit():
                continue
            cids.append(line)

        if max_records is not None and max_records > 0:
            cids = cids[:max_records]

        if cids:
            return cids

    if last_error is not None:
        detail = ""
        if isinstance(last_error, urllib.error.HTTPError):
            try:
                detail = last_error.read().decode("utf-8", errors="ignore")
            except Exception:
                detail = ""
            raise RuntimeError(
                "Falha na busca de similares no PubChem via fastsimilarity_2d. "
                f"HTTP {last_error.code}. Threshold={threshold}. "
                "Teste com threshold maior, outro SMILES ou use pubchem-file:. "
                f"Detalhe: {detail[:500]}"
            ) from last_error

    raise RuntimeError(
        "Nenhum CID similar retornado pelo PubChem. "
        f"Threshold={threshold}. Tente reduzir o threshold ou usar outro SMILES."
    )


def _write_summary_csv(results: list[Dock6LibraryResult], csv_file: Path) -> None:
    csv_file.parent.mkdir(parents=True, exist_ok=True)
    with csv_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "ligand_label",
                "grid_score",
                "rmsd",
                "status",
                "scored_mol2",
                "dock_out",
                "run_dir",
                "ligand_file",
                "error",
            ],
        )
        writer.writeheader()
        for result in results:
            writer.writerow(asdict(result))


def dock_library(
    pdb_id: str,
    chain: str,
    ligand_library: str | Path,
    output_dir: str | Path,
    dock6_param_dir: str | Path,
    *,
    ligand_ref: str,
    pdb2pqr_exe: str = "pdb2pqr",
    obabel_exe: str = "obabel",
    dms_exe: str = "dms",
    sphgen_exe: str = "sphgen",
    sphere_selector_exe: str = "sphere_selector",
    showbox_exe: str = "showbox",
    grid_exe: str = "grid",
    dock6_exe: str = "dock6",
    ph: float = 7.4,
    sphere_radius: float = 10.0,
    calculate_rmsd_against_ref: bool = False,
    continue_on_error: bool = True,
) -> list[Dock6LibraryResult]:
    """
    Faz triagem de uma biblioteca de ligantes usando DOCK6.

    `ligand_ref` é o ligante cristalográfico usado para definir o sítio ativo
    via sphere_selector. Os ligantes da biblioteca são avaliados um a um.
    """
    out = Path(output_dir)
    out_ligands = out / "ligands"
    prepared_dir = out / "prepared_site"

    prepared_site = prepare_dock6_inputs(
        pdb_id=pdb_id,
        chain_id=chain,
        ligand_id=ligand_ref,
        output_dir=prepared_dir,
        pdb2pqr_exe=pdb2pqr_exe,
        obabel_exe=obabel_exe,
        dms_exe=dms_exe,
        ph=ph,
    )

    ligand_library = _expand_library_spec_from_file(ligand_library)
    ligand_library_text = str(ligand_library).strip().lower()
    if ligand_library_text in {"self", "ref", "reference", "crystal", "crystallographic"}:
        # Modo conveniente para testar a biblioteca com o próprio ligante cristalográfico.
        # Isso evita depender de um SDF gerado por outro pipeline.
        ligand_files = [Path(prepared_site.ligand_dock_mol2)]
    else:
        ligand_files = _prepare_ligand_files(ligand_library, out_ligands, obabel_exe=obabel_exe)

    if not ligand_files:
        raise RuntimeError("Nenhum ligante foi preparado para o DOCK6.")

    raw_results: list[Dock6LibraryResult] = []

    for index, ligand_file in enumerate(ligand_files, start=1):
        label = _safe_label(ligand_file.stem, fallback=f"ligand_{index:04d}")
        try:
            run_dir = dock6_docking.run_dock6_pipeline_for_ligand(
                ligand_label=label,
                ligand_source=ligand_file,
                receptor_pdb=Path(prepared_site.receptor_clean_pdb),
                receptor_ms=Path(prepared_site.receptor_ms),
                receptor_mol2_source=Path(prepared_site.receptor_fixed_mol2),
                ligand_ref_mol2_source=Path(prepared_site.ligand_ref_mol2),
                ligand_dock_mol2_source=ligand_file,
                base_output_dir=out,
                dock6_param_dir=Path(dock6_param_dir),
                obabel_exe=obabel_exe,
                sphgen_exe=sphgen_exe,
                sphere_selector_exe=sphere_selector_exe,
                showbox_exe=showbox_exe,
                grid_exe=grid_exe,
                dock6_exe=dock6_exe,
                sphere_radius=sphere_radius,
                prepare_only=False,
            )

            analysis_result = dock6_analysis.analyze_dock6_result(
                dock6_run_dir=run_dir,
                crystal_ligand_mol2=Path(prepared_site.ligand_ref_mol2) if calculate_rmsd_against_ref else None,
            )

            raw_results.append(
                Dock6LibraryResult(
                    ligand_label=label,
                    ligand_file=str(ligand_file),
                    run_dir=str(run_dir),
                    scored_mol2=analysis_result.get("scored_mol2"),
                    dock_out=analysis_result.get("dock_out"),
                    grid_score=analysis_result.get("grid_score"),
                    rmsd=analysis_result.get("rmsd"),
                    rank=None,
                    status="ok",
                )
            )

        except Exception as exc:
            failure = Dock6LibraryResult(
                ligand_label=label,
                ligand_file=str(ligand_file),
                run_dir=None,
                scored_mol2=None,
                dock_out=None,
                grid_score=None,
                rmsd=None,
                rank=None,
                status="error",
                error=str(exc),
            )
            raw_results.append(failure)
            if not continue_on_error:
                raise

    successful = sorted(
        [result for result in raw_results if result.grid_score is not None],
        key=lambda item: item.grid_score if item.grid_score is not None else float("inf"),
    )
    ranks = {result.ligand_label: rank for rank, result in enumerate(successful, start=1)}

    ranked_results = [
        Dock6LibraryResult(
            ligand_label=result.ligand_label,
            ligand_file=result.ligand_file,
            run_dir=result.run_dir,
            scored_mol2=result.scored_mol2,
            dock_out=result.dock_out,
            grid_score=result.grid_score,
            rmsd=result.rmsd,
            rank=ranks.get(result.ligand_label),
            status=result.status,
            error=result.error,
        )
        for result in raw_results
    ]

    _write_summary_csv(ranked_results, out / "dock6_results_summary.csv")
    (out / "dock6_results_summary.json").write_text(
        json.dumps([asdict(result) for result in ranked_results], indent=2),
        encoding="utf-8",
    )

    return ranked_results


def search_and_dock(
    pdb_id: str,
    chain: str,
    query_smiles: str,
    output_dir: str | Path,
    dock6_param_dir: str | Path,
    *,
    ligand_ref: str,
    threshold: int = 90,
    max_ligands: int = 20,
    **kwargs: Any,
) -> list[Dock6LibraryResult]:
    cids = _search_pubchem_by_smiles(query_smiles, threshold=threshold, max_records=max_ligands)
    if not cids:
        raise RuntimeError("Nenhum CID retornado pela busca PubChem.")

    return dock_library(
        pdb_id=pdb_id,
        chain=chain,
        ligand_library="pubchem:" + ",".join(cids[:max_ligands]),
        output_dir=output_dir,
        dock6_param_dir=dock6_param_dir,
        ligand_ref=ligand_ref,
        **kwargs,
    )


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--dock6-param-dir", required=True)
    parser.add_argument("--pdb2pqr-exe", default="pdb2pqr")
    parser.add_argument("--obabel-exe", default="obabel")
    parser.add_argument("--dms-exe", default="dms")
    parser.add_argument("--sphgen-exe", default="sphgen")
    parser.add_argument("--sphere-selector-exe", default="sphere_selector")
    parser.add_argument("--showbox-exe", default="showbox")
    parser.add_argument("--grid-exe", default="grid")
    parser.add_argument("--dock6-exe", default="dock6")
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--sphere-radius", type=float, default=10.0)
    parser.add_argument("--calculate-rmsd-against-ref", action="store_true")
    parser.add_argument("--stop-on-error", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Triagem de biblioteca com DOCK6.")
    parser.add_argument("pdb_id")
    parser.add_argument("chain")
    parser.add_argument("ligand_ref")
    parser.add_argument("ligand_library")
    parser.add_argument("output_dir")
    _add_common_args(parser)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        results = dock_library(
            pdb_id=args.pdb_id,
            chain=args.chain,
            ligand_ref=args.ligand_ref,
            ligand_library=args.ligand_library,
            output_dir=args.output_dir,
            dock6_param_dir=args.dock6_param_dir,
            pdb2pqr_exe=args.pdb2pqr_exe,
            obabel_exe=args.obabel_exe,
            dms_exe=args.dms_exe,
            sphgen_exe=args.sphgen_exe,
            sphere_selector_exe=args.sphere_selector_exe,
            showbox_exe=args.showbox_exe,
            grid_exe=args.grid_exe,
            dock6_exe=args.dock6_exe,
            ph=args.ph,
            sphere_radius=args.sphere_radius,
            calculate_rmsd_against_ref=args.calculate_rmsd_against_ref,
            continue_on_error=not args.stop_on_error,
        )
        print(json.dumps([asdict(result) for result in results], indent=2))
        print(f"Resumo CSV: {Path(args.output_dir) / 'dock6_results_summary.csv'}")
        return 0
    except Exception as exc:
        print(f"Erro durante triagem DOCK6: {exc}", file=sys.stderr)
        import traceback

        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())




def _read_ids_from_text_file(file_path: Path) -> list[str]:
    """
    Lê IDs de um arquivo texto.

    Aceita:
    - um ID por linha;
    - IDs separados por vírgula;
    - comentários iniciados por #.
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo de IDs não encontrado: {file_path}")

    ids: list[str] = []
    for raw_line in file_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        for item in line.split(","):
            item = item.strip()
            if item:
                ids.append(item)

    if not ids:
        raise ValueError(f"Nenhum ID encontrado em: {file_path}")

    return ids


def _expand_library_spec_from_file(ligand_library: str) -> str:
    """
    Expande specs baseadas em arquivo para o formato já suportado pelo pipeline.

    Exemplos:
    - pubchem-file:data/cids.txt -> pubchem:2244,2662,...
    - chembl-file:data/chembl_ids.txt -> chembl:CHEMBL25,...
    - zinc-file:data/zinc_ids.txt -> zinc:ZINC000...
    """
    lowered = ligand_library.lower()

    prefixes = {
        "pubchem-file:": "pubchem:",
        "chembl-file:": "chembl:",
        "zinc-file:": "zinc:",
    }

    for file_prefix, target_prefix in prefixes.items():
        if lowered.startswith(file_prefix):
            file_path = Path(ligand_library.split(":", 1)[1]).expanduser()
            ids = _read_ids_from_text_file(file_path)
            return target_prefix + ",".join(ids)

    return ligand_library


def prepare_library_only(
    ligand_library: str,
    output_dir: Path | str,
    *,
    obabel_exe: str = "obabel",
) -> list[Path]:
    """
    Baixa/converte uma biblioteca de ligantes sem executar docking.

    Aceita os mesmos formatos do modo library:
    arquivo .sdf, .mol2, .smi, diretório, pubchem:, chembl: ou zinc:.
    """
    ligand_library = _expand_library_spec_from_file(ligand_library)

    output_dir = Path(output_dir)
    ligands_dir = output_dir / "ligands"
    ligands_dir.mkdir(parents=True, exist_ok=True)

    ligand_files = _prepare_ligand_files(ligand_library, ligands_dir, obabel_exe=obabel_exe)

    import json
    manifest = []
    for index, ligand_file in enumerate(ligand_files, start=1):
        manifest.append(
            {
                "rank": index,
                "name": ligand_file.stem,
                "path": str(ligand_file),
                "format": ligand_file.suffix.lstrip(".").lower(),
            }
        )

    (output_dir / "ligands_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    with (output_dir / "ligands_manifest.csv").open("w", encoding="utf-8") as handle:
        handle.write("rank,name,path,format\n")
        for item in manifest:
            handle.write(f'{item["rank"]},{item["name"]},{item["path"]},{item["format"]}\n')

    return ligand_files



def prepare_similar_library_only(
    query_smiles: str,
    output_dir: Path | str,
    *,
    threshold: int = 90,
    max_ligands: int = 1000,
    obabel_exe: str = "obabel",
) -> list[Path]:
    """
    Busca compostos similares no PubChem por SMILES e prepara a biblioteca,
    sem executar docking.

    Exemplo:
    prepare_similar_library_only("CC(=O)OC1=CC=CC=C1C(=O)O", "data/lib", max_ligands=1000)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cids = _search_pubchem_by_smiles(query_smiles, threshold=threshold, max_records=max_ligands)

    if not cids:
        raise RuntimeError("Nenhum CID similar encontrado no PubChem.")

    (output_dir / "pubchem_similar_cids.txt").write_text(
        "\n".join(cids) + "\n",
        encoding="utf-8",
    )

    return prepare_library_only(
        "pubchem:" + ",".join(cids),
        output_dir,
        obabel_exe=obabel_exe,
    )
