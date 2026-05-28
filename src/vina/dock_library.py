from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import List, Optional, Tuple

from .. import base
from .redocking import (
    _convert_to_pdbqt,
    _ensure_obabel,
    _kabsch_rmsd,
    _parse_pdbqt_poses,
    _write_ligand_ref_from_mmcif,
    _write_ligand_ref_from_pdb,
)


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


def _download_library_from_spec(spec: str, out_ligands: Path) -> tuple[Path, List[str]]:
    parts = spec.split(":", 1)
    if len(parts) != 2:
        raise ValueError("Spec de biblioteca inválida. Use 'pubchem:', 'chembl:' ou 'zinc:'")

    source, ids = parts[0].lower(), parts[1]
    ids_list = [item.strip() for item in ids.split(",") if item.strip()]
    combined = out_ligands / f"downloaded_{source}.sdf"
    out_ligands.mkdir(parents=True, exist_ok=True)

    with combined.open("wb") as outfh:
        for identifier in ids_list:
            tmp = out_ligands / f"{source}_{identifier}.sdf"
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


def _split_sdf_records(sdf_content: str) -> List[str]:
    records = [rec.strip() for rec in re.split(r"\$\$\$\$\s*(?:\n|$)", sdf_content) if rec.strip()]
    if len(records) <= 1:
        records = [rec.strip() for rec in re.split(r"(?=\nM  END\s*\n)", sdf_content) if rec.strip()]
    return records


def _prepare_ligand_files(library: Path, out_ligands: Path) -> List[Path]:
    out_ligands.mkdir(parents=True, exist_ok=True)
    obabel = _ensure_obabel()
    lig_files: List[Path] = []

    if library.is_dir():
        for p in sorted(library.iterdir()):
            if not p.is_file():
                continue
            if p.suffix.lower() == ".pdbqt":
                lig_files.append(p)
            else:
                dest = out_ligands / f"{p.stem}.pdbqt"
                _convert_to_pdbqt(p, dest)
                lig_files.append(dest)
        return lig_files

    if library.suffix.lower() in {".sdf", ".mol2", ".smi"}:
        temp_sdf = library
        if library.suffix.lower() != ".sdf":
            temp_sdf = out_ligands / "temp_combined.sdf"
            cmd = [obabel, str(library), "-O", str(temp_sdf)]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        content = temp_sdf.read_text(encoding="utf-8", errors="ignore")
        records = _split_sdf_records(content)
        if not records:
            raise RuntimeError(f"Nenhum registro SDF encontrado em {library}")

        for idx, record in enumerate(records, 1):
            mol_sdf = out_ligands / f"mol_{idx}.sdf"
            mol_sdf.write_text(record.rstrip() + "\n$$$$\n", encoding="utf-8")
            mol_pdbqt = out_ligands / f"ligand_{idx:04d}.pdbqt"
            _convert_to_pdbqt(mol_sdf, mol_pdbqt)
            lig_files.append(mol_pdbqt)

        return lig_files

    if library.is_file():
        dest = out_ligands / f"{library.stem}.pdbqt"
        _convert_to_pdbqt(library, dest)
        return [dest]

    raise FileNotFoundError(f"Biblioteca de ligantes não encontrada: {library}")


def _search_pubchem_by_smiles(smiles: str, threshold: int = 90, max_records: int = 50) -> List[str]:
    esc = urllib.parse.quote(smiles, safe="")
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/"
        f"smiles/{esc}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "DF-AUTOMATIZADO/1.0"})
    try:
        with urllib.request.urlopen(req) as resp:
            data = json.load(resp)
        cids = data.get("IdentifierList", {}).get("CID", [])
        return [str(cid) for cid in cids]
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Erro HTTP na busca PubChem: {exc}") from exc


def _parse_vina_scores(pdbqt_file: Path) -> List[float]:
    scores: List[float] = []
    if not pdbqt_file.exists():
        return scores

    with pdbqt_file.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if "VINA RESULT" not in line and "affinity" not in line:
                continue
            floats = re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)
            if not floats:
                continue
            try:
                scores.append(float(floats[0]))
            except ValueError:
                continue
    return scores


def search_and_dock(
    pdb_id: str,
    chain: str,
    query_smiles: str,
    output_dir: str | Path,
    vina_exe: str = "vina",
    ligand_ref: Optional[str] = None,
    threshold: int = 90,
    max_ligands: int = 20,
) -> List[Tuple[Path, Optional[float], Optional[float]]]:
    out = Path(output_dir)
    cids = _search_pubchem_by_smiles(query_smiles, threshold=threshold, max_records=max_ligands)
    if not cids:
        raise RuntimeError("Nenhum CID retornado pela busca PubChem")

    results = dock_library(
        pdb_id=pdb_id,
        chain=chain,
        ligand_library="pubchem:" + ",".join(cids[:max_ligands]),
        output_dir=out,
        vina_exe=vina_exe,
        ligand_ref=ligand_ref,
    )

    summarized: List[Tuple[Path, Optional[float], Optional[float]]] = []
    for out_file, rmsd in results:
        scores = _parse_vina_scores(out_file)
        summarized.append((out_file, min(scores) if scores else None, rmsd))
    return summarized


def dock_library(
    pdb_id: str,
    chain: str,
    ligand_library: str | Path,
    output_dir: str | Path,
    vina_exe: str = "vina",
    ligand_ref: Optional[str] = None,
    center: Optional[Tuple[float, float, float]] = None,
    size: Optional[Tuple[float, float, float]] = None,
) -> List[Tuple[Path, Optional[float]]]:
    out = Path(output_dir)
    out_ligands = out / "ligands"
    out_vina = out / "vina_out"
    out_vina.mkdir(parents=True, exist_ok=True)

    raw = base.download_structure(pdb_id, out)
    if raw.suffix.lower() == ".cif":
        receptor = base.extract_chain_from_mmcif(raw, chain, out, pdb_id)
    else:
        receptor = base.extract_chain_from_pdb(raw, chain, out, pdb_id)

    repo_root = Path(__file__).resolve().parents[2]
    backup_receptor_pdbqt = repo_root / "Backup" / "vina" / "receptor.pdbqt"
    if backup_receptor_pdbqt.exists() and pdb_id.upper() == "9THJ":
        receptor_pdbqt = backup_receptor_pdbqt
    else:
        receptor_pdbqt = out / "receptor" / f"{pdb_id.lower()}_{chain}.pdbqt"
        _convert_to_pdbqt(receptor, receptor_pdbqt)

    if center is None or size is None:
        if not ligand_ref:
            raise ValueError(
                "Para calcular a caixa automaticamente é necessário informar `ligand_ref` ou passar `center` e `size`"
            )

        lig_ref_pdb = out_ligands / f"{ligand_ref}_ref.pdb"
        if raw.suffix.lower() == ".cif":
            _write_ligand_ref_from_mmcif(raw, ligand_ref, lig_ref_pdb)
        else:
            _write_ligand_ref_from_pdb(raw, ligand_ref, lig_ref_pdb)

        coords: List[Tuple[float, float, float]] = []
        with lig_ref_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                try:
                    if len(line) >= 54:
                        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                    else:
                        parts = line.split()
                        if len(parts) >= 6:
                            coords.append((float(parts[-3]), float(parts[-2]), float(parts[-1])))
                except Exception:
                    continue

        if not coords:
            raise RuntimeError("Não foi possível ler coordenadas do ligante de referência para definir caixa")

        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]
        cx = sum(xs) / len(xs)
        cy = sum(ys) / len(ys)
        cz = sum(zs) / len(zs)
        dx = max(xs) - min(xs)
        dy = max(ys) - min(ys)
        dz = max(zs) - min(zs)
        center = (cx, cy, cz)
        size = (max(8.0, dx + 10.0), max(8.0, dy + 10.0), max(8.0, dz + 10.0))

    pubchem_ids: List[str] = []
    ligand_library_path = Path(ligand_library) if not isinstance(ligand_library, str) else None
    if isinstance(ligand_library, str):
        low = ligand_library.lower()
        if low.startswith(("pubchem:", "chembl:", "zinc:")):
            ligand_library_path, pubchem_ids = _download_library_from_spec(ligand_library, out_ligands)
        else:
            ligand_library_path = Path(ligand_library)

    lig_pdbqt_files = _prepare_ligand_files(ligand_library_path, out_ligands)

    results: List[Tuple[Path, Optional[float]]] = []
    summary_rows: List[Tuple[str, Optional[float], Optional[float], str]] = []

    for idx, lig in enumerate(lig_pdbqt_files):
        out_file = out_vina / f"{pdb_id.lower()}_{chain}_{lig.stem}_out.pdbqt"
        cmd = [
            vina_exe,
            "--receptor",
            str(receptor_pdbqt),
            "--ligand",
            str(lig),
            "--center_x",
            str(center[0]),
            "--center_y",
            str(center[1]),
            "--center_z",
            str(center[2]),
            "--size_x",
            str(size[0]),
            "--size_y",
            str(size[1]),
            "--size_z",
            str(size[2]),
            "--out",
            str(out_file),
        ]
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError:
            print(f"Aviso: docking falhou para {lig}. Pulando...", file=sys.stderr)
            continue

        rmsd: Optional[float] = None
        try:
            if ligand_ref:
                poses = _parse_pdbqt_poses(out_file)
                lig_ref_pdb = out_ligands / f"{ligand_ref}_ref.pdb"
                ref_coords: List[Tuple[float, float, float]] = []
                with lig_ref_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        if not line.startswith(("ATOM  ", "HETATM")):
                            continue
                        try:
                            if len(line) >= 54:
                                ref_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                            else:
                                parts = line.split()
                                if len(parts) >= 6:
                                    ref_coords.append((float(parts[-3]), float(parts[-2]), float(parts[-1])))
                        except Exception:
                            continue

                best = float("inf")
                for pose in poses:
                    try:
                        if len(pose) == len(ref_coords):
                            rms = _kabsch_rmsd(pose, ref_coords)
                        else:
                            rms = _kabsch_rmsd(pose, ref_coords[: len(pose)])
                    except Exception:
                        rms = float("inf")
                    if rms < best:
                        best = rms
                if best != float("inf"):
                    rmsd = best
        except Exception:
            rmsd = None

        results.append((out_file, rmsd))

        scores = _parse_vina_scores(out_file)
        best_score = min(scores) if scores else None
        pubchem_id = pubchem_ids[idx] if idx < len(pubchem_ids) else ""
        summary_rows.append((pubchem_id, best_score, rmsd, out_file.name))

    summary_rows_sorted = sorted(
        summary_rows,
        key=lambda row: row[1] if row[1] is not None else float("inf"),
    )

    csv_file = out / "vina_results_summary.csv"
    with csv_file.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["pubchem_id", "score", "rmsd", "rank"])
        for rank, (pubchem_id, best_score, rmsd, _ligand_file) in enumerate(summary_rows_sorted, 1):
            writer.writerow([
                pubchem_id,
                "" if best_score is None else f"{best_score:.3f}",
                "" if rmsd is None else f"{rmsd:.3f}",
                rank,
            ])

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Executa docking Vina para uma biblioteca de ligantes")
    parser.add_argument("pdb_id")
    parser.add_argument("chain")
    parser.add_argument("ligand_library", help="Diretório de ligantes ou arquivo SDF ou spec pubchem:/chembl:/zinc:")
    parser.add_argument("output_dir")
    parser.add_argument("--vina-exe", default="vina")
    parser.add_argument("--ligand-ref", help="ID do ligante de referência presente na estrutura bruta para definir caixa")
    args = parser.parse_args()
    res = dock_library(args.pdb_id, args.chain, args.ligand_library, args.output_dir, args.vina_exe, ligand_ref=args.ligand_ref)
    for p, r in res:
        if r is None:
            print(f"{p}: RMSD=N/A")
        else:
            print(f"{p}: RMSD={r:.3f} Å")
