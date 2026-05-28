from __future__ import annotations
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple, Optional
import urllib.request
import urllib.parse
import urllib.error

from .. import base

# Reuse helper functions from redocking for conversion and pose parsing / RMSD
from .redocking import (
    _ensure_obabel,
    _convert_to_pdbqt,
    _parse_pdbqt_poses,
    _kabsch_rmsd,
    _write_ligand_ref_from_pdb,
    _write_ligand_ref_from_mmcif,
)


def _split_combined_pdbqt(combined: Path, out_dir: Path, prefix: str = "ligand_") -> List[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    parts: List[Path] = []
    with combined.open("r", encoding="utf-8", errors="ignore") as fh:
        buf: List[str] = []
        idx = 0
        for line in fh:
            buf.append(line)
            if line.strip() == "END":
                idx += 1
                fname = out_dir / f"{prefix}{idx:04d}.pdbqt"
                fname.write_text("".join(buf), encoding="utf-8")
                parts.append(fname)
                buf = []
        if buf:
            idx += 1
            fname = out_dir / f"{prefix}{idx:04d}.pdbqt"
            fname.write_text("".join(buf), encoding="utf-8")
            parts.append(fname)
    return parts


def _download_url_to_file(url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "DF-AUTOMATIZADO/1.0"})
    with urllib.request.urlopen(req) as resp:
        data = resp.read()
    out_path.write_bytes(data)


def _download_pubchem_cid(cid: str, out_path: Path) -> None:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{urllib.parse.quote(cid)}/SDF"
    _download_url_to_file(url, out_path)


def _download_chembl_id(chembl_id: str, out_path: Path) -> None:
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{urllib.parse.quote(chembl_id)}.sdf"
    try:
        _download_url_to_file(url, out_path)
    except urllib.error.HTTPError:
        url2 = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{urllib.parse.quote(chembl_id)}"
        _download_url_to_file(url2, out_path)


def _download_zinc_id(zinc_id: str, out_path: Path) -> None:
    tried = []
    patterns = [
        f"https://zinc15.docking.org/substances/{urllib.parse.quote(zinc_id)}/?format=sdf",
        f"https://zinc15.docking.org/substances/{urllib.parse.quote(zinc_id)}.sdf",
        f"https://zinc15.docking.org/substances/{urllib.parse.quote(zinc_id)}/download?format=sdf",
    ]
    for url in patterns:
        tried.append(url)
        try:
            _download_url_to_file(url, out_path)
            return
        except Exception:
            continue
    raise RuntimeError(f"Falha ao baixar {zinc_id} do ZINC. URLs tentadas: {tried}. Forneça os ligantes localmente ou verifique a API do ZINC.")


def _download_library_from_spec(spec: str, out_ligands: Path) -> Path:
    out_ligands.mkdir(parents=True, exist_ok=True)
    parts = spec.split(":", 1)
    if len(parts) != 2:
        raise ValueError("Spec de biblioteca inválida. Use 'pubchem:', 'chembl:' ou 'zinc:' prefixos")
    source, ids = parts[0].lower(), parts[1]
    ids_list = [i.strip() for i in ids.split(",") if i.strip()]
    combined = out_ligands / f"downloaded_{source}.sdf"
    with combined.open("wb") as outfh:
        for idx in ids_list:
            tmp = out_ligands / f"{source}_{idx}.sdf"
            try:
                if source == "pubchem":
                    _download_pubchem_cid(idx, tmp)
                elif source == "chembl":
                    _download_chembl_id(idx, tmp)
                elif source == "zinc":
                    _download_zinc_id(idx, tmp)
                else:
                    raise ValueError(f"Fonte desconhecida: {source}")
            except Exception as e:
                raise RuntimeError(f"Falha ao baixar {idx} de {source}: {e}")
            outfh.write(tmp.read_bytes())
            outfh.write(b"\n")
    return combined


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
                continue
            dest = out_ligands / f"{p.stem}.pdbqt"
            _convert_to_pdbqt(p, dest)
            lig_files.append(dest)
        return lig_files

    if library.suffix.lower() in {".sdf", ".mol2", ".smi"}:
        # Se o arquivo tem múltiplas moléculas, vamos extrair e converter cada uma
        # Estratégia: converter arquivo inteiro para SDF primeiro se necessário, depois para PDBQT
        import re
        temp_sdf = out_ligands / f"temp_combined.sdf"
        if library.suffix.lower() != ".sdf":
            cmd_to_sdf = [obabel, str(library), "-O", str(temp_sdf)]
            try:
                subprocess.run(cmd_to_sdf, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except Exception:
                temp_sdf = library
        else:
            temp_sdf = library
        
        # Agora parser o SDF para extrair moléculas individuais
        with temp_sdf.open("r", encoding="utf-8", errors="ignore") as fh:
            content = fh.read()
        
        # Split por "$$$$" (separador de moléculas em SDF) ou por "M  END"
        molecules = re.split(r"\$\$\$\$\s*\n", content)
        if len(molecules) <= 1:
            # Fallback para M  END
            molecules = re.split(r"M  END\s*\n", content)
        
        mol_count = 0
        for mol_str in molecules:
            if not mol_str.strip():
                continue
            # restaurar separador apropriado
            if "M  END" not in mol_str and "$$$$" not in mol_str:
                mol_str = mol_str.strip() + "\nM  END\n"
            mol_count += 1
            mol_sdf = out_ligands / f"mol_{mol_count}.sdf"
            mol_sdf.write_text(mol_str, encoding="utf-8")
            # converter para PDBQT
            mol_pdbqt = out_ligands / f"ligand_{mol_count:04d}.pdbqt"
            try:
                _convert_to_pdbqt(mol_sdf, mol_pdbqt)
                lig_files.append(mol_pdbqt)
            except Exception as e:
                # skip molecules que não convertem
                continue
        
        # Cleanup temp file if created
        if temp_sdf != library and temp_sdf.exists():
            try:
                temp_sdf.unlink()
            except Exception:
                pass
        
        return lig_files if lig_files else []

    if library.is_file():
        dest = out_ligands / f"{library.stem}.pdbqt"
        _convert_to_pdbqt(library, dest)
        return [dest]

    raise FileNotFoundError(f"Biblioteca de ligantes não encontrada: {library}")


def _search_pubchem_by_smiles(smiles: str, threshold: int = 90, max_records: int = 50) -> List[str]:
    """Pesquisa por similaridade no PubChem usando SMILES, retorna lista de CIDs (strings).

    Usa o endpoint PUG-REST similarity. Se a chamada falhar, lança RuntimeError.
    """
    # Encode smiles for URL
    esc = urllib.parse.quote(smiles, safe='')
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{esc}/cids/JSON?Threshold={threshold}&MaxRecords={max_records}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "DF-AUTOMATIZADO/1.0"})
        with urllib.request.urlopen(req) as resp:
            import json

            data = json.load(resp)
        cids = data.get("IdentifierList", {}).get("CID", [])
        return [str(c) for c in cids]
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Erro HTTP na busca PubChem: {e}")
    except Exception as e:
        raise RuntimeError(f"Erro ao consultar PubChem: {e}")


def _parse_vina_scores(pdbqt_file: Path) -> List[float]:
    """Extrai as afinidades reportadas pelo Vina no arquivo PDBQT de saída.

    Retorna lista de scores (float) na mesma ordem das poses reportadas.
    """
    import re

    scores: List[float] = []
    if not pdbqt_file.exists():
        return scores
    with pdbqt_file.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            # Vina usualmente escreve: "REMARK VINA RESULT:    -7.4     0.000     0.000"
            if line.startswith("REMARK VINA RESULT") or "VINA RESULT" in line:
                # extract floats
                floats = re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)
                if floats:
                    try:
                        scores.append(float(floats[0]))
                    except Exception:
                        continue
            # fallback: some versions write: "REMARK VINA: affinity = -7.4"
            if "affinity" in line and "=" in line:
                floats = re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)
                if floats:
                    try:
                        scores.append(float(floats[0]))
                    except Exception:
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
    """Faz busca por similaridade no PubChem (SMILES) e executa docking nos top N ligantes.

    Retorna lista de tuplas `(out_pdbqt_path, best_score, rmsd)`.
    """
    out = Path(output_dir)
    out_ligands = out / "ligands"
    out_ligands.mkdir(parents=True, exist_ok=True)

    cids = _search_pubchem_by_smiles(query_smiles, threshold=threshold, max_records=max_ligands)
    if not cids:
        raise RuntimeError("Nenhum CID retornado pela busca PubChem")
    spec = "pubchem:" + ",".join(cids[:max_ligands])
    # Reuse dock_library path by calling it with ligand_library spec
    results = dock_library(pdb_id, chain, spec, out, vina_exe=vina_exe, ligand_ref=ligand_ref)

    # For each result, compute best score and attach RMSD (results currently (path, rmsd))
    summarized: List[Tuple[Path, Optional[float], Optional[float]]] = []
    for out_file, rmsd in results:
        scores = _parse_vina_scores(out_file)
        best_score = min(scores) if scores else None
        summarized.append((out_file, best_score, rmsd))
    # write summary CSV
    import csv

    csv_file = out / "vina_results_summary.csv"
    with csv_file.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ligand_file", "best_score", "rmsd"])
        for p, s, r in summarized:
            writer.writerow([str(p), "" if s is None else f"{s:.3f}", "" if r is None else f"{r:.3f}"])

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
        if ligand_ref:
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
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append((x, y, z))
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
        else:
            raise ValueError("Para calcular a caixa automaticamente é necessário informar `ligand_ref` ou passar `center` e `size`")

    # suporta specs como 'pubchem:CID1,CID2', 'chembl:CHEMBL_ID', 'zinc:ZINC001'
    ligand_library_path = Path(ligand_library) if not isinstance(ligand_library, str) else None
    if isinstance(ligand_library, str):
        low = ligand_library.lower()
        if low.startswith("pubchem:") or low.startswith("chembl:") or low.startswith("zinc:"):
            ligand_library_path = _download_library_from_spec(ligand_library, out_ligands)
        else:
            ligand_library_path = Path(ligand_library)

    lig_pdbqt_files = _prepare_ligand_files(ligand_library_path, out_ligands)

    results: List[Tuple[Path, Optional[float]]] = []

    for lig in lig_pdbqt_files:
        ligand_name = lig.stem
        out_file = out_vina / f"{pdb_id.lower()}_{chain}_{ligand_name}_out.pdbqt"
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
        subprocess.run(cmd, check=True)

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
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                ref_coords.append((x, y, z))
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
                            r = _kabsch_rmsd(pose, ref_coords)
                        else:
                            r = _kabsch_rmsd(pose, ref_coords[: len(pose)])
                    except Exception:
                        r = float("inf")
                    if r < best:
                        best = r
                if best != float("inf"):
                    rmsd = best
        except Exception:
            rmsd = None

        results.append((out_file, rmsd))

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
def _download_pubchem_cid(cid: str, out_path: Path) -> None:
    # Use PubChem PUG REST to get SDF
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{urllib.parse.quote(cid)}/SDF"
    _download_url_to_file(url, out_path)


def _download_chembl_id(chembl_id: str, out_path: Path) -> None:
    # EBI ChEMBL API - try molecule endpoint returning SDF
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{urllib.parse.quote(chembl_id)}.sdf"
    try:
        _download_url_to_file(url, out_path)
    except urllib.error.HTTPError:
        # older API path fallback
        url2 = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{urllib.parse.quote(chembl_id)}"
        _download_url_to_file(url2, out_path)


def _download_zinc_id(zinc_id: str, out_path: Path) -> None:
    # ZINC access patterns vary; try common endpoints and provide informative error on failure
    tried = []
    patterns = [
        f"https://zinc15.docking.org/substances/{urllib.parse.quote(zinc_id)}/?format=sdf",
        f"https://zinc15.docking.org/substances/{urllib.parse.quote(zinc_id)}.sdf",
        f"https://zinc15.docking.org/substances/{urllib.parse.quote(zinc_id)}/download?format=sdf",
    ]
    for url in patterns:
        tried.append(url)
        try:
            _download_url_to_file(url, out_path)
            return
        except Exception:
            continue
    raise RuntimeError(f"Falha ao baixar {zinc_id} do ZINC. URLs tentadas: {tried}. Forneça os ligantes localmente ou verifique a API do ZINC.")


def _download_library_from_spec(spec: str, out_ligands: Path) -> Path:
    """Interpreta specs como `pubchem:CID1,CID2` ou `chembl:ID1,ID2` ou `zinc:ZINCID1`.
    Baixa SDFs individuais e concatena em um único SDF combinado retornando o caminho do arquivo.
    """
    out_ligands.mkdir(parents=True, exist_ok=True)
    parts = spec.split(":", 1)
    if len(parts) != 2:
        raise ValueError("Spec de biblioteca inválida. Use 'pubchem:', 'chembl:' ou 'zinc:' prefixos")
    source, ids = parts[0].lower(), parts[1]
    ids_list = [i.strip() for i in ids.split(",") if i.strip()]
    combined = out_ligands / f"downloaded_{source}.sdf"
    with combined.open("wb") as outfh:
        for idx in ids_list:
            tmp = out_ligands / f"{source}_{idx}.sdf"
            try:
                if source == "pubchem":
                    _download_pubchem_cid(idx, tmp)
                elif source == "chembl":
                    _download_chembl_id(idx, tmp)
                elif source == "zinc":
                    _download_zinc_id(idx, tmp)
                else:
                    raise ValueError(f"Fonte desconhecida: {source}")
            except Exception as e:
                raise RuntimeError(f"Falha ao baixar {idx} de {source}: {e}")
            # append to combined
            outfh.write(tmp.read_bytes())
            # ensure separator
            outfh.write(b"\n")
    return combined


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
    """Executa docking com AutoDock Vina para todos os ligantes de uma biblioteca.

    Retorna uma lista de tuplas `(out_pdbqt_path, rmsd)` onde `rmsd` é `None` quando
    não foi possível calcular (por exemplo se não foi fornecida referência).

    - `ligand_library` pode ser um diretório ou um arquivo (SDF multi-molécula também suportado).
    - `ligand_ref` (opcional): ID do ligante presente no arquivo de estrutura bruta para calcular box.
    - `center` e `size` podem ser informados manualmente; caso contrário, `ligand_ref` é necessário.
    """
    out = Path(output_dir)
    out_ligands = out / "ligands"
    out_vina = out / "vina_out"
    out_vina.mkdir(parents=True, exist_ok=True)

    # baixa estrutura bruta e extrai cadeia
    raw = base.download_structure(pdb_id, out)
    if raw.suffix.lower() == ".cif":
        receptor = base.extract_chain_from_mmcif(raw, chain, out, pdb_id)
    else:
        receptor = base.extract_chain_from_pdb(raw, chain, out, pdb_id)

    # receptor pdbqt (usa backup se disponível como no redocking)
    repo_root = Path(__file__).resolve().parents[2]
    backup_receptor_pdbqt = repo_root / "Backup" / "vina" / "receptor.pdbqt"
    if backup_receptor_pdbqt.exists() and pdb_id.upper() == "9THJ":
        receptor_pdbqt = backup_receptor_pdbqt
    else:
        receptor_pdbqt = out / "receptor" / f"{pdb_id.lower()}_{chain}.pdbqt"
        _convert_to_pdbqt(receptor, receptor_pdbqt)

    # calcula caixa se não dada
    if center is None or size is None:
        if ligand_ref:
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
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            coords.append((x, y, z))
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
        else:
            raise ValueError("Para calcular a caixa automaticamente é necessário informar `ligand_ref` ou passar `center` e `size`")

    # prepara arquivos de ligante
    # suporta specs como 'pubchem:CID1,CID2', 'chembl:CHEMBL_ID', 'zinc:ZINC001'
    ligand_library_path = Path(ligand_library) if not isinstance(ligand_library, str) else None
    if isinstance(ligand_library, str):
        low = ligand_library.lower()
        if low.startswith("pubchem:") or low.startswith("chembl:") or low.startswith("zinc:"):
            # baixa e cria um SDF combinado
            ligand_library_path = _download_library_from_spec(ligand_library, out_ligands)
        else:
            ligand_library_path = Path(ligand_library)

    lig_pdbqt_files = _prepare_ligand_files(ligand_library_path, out_ligands)

    results: List[Tuple[Path, Optional[float]]] = []

    # executa docking para cada ligante
    for lig in lig_pdbqt_files:
        ligand_name = lig.stem
        out_file = out_vina / f"{pdb_id.lower()}_{chain}_{ligand_name}_out.pdbqt"
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
        except subprocess.CalledProcessError as e:
            print(f"Aviso: docking falhou para {lig} (pode estar mal formatado). Pulando...", file=sys.stderr)
            continue

        # tenta calcular RMSD em relação ao ligand_ref (se disponível)
        rmsd: Optional[float] = None
        try:
            if ligand_ref:
                # extrai poses do arquivo de saída e compara com referência
                poses = _parse_pdbqt_poses(out_file)
                # read reference coords from previously created lig_ref_pdb
                lig_ref_pdb = out_ligands / f"{ligand_ref}_ref.pdb"
                ref_coords: List[Tuple[float, float, float]] = []
                with lig_ref_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        if not line.startswith(("ATOM  ", "HETATM")):
                            continue
                        try:
                            if len(line) >= 54:
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                ref_coords.append((x, y, z))
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
                            r = _kabsch_rmsd(pose, ref_coords)
                        else:
                            r = _kabsch_rmsd(pose, ref_coords[: len(pose)])
                    except Exception:
                        r = float("inf")
                    if r < best:
                        best = r
                if best != float("inf"):
                    rmsd = best
        except Exception:
            rmsd = None

        results.append((out_file, rmsd))

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Executa docking Vina para uma biblioteca de ligantes")
    parser.add_argument("pdb_id")
    parser.add_argument("chain")
    parser.add_argument("ligand_library", help="Diretório de ligantes ou arquivo SDF")
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
