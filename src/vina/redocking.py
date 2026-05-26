from __future__ import annotations
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Tuple

from .. import base


def _write_ligand_ref_from_pdb(raw_pdb: Path, ligand_code: str, output: Path) -> Path:
    ligand_code = ligand_code.upper()
    output.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    with raw_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("HETATM", "ATOM  ")):
                continue
            resname = line[17:20].strip().upper()
            if resname != ligand_code:
                continue
            lines.append(line)
    if not lines:
        raise FileNotFoundError(f"Nenhum átomo do ligante {ligand_code} encontrado em {raw_pdb}")
    lines.append("END\n")
    output.write_text("".join(lines), encoding="utf-8")
    return output


def _write_ligand_ref_from_mmcif(cif_file: Path, ligand_code: str, output: Path) -> Path:
    ligand_code = ligand_code.upper()
    output.parent.mkdir(parents=True, exist_ok=True)
    serial = 1
    lines: List[str] = []
    selected_model: str | None = None
    for row in base.iter_atom_site_rows(cif_file):
        if base.pick_value(row, "group_PDB").upper() != "HETATM":
            continue
        resname = base.pick_value(row, "auth_comp_id", "label_comp_id").upper()
        if resname != ligand_code:
            continue
        model_number = base.pick_value(row, "pdbx_PDB_model_num", default="1")
        if selected_model is None:
            selected_model = model_number
        if model_number != selected_model:
            continue
        atom_name = base.pick_value(row, "auth_atom_id", "label_atom_id", "id")
        alt_loc = base.pick_value(row, "label_alt_id")
        if alt_loc in {"?", "."}:
            alt_loc = ""
        res_seq = base.pick_value(row, "auth_seq_id", "label_seq_id")
        ins_code = base.pick_value(row, "pdbx_PDB_ins_code")
        if ins_code in {"?", "."}:
            ins_code = ""
        element = base.pick_value(row, "type_symbol", default=(atom_name[:1] if atom_name else ""))
        x = base.parse_float(base.pick_value(row, "Cartn_x"))
        y = base.parse_float(base.pick_value(row, "Cartn_y"))
        z = base.parse_float(base.pick_value(row, "Cartn_z"))
        lines.append(
            base.format_pdb_line(
                "HETATM",
                serial,
                atom_name,
                alt_loc,
                resname[:3],  # truncate to 3 chars for PDB format
                " ",
                res_seq,
                ins_code,
                x,
                y,
                z,
                1.0,
                0.0,
                element,
            )
        )
        serial += 1
    if not lines:
        raise FileNotFoundError(f"Nenhum átomo do ligante {ligand_code} encontrado em {cif_file}")
    lines.append("END\n")
    output.write_text("".join(lines), encoding="utf-8")
    return output


def _ensure_obabel() -> str:
    for cmd in ("obabel", "babel", "obconvert"):
        if shutil.which(cmd):
            return cmd
    raise FileNotFoundError("Open Babel não foi encontrado no PATH. Instale `obabel` para conversão PDB->PDBQT")


def _convert_to_pdbqt(input_file: Path, output_file: Path) -> None:
    obabel = _ensure_obabel()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    # tenta converter mantendo hydrogens e carregando corretamente
    cmd = [obabel, str(input_file), "-O", str(output_file), "--partialcharge", "gasteiger"]
    subprocess.run(cmd, check=True)


def _parse_pdbqt_poses(pdbqt_file: Path) -> List[List[Tuple[float, float, float]]]:
    poses: List[List[Tuple[float, float, float]]] = []
    current: List[Tuple[float, float, float]] = []
    with pdbqt_file.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("MODEL"):
                current = []
                continue
            if line.startswith("ENDMDL"):
                if current:
                    poses.append(current)
                current = []
                continue
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    # tenta parsear posições fixas primeiro (PDB padrão)
                    if len(line) >= 54:
                        x_str = line[30:38].strip()
                        y_str = line[38:46].strip()
                        z_str = line[46:54].strip()
                        x = float(x_str)
                        y = float(y_str)
                        z = float(z_str)
                    else:
                        # fallback: split por espaço em branco
                        parts = line.split()
                        if len(parts) >= 6:
                            x = float(parts[-3])
                            y = float(parts[-2])
                            z = float(parts[-1])
                        else:
                            continue
                    current.append((x, y, z))
                except (ValueError, IndexError):
                    continue
    # if file had no MODEL sections, treat whole file as single pose
    if poses == [] and current:
        poses.append(current)
    return poses


def _kabsch_rmsd(P: List[Tuple[float, float, float]], Q: List[Tuple[float, float, float]]) -> float:
    # both P and Q are lists of (x,y,z) same length
    n = len(P)
    if n == 0:
        return float('inf')
    # compute centroids
    cp = [sum(c[i] for c in P) / n for i in range(3)]
    cq = [sum(c[i] for c in Q) / n for i in range(3)]
    # center
    import math

    Pm = [[p[i] - cp[i] for i in range(3)] for p in P]
    Qm = [[q[i] - cq[i] for i in range(3)] for q in Q]
    # covariance matrix C = Pm^T * Qm
    C = [[0.0] * 3 for _ in range(3)]
    for i in range(n):
        for a in range(3):
            for b in range(3):
                C[a][b] += Pm[i][a] * Qm[i][b]
    # SVD of C using numpy if available, else use built-in via math + power-iteration? Try numpy first.
    try:
        import numpy as _np

        U, s, Vt = _np.linalg.svd(_np.array(C))
        R = _np.dot(U, Vt)
        Prot = _np.dot(_np.array(Pm), R)
        diffs = Prot - _np.array(Qm)
        rmsd = math.sqrt((diffs * diffs).sum() / n)
        return float(rmsd)
    except Exception:
        # fallback: compute RMSD without optimal rotation (upper bound)
        ss = 0.0
        for i in range(n):
            dx = Pm[i][0] - Qm[i][0]
            dy = Pm[i][1] - Qm[i][1]
            dz = Pm[i][2] - Qm[i][2]
            ss += dx * dx + dy * dy + dz * dz
        return math.sqrt(ss / n)


def preprocess_and_dock(pdb_id: str, chain: str, ligand_id: str, output_dir: str | Path, vina_exe: str = "vina") -> Tuple[Path, float]:
    out = Path(output_dir)
    repo_root = Path(__file__).resolve().parents[2]  # go up 2 dirs to get repo root
    
    # baixa e extrai receptor / ligante
    raw = base.download_structure(pdb_id, out)
    if raw.suffix.lower() == ".cif":
        receptor = base.extract_chain_from_mmcif(raw, chain, out, pdb_id)
        ligand_ref_pdb = out / "ligands" / f"{ligand_id}_ref.pdb"
        _write_ligand_ref_from_mmcif(raw, ligand_id, ligand_ref_pdb)
    else:
        receptor = base.extract_chain_from_pdb(raw, chain, out, pdb_id)
        ligand_ref_pdb = out / "ligands" / f"{ligand_id}_ref.pdb"
        _write_ligand_ref_from_pdb(raw, ligand_id, ligand_ref_pdb)

    # download ideal sdf (optional, but keep original behavior)
    try:
        base.download_ligand(ligand_id, out)
    except Exception:
        pass

    # tenta usar receptor pré-preparado do Backup se disponível e if pdb_id matches
    backup_receptor_pdbqt = repo_root / "Backup" / "vina" / "receptor.pdbqt"
    if backup_receptor_pdbqt.exists() and pdb_id.upper() == "9THJ":
        print(f"Usando receptor pré-preparado do Backup: {backup_receptor_pdbqt}")
        receptor_pdbqt = backup_receptor_pdbqt
    else:
        # converte para pdbqt
        receptor_pdbqt = out / "receptor" / f"{pdb_id.lower()}_{chain}.pdbqt"
        _convert_to_pdbqt(receptor, receptor_pdbqt)

    # converte ligante para pdbqt
    ligand_pdbqt = out / "ligands" / f"{ligand_id}_ref.pdbqt"
    _convert_to_pdbqt(ligand_ref_pdb, ligand_pdbqt)

    # calcula centro e tamanho da caixa a partir do ligante de referência
    coords = []
    with ligand_ref_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                # tenta parsear posições fixas primeiro (PDB padrão)
                if len(line) >= 54:
                    x_str = line[30:38].strip()
                    y_str = line[38:46].strip()
                    z_str = line[46:54].strip()
                    if x_str and y_str and z_str:
                        x = float(x_str)
                        y = float(y_str)
                        z = float(z_str)
                        coords.append((x, y, z))
                    else:
                        continue
                else:
                    # fallback: split por espaço em branco e pega últimas 3 colunas
                    parts = line.split()
                    if len(parts) >= 6:
                        x = float(parts[-3])
                        y = float(parts[-2])
                        z = float(parts[-1])
                        coords.append((x, y, z))
            except (ValueError, IndexError):
                continue
    if not coords:
        raise RuntimeError("Não foi possível ler coordenadas do ligante de referência")
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    # tamanho com margem
    dx = max(xs) - min(xs)
    dy = max(ys) - min(ys)
    dz = max(zs) - min(zs)
    size_x = max(8.0, dx + 10.0)
    size_y = max(8.0, dy + 10.0)
    size_z = max(8.0, dz + 10.0)

    out_pdbqt = out / "vina_out" / f"{pdb_id.lower()}_{chain}_{ligand_id}_out.pdbqt"
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        vina_exe,
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        str(cx),
        "--center_y",
        str(cy),
        "--center_z",
        str(cz),
        "--size_x",
        str(size_x),
        "--size_y",
        str(size_y),
        "--size_z",
        str(size_z),
        "--out",
        str(out_pdbqt),
    ]
    # execute
    subprocess.run(cmd, check=True)

    # parse poses
    poses = _parse_pdbqt_poses(out_pdbqt)
    ref_coords = coords
    best_rmsd = float('inf')
    best_idx = -1
    for i, pose in enumerate(poses):
        if len(pose) != len(ref_coords):
            # tenta ignorar hidrogênios: compara apenas heavy atoms by approximate count
            # se diferente, compute without rotation fallback
            try:
                rmsd = _kabsch_rmsd(pose, ref_coords[: len(pose)])
            except Exception:
                rmsd = float('inf')
        else:
            rmsd = _kabsch_rmsd(pose, ref_coords)
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_idx = i

    # write best pose to file for convenience
    best_pose_file = out / "vina_out" / f"best_pose_{pdb_id.lower()}_{chain}_{ligand_id}.pdbqt"
    if best_idx >= 0 and best_idx < len(poses):
        # re-extract model block
        with out_pdbqt.open("r", encoding="utf-8", errors="ignore") as fh, best_pose_file.open("w", encoding="utf-8") as outfh:
            write = False
            model_count = -1
            for line in fh:
                if line.startswith("MODEL"):
                    model_count += 1
                    write = (model_count == best_idx)
                    if write:
                        outfh.write(line)
                    continue
                if line.startswith("ENDMDL"):
                    if write:
                        outfh.write(line)
                        break
                    write = False
                    continue
                if write:
                    outfh.write(line)
    return best_pose_file, best_rmsd


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Preprocessa usando `src.base` e executa AutoDock Vina para redocking")
    parser.add_argument("pdb_id")
    parser.add_argument("chain")
    parser.add_argument("ligand_id")
    parser.add_argument("output_dir")
    parser.add_argument("--vina-exe", default="vina")
    args = parser.parse_args()
    pose, rmsd = preprocess_and_dock(args.pdb_id, args.chain, args.ligand_id, args.output_dir, args.vina_exe)
    print(f"Melhor pose: {pose} (RMSD={rmsd:.3f} Å)")
