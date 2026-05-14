# src/analysis.py: análise de docking (RMSD e scores)

from __future__ import annotations
import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

@dataclass(frozen=True)
class AtomCoords:
    x: float
    y: float
    z: float

def parse_pdb_atoms(pdb_file: Path) -> list[tuple[str, AtomCoords]]:
    atoms: list[tuple[str, AtomCoords]] = []
    with pdb_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((atom_name, AtomCoords(x, y, z)))
            except (ValueError, IndexError):
                continue
    return atoms


def extract_score_from_vina_log(log_file: Path) -> float | None:
    with log_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            match = re.match(r"\s*1\s+(-\d+\.\d+)", line)
            if match:
                return float(match.group(1))
    return None


def kabsch_align(atoms1: list[tuple[str, AtomCoords]], atoms2: list[tuple[str, AtomCoords]]):
    try:
        import numpy as np
    except ImportError:
        return None
    if not atoms1 or not atoms2:
        return None
    coords1 = np.array([[c.x, c.y, c.z] for _, c in atoms1], dtype=float)
    coords2 = np.array([[c.x, c.y, c.z] for _, c in atoms2], dtype=float)
    if coords1.shape[0] != coords2.shape[0]:
        min_n = min(coords1.shape[0], coords2.shape[0])
        coords1 = coords1[:min_n]
        coords2 = coords2[:min_n]
    if coords1.shape[0] < 3:
        return None
    centroid1 = coords1.mean(axis=0)
    centroid2 = coords2.mean(axis=0)
    coords1_centered = coords1 - centroid1
    coords2_centered = coords2 - centroid2
    cov_matrix = np.dot(coords1_centered.T, coords2_centered)
    U, _, VT = np.linalg.svd(cov_matrix)
    rotation = np.dot(VT.T, U.T)
    if np.linalg.det(rotation) < 0:
        VT[-1, :] *= -1
        rotation = np.dot(VT.T, U.T)
    coords1_aligned = np.dot(coords1_centered, rotation.T) + centroid2
    aligned1 = [AtomCoords(float(c[0]), float(c[1]), float(c[2])) for c in coords1_aligned]
    coords2_list = [AtomCoords(float(c[0]), float(c[1]), float(c[2])) for c in coords2]
    return (aligned1, coords2_list)


def simple_rmsd(atoms1: list[tuple[str, AtomCoords]], atoms2: list[tuple[str, AtomCoords]]) -> float | None:
    if not atoms1 or not atoms2:
        return None
    result = kabsch_align(atoms1, atoms2)
    if result is not None:
        aligned1, coords2 = result
        sum_sq = 0.0
        for coord1, coord2 in zip(aligned1, coords2):
            dist_sq = ((coord1.x - coord2.x) ** 2 + (coord1.y - coord2.y) ** 2 + (coord1.z - coord2.z) ** 2)
            sum_sq += dist_sq
        return math.sqrt(sum_sq / len(aligned1))
    sorted_atoms1 = sorted(atoms1, key=lambda x: x[0])
    sorted_atoms2 = sorted(atoms2, key=lambda x: x[0])
    min_count = min(len(sorted_atoms1), len(sorted_atoms2))
    if min_count < 3:
        return None
    sum_sq = 0.0
    for i in range(min_count):
        _, coord1 = sorted_atoms1[i]
        _, coord2 = sorted_atoms2[i]
        dist_sq = ((coord1.x - coord2.x) ** 2 + (coord1.y - coord2.y) ** 2 + (coord1.z - coord2.z) ** 2)
        sum_sq += dist_sq
    return math.sqrt(sum_sq / min_count)


def analyze_docking_result(docking_dir: Path, ligand_crystal_file: Path | None = None) -> dict:
    result = {"docking_dir": str(docking_dir), "best_score": None, "rmsd": None, "num_atoms_docked": None, "num_atoms_crystal": None}
    log_file = docking_dir / "log.txt"
    if log_file.exists():
        result["best_score"] = extract_score_from_vina_log(log_file)
    if ligand_crystal_file and ligand_crystal_file.exists():
        out_pdbqt = docking_dir / "out.pdbqt"
        if out_pdbqt.exists():
            try:
                atoms_crystal = parse_pdb_atoms(ligand_crystal_file)
                atoms_docked = parse_pdb_atoms(out_pdbqt)
                result["num_atoms_crystal"] = len(atoms_crystal)
                result["num_atoms_docked"] = len(atoms_docked)
                result["rmsd"] = simple_rmsd(atoms_crystal, atoms_docked)
            except Exception as exc:
                result["error"] = str(exc)
    return result


def analyze_multiple_dockings(docking_base_dir: Path, ligand_crystal_dir: Path | None = None) -> list[dict]:
    results = []
    if not docking_base_dir.exists():
        return results
    for ligand_dir in sorted(docking_base_dir.iterdir()):
        if not ligand_dir.is_dir():
            continue
        ligand_name = ligand_dir.name
        crystal_file = None
        if ligand_crystal_dir:
            for suffix in (".pdb", ".sdf", ".mol2"):
                candidate = ligand_crystal_dir / f"{ligand_name}{suffix}"
                if candidate.exists():
                    crystal_file = candidate
                    break
        result = analyze_docking_result(ligand_dir, crystal_file)
        result["ligand"] = ligand_name
        results.append(result)
    return results


def format_results_table(results: list[dict]) -> str:
    if not results:
        return "Nenhum resultado para analisar."
    lines = []
    lines.append("=" * 100)
    lines.append(f"{'Ligante':<20} {'Score (kcal/mol)':<20} {'RMSD (Å)':<20} {'Átomos (Docked/Crystal)':<20}")
    lines.append("-" * 100)
    for result in results:
        ligand = result.get("ligand", "N/A")
        score = result.get("best_score")
        rmsd = result.get("rmsd")
        n_docked = result.get("num_atoms_docked")
        n_crystal = result.get("num_atoms_crystal")
        score_str = f"{score:.2f}" if score is not None else "N/A"
        rmsd_str = f"{rmsd:.2f}" if rmsd is not None else "N/A"
        atoms_str = f"{n_docked}/{n_crystal}" if n_docked is not None else "N/A"
        lines.append(f"{ligand:<20} {score_str:<20} {rmsd_str:<20} {atoms_str:<20}")
    lines.append("=" * 100)
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--docking-dir")
    parser.add_argument("--ligand-crystal")
    parser.add_argument("--ligand-crystal-dir")
    parser.add_argument("--output-json")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if not args.docking_dir:
            print("Erro: --docking-dir é obrigatório", file=__import__("sys").stderr)
            return 1
        docking_dir = Path(args.docking_dir)
        ligand_crystal = Path(args.ligand_crystal) if args.ligand_crystal else None
        ligand_crystal_dir = Path(args.ligand_crystal_dir) if args.ligand_crystal_dir else None
        if (docking_dir / "out.pdbqt").exists():
            result = analyze_docking_result(docking_dir, ligand_crystal)
            print(json.dumps(result, indent=2))
            if args.output_json:
                Path(args.output_json).write_text(json.dumps([result], indent=2))
        else:
            results = analyze_multiple_dockings(docking_dir, ligand_crystal_dir)
            print(format_results_table(results))
            if args.output_json:
                Path(args.output_json).write_text(json.dumps(results, indent=2))
        return 0
    except Exception as exc:
        print(f"Erro: {exc}", file=__import__("sys").stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
