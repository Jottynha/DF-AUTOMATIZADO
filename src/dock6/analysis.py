from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class AtomCoords:
    atom_name: str
    x: float
    y: float
    z: float


def parse_mol2_atoms(mol2_file: Path) -> list[AtomCoords]:
    atoms: list[AtomCoords] = []
    in_atom_section = False

    with mol2_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()

            if stripped.startswith("@<TRIPOS>ATOM"):
                in_atom_section = True
                continue

            if stripped.startswith("@<TRIPOS>") and in_atom_section:
                break

            if not in_atom_section or not stripped:
                continue

            parts = stripped.split()
            if len(parts) < 6:
                continue

            try:
                atom_name = parts[1]
                x = float(parts[2])
                y = float(parts[3])
                z = float(parts[4])
            except ValueError:
                continue

            atoms.append(AtomCoords(atom_name=atom_name, x=x, y=y, z=z))

    return atoms


def extract_grid_score_from_mol2(mol2_file: Path) -> float | None:
    if not mol2_file.exists():
        return None

    patterns = [
        r"Grid_Score:\s*([-+]?\d+(?:\.\d+)?)",
        r"Grid Score:\s*([-+]?\d+(?:\.\d+)?)",
        r"##########\s*Grid_Score:\s*([-+]?\d+(?:\.\d+)?)",
    ]

    text = mol2_file.read_text(encoding="utf-8", errors="ignore")

    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return float(match.group(1))

    return None


def extract_grid_score_from_dock_out(dock_out: Path) -> float | None:
    if not dock_out.exists():
        return None

    text = dock_out.read_text(encoding="utf-8", errors="ignore")

    patterns = [
        r"Grid_Score:\s*([-+]?\d+(?:\.\d+)?)",
        r"Grid Score:\s*([-+]?\d+(?:\.\d+)?)",
    ]

    values: list[float] = []

    for pattern in patterns:
        for match in re.finditer(pattern, text):
            values.append(float(match.group(1)))

    if not values:
        return None

    return min(values)


def kabsch_rmsd(atoms1: list[AtomCoords], atoms2: list[AtomCoords]) -> float | None:
    try:
        import numpy as np
    except ImportError:
        return simple_rmsd_without_alignment(atoms1, atoms2)

    if not atoms1 or not atoms2:
        return None

    n = min(len(atoms1), len(atoms2))
    if n < 3:
        return None

    coords1 = np.array([[a.x, a.y, a.z] for a in atoms1[:n]], dtype=float)
    coords2 = np.array([[a.x, a.y, a.z] for a in atoms2[:n]], dtype=float)

    centroid1 = coords1.mean(axis=0)
    centroid2 = coords2.mean(axis=0)

    coords1_centered = coords1 - centroid1
    coords2_centered = coords2 - centroid2

    covariance = coords1_centered.T @ coords2_centered

    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T

    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T

    coords1_aligned = coords1_centered @ rotation.T + centroid2

    diff = coords1_aligned - coords2
    rmsd = math.sqrt(float((diff * diff).sum()) / n)

    return rmsd


def simple_rmsd_without_alignment(
    atoms1: list[AtomCoords],
    atoms2: list[AtomCoords],
) -> float | None:
    if not atoms1 or not atoms2:
        return None

    n = min(len(atoms1), len(atoms2))
    if n < 3:
        return None

    sum_sq = 0.0

    for atom1, atom2 in zip(atoms1[:n], atoms2[:n]):
        sum_sq += (
            (atom1.x - atom2.x) ** 2
            + (atom1.y - atom2.y) ** 2
            + (atom1.z - atom2.z) ** 2
        )

    return math.sqrt(sum_sq / n)


def analyze_dock6_result(
    dock6_run_dir: Path,
    crystal_ligand_mol2: Path | None = None,
) -> dict:
    scored_mol2 = dock6_run_dir / "scored.mol2"
    dock6_run_scored = dock6_run_dir / "dock6_run_scored.mol2"
    dock_out = dock6_run_dir / "dock6.out"

    if not scored_mol2.exists() and dock6_run_scored.exists():
        scored_mol2 = dock6_run_scored

    result = {
        "dock6_run_dir": str(dock6_run_dir),
        "scored_mol2": str(scored_mol2) if scored_mol2.exists() else None,
        "dock_out": str(dock_out) if dock_out.exists() else None,
        "grid_score": None,
        "rmsd": None,
        "num_atoms_docked": None,
        "num_atoms_crystal": None,
    }

    if scored_mol2.exists():
        result["grid_score"] = extract_grid_score_from_mol2(scored_mol2)

    if result["grid_score"] is None:
        result["grid_score"] = extract_grid_score_from_dock_out(dock_out)

    if scored_mol2.exists():
        atoms_docked = parse_mol2_atoms(scored_mol2)
        result["num_atoms_docked"] = len(atoms_docked)
    else:
        atoms_docked = []

    if crystal_ligand_mol2 and crystal_ligand_mol2.exists():
        atoms_crystal = parse_mol2_atoms(crystal_ligand_mol2)
        result["num_atoms_crystal"] = len(atoms_crystal)

        if atoms_docked and atoms_crystal:
            result["rmsd"] = kabsch_rmsd(atoms_crystal, atoms_docked)

    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Analisa resultado de docking DOCK6.")
    parser.add_argument("--dock6-run-dir", required=True)
    parser.add_argument("--crystal-ligand-mol2")
    parser.add_argument("--output-json")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    dock6_run_dir = Path(args.dock6_run_dir)
    crystal_ligand_mol2 = (
        Path(args.crystal_ligand_mol2) if args.crystal_ligand_mol2 else None
    )

    result = analyze_dock6_result(dock6_run_dir, crystal_ligand_mol2)

    print(json.dumps(result, indent=2))

    if args.output_json:
        Path(args.output_json).write_text(json.dumps(result, indent=2), encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())