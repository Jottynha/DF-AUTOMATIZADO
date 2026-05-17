
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict
from pathlib import Path

try:
    import base as base_module
except ImportError:
    from .. import base as base_module  # type: ignore


@dataclass(frozen=True)
class Dock6PreparedInputs:
    output_dir: str
    raw_structure: str
    receptor_clean_pdb: str
    receptor_pqr: str
    receptor_fixed_pdb: str
    receptor_fixed_mol2: str
    ligand_ref_pdb: str
    ligand_ref_mol2: str
    ligand_dock_pdb: str
    ligand_dock_mol2: str
    receptor_ms: str


def ensure_executable(command: str) -> None:
    if shutil.which(command) is None:
        raise FileNotFoundError(f"Não encontrei o executável `{command}` no PATH")


def run_command(command: list[str], cwd: Path | None = None, log_file: Path | None = None) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    log_text = "\n".join(part for part in [completed.stdout.strip(), completed.stderr.strip()] if part)

    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        log_file.write_text(log_text, encoding="utf-8")

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

    return completed


def convert_file(
    input_file: Path,
    output_file: Path,
    obabel_exe: str,
    input_format: str | None = None,
    output_format: str | None = None,
    extra_args: list[str] | None = None,
) -> None:
    if not input_file.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {input_file}")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    in_format = input_format or input_file.suffix.lstrip(".").lower()
    out_format = output_format or output_file.suffix.lstrip(".").lower()

    command = [obabel_exe, f"-i{in_format}", str(input_file), f"-o{out_format}", "-O", str(output_file)]

    if extra_args:
        command.extend(extra_args)

    run_command(command, log_file=output_file.with_suffix(output_file.suffix + ".log"))


def write_clean_receptor_pdb(raw_structure: Path, chain_id: str, output_file: Path, pdb_id: str) -> Path:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if raw_structure.suffix.lower() == ".cif":
        temp_file = base_module.extract_chain_from_mmcif(raw_structure, chain_id, output_file.parent, pdb_id)
    else:
        temp_file = base_module.extract_chain_from_pdb(raw_structure, chain_id, output_file.parent, pdb_id)

    clean_lines: list[str] = []
    with temp_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM  "):
                continue
            element = line[76:78].strip().upper()
            atom_name = line[12:16].strip().upper()
            if element == "H" or atom_name.startswith("H"):
                continue
            clean_lines.append(line)

    clean_lines.append("END\n")
    output_file.write_text("".join(clean_lines), encoding="utf-8")

    if temp_file != output_file and temp_file.exists():
        try:
            temp_file.unlink()
        except OSError:
            pass

    return output_file


def extract_ligand_ref_from_pdb(source_file: Path, ligand_id: str, output_file: Path) -> Path:
    ligand_id_upper = ligand_id.upper()
    output_file.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    with source_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("HETATM"):
                continue
            residue_name = line[17:20].strip().upper()
            if residue_name != ligand_id_upper[:3]:
                continue
            lines.append(line)

    if not lines:
        raise ValueError(f"Não encontrei o ligante {ligand_id} no PDB {source_file}")

    lines.append("END\n")
    output_file.write_text("".join(lines), encoding="utf-8")
    return output_file


def extract_ligand_ref_from_mmcif(source_file: Path, ligand_id: str, output_file: Path) -> Path:
    ligand_id_upper = ligand_id.upper()
    output_file.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    serial = 1
    selected_model: str | None = None

    for row in base_module.iter_atom_site_rows(source_file):
        group_pdb = base_module.pick_value(row, "group_PDB").upper()
        if group_pdb != "HETATM":
            continue

        comp_id = base_module.pick_value(row, "auth_comp_id", "label_comp_id").upper()
        if comp_id != ligand_id_upper:
            continue

        model_number = base_module.pick_value(row, "pdbx_PDB_model_num", default="1")
        if selected_model is None:
            selected_model = model_number
        if model_number != selected_model:
            continue

        atom_name = base_module.pick_value(row, "auth_atom_id", "label_atom_id", "id")
        element = base_module.pick_value(row, "type_symbol", default=(atom_name[:1] if atom_name else ""))
        chain = base_module.pick_value(row, "auth_asym_id", "label_asym_id", default="L")[:1]
        res_seq = base_module.pick_value(row, "auth_seq_id", "label_seq_id", default="1")
        x = base_module.parse_float(base_module.pick_value(row, "Cartn_x"))
        y = base_module.parse_float(base_module.pick_value(row, "Cartn_y"))
        z = base_module.parse_float(base_module.pick_value(row, "Cartn_z"))
        occupancy = base_module.parse_float(base_module.pick_value(row, "occupancy"), default=1.0)
        temp_factor = base_module.parse_float(base_module.pick_value(row, "B_iso_or_equiv"))

        lines.append(
            base_module.format_pdb_line(
                "HETATM",
                serial,
                atom_name,
                "",
                comp_id[:3],
                chain or "L",
                res_seq,
                "",
                x,
                y,
                z,
                occupancy,
                temp_factor,
                element,
            )
        )
        serial += 1

    if not lines:
        raise ValueError(f"Não encontrei o ligante {ligand_id} no mmCIF {source_file}")

    lines.append("END\n")
    output_file.write_text("".join(lines), encoding="utf-8")
    return output_file


def extract_ligand_ref(raw_structure: Path, ligand_id: str, output_file: Path) -> Path:
    if raw_structure.suffix.lower() == ".cif":
        return extract_ligand_ref_from_mmcif(raw_structure, ligand_id, output_file)
    return extract_ligand_ref_from_pdb(raw_structure, ligand_id, output_file)


def run_pdb2pqr(receptor_clean_pdb: Path, receptor_pqr: Path, pdb2pqr_exe: str, ph: float) -> None:
    receptor_pqr.parent.mkdir(parents=True, exist_ok=True)

    command = [
        pdb2pqr_exe,
        "--ff=AMBER",
        f"--with-ph={ph}",
        str(receptor_clean_pdb),
        str(receptor_pqr),
    ]

    run_command(command, log_file=receptor_pqr.with_suffix(".pdb2pqr.log"))


def run_dms(receptor_clean_pdb: Path, receptor_ms: Path, dms_exe: str) -> None:
    receptor_ms.parent.mkdir(parents=True, exist_ok=True)

    command = [
        dms_exe,
        str(receptor_clean_pdb),
        "-n",
        "-w",
        "1.4",
        "-v",
        "-o",
        str(receptor_ms),
    ]

    run_command(command, log_file=receptor_ms.with_suffix(".dms.log"))


def parse_pqr_atoms(pqr_file: Path) -> list[dict[str, object]]:
    """
    Lê átomos de um PQR gerado pelo pdb2pqr.

    O PQR é parecido com PDB, mas inclui carga e raio no final.
    Tentamos primeiro por split, porque o pdb2pqr pode variar o espaçamento.
    """
    atoms: list[dict[str, object]] = []

    with pqr_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            parts = line.split()
            if len(parts) < 10:
                continue

            try:
                atom_name = parts[2]
                res_name = parts[3]

                # Formatos comuns:
                # ATOM serial atom res chain resseq x y z charge radius
                # ATOM serial atom res resseq x y z charge radius
                if len(parts) >= 11 and not _looks_like_number(parts[4]):
                    chain_id = parts[4]
                    res_seq = parts[5]
                else:
                    chain_id = ""
                    res_seq = parts[4]

                charge = float(parts[-2])
            except (ValueError, IndexError):
                continue

            atoms.append(
                {
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "chain_id": chain_id,
                    "res_seq": res_seq,
                    "charge": charge,
                }
            )

    return atoms


def _looks_like_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def normalize_receptor_mol2_for_dock6(pqr_file: Path, mol2_file: Path) -> None:
    """
    Corrige o receptor_fixed.mol2 gerado pelo Open Babel para um formato
    mais compatível com o grid do DOCK6.

    O Open Babel tende a escrever:
      SMALL
      GASTEIGER
      0 substructures
      cargas zeradas ou inadequadas

    O grid do DOCK6, no fluxo da disciplina, espera um receptor em padrão
    de proteína, parecido com:
      PROTEIN
      AMBER ff14SB
      subestruturas/resíduos definidos

    Esta função:
    - troca o nome da molécula para receptor_fixed.pdb;
    - muda SMALL/GASTEIGER para PROTEIN/AMBER ff14SB;
    - usa cargas do PQR gerado pelo pdb2pqr;
    - cria subestruturas por resíduo;
    - reescreve a seção @<TRIPOS>ATOM.
    """
    pqr_atoms = parse_pqr_atoms(pqr_file)
    if not pqr_atoms:
        raise ValueError(f"Não consegui ler átomos/cargas do PQR: {pqr_file}")

    lines = mol2_file.read_text(encoding="utf-8", errors="ignore").splitlines()

    try:
        atom_start = lines.index("@<TRIPOS>ATOM")
    except ValueError as exc:
        raise ValueError(f"MOL2 sem seção @<TRIPOS>ATOM: {mol2_file}") from exc

    next_section = len(lines)
    for i in range(atom_start + 1, len(lines)):
        if lines[i].startswith("@<TRIPOS>"):
            next_section = i
            break

    atom_lines = [line for line in lines[atom_start + 1:next_section] if line.strip()]
    if len(atom_lines) != len(pqr_atoms):
        raise ValueError(
            f"Número de átomos diferente entre PQR ({len(pqr_atoms)}) "
            f"e MOL2 ({len(atom_lines)})."
        )

    # Preserva a seção BOND e ignora a SUBSTRUCTURE do Open Babel, se houver.
    bond_lines: list[str] = []
    if "@<TRIPOS>BOND" in lines:
        bond_start = lines.index("@<TRIPOS>BOND")
        bond_end = len(lines)
        for i in range(bond_start + 1, len(lines)):
            if lines[i].startswith("@<TRIPOS>"):
                bond_end = i
                break
        bond_lines = [line for line in lines[bond_start + 1:bond_end] if line.strip()]

    residue_map: dict[tuple[str, str, str], int] = {}
    residue_names: dict[int, str] = {}
    new_atom_lines: list[str] = []

    for idx, (mol2_atom_line, pqr_atom) in enumerate(zip(atom_lines, pqr_atoms), start=1):
        parts = mol2_atom_line.split()
        if len(parts) < 6:
            raise ValueError(f"Linha MOL2 inválida: {mol2_atom_line}")

        atom_id = int(parts[0])
        atom_name = str(pqr_atom["atom_name"])
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])
        atom_type = parts[5]

        res_name = str(pqr_atom["res_name"])
        res_seq = str(pqr_atom["res_seq"])
        chain_id = str(pqr_atom["chain_id"])
        charge = float(pqr_atom["charge"])

        residue_key = (chain_id, res_seq, res_name)
        if residue_key not in residue_map:
            residue_map[residue_key] = len(residue_map) + 1
            residue_names[residue_map[residue_key]] = res_name

        subst_id = residue_map[residue_key]

        new_atom_lines.append(
            f"{atom_id:7d} {atom_name:<8} {x:10.4f} {y:10.4f} {z:10.4f} "
            f"{atom_type:<6} {subst_id:5d} {res_name:<8} {charge:10.4f}"
        )

    new_substructure_lines: list[str] = []
    for subst_id in sorted(residue_names):
        res_name = residue_names[subst_id]
        new_substructure_lines.append(
            f"{subst_id:6d} {res_name:<8} {1:5d} RESIDUE           4 ****  ****    0 ROOT"
        )

    atom_count = len(new_atom_lines)
    bond_count = len(bond_lines)
    subst_count = len(new_substructure_lines)

    new_lines = [
        "@<TRIPOS>MOLECULE",
        "receptor_fixed.pdb",
        f"{atom_count} {bond_count} {subst_count} 0 0",
        "PROTEIN",
        "AMBER ff14SB",
        "",
        "",
        "@<TRIPOS>ATOM",
        *new_atom_lines,
    ]

    if bond_lines:
        new_lines.extend(["@<TRIPOS>BOND", *bond_lines])

    new_lines.extend(["@<TRIPOS>SUBSTRUCTURE", *new_substructure_lines, ""])

    mol2_file.write_text("\n".join(new_lines), encoding="utf-8")


def prepare_dock6_inputs(
    pdb_id: str,
    chain_id: str,
    ligand_id: str,
    output_dir: Path,
    pdb2pqr_exe: str = "pdb2pqr",
    obabel_exe: str = "obabel",
    dms_exe: str = "dms",
    ph: float = 7.4,
) -> Dock6PreparedInputs:
    ensure_executable(pdb2pqr_exe)
    ensure_executable(obabel_exe)
    ensure_executable(dms_exe)

    output_dir.mkdir(parents=True, exist_ok=True)

    raw_structure = base_module.download_structure(pdb_id, output_dir)

    receptor_clean_pdb = output_dir / "receptor_clean.pdb"
    receptor_pqr = output_dir / "receptor_pH74.pqr"
    receptor_fixed_pdb = output_dir / "receptor_fixed.pdb"
    receptor_fixed_mol2 = output_dir / "receptor_fixed.mol2"

    ligand_ref_pdb = output_dir / "ligand_ref.pdb"
    ligand_ref_mol2 = output_dir / "ligand_ref.mol2"
    ligand_dock_pdb = output_dir / "ligand_dock.pdb"
    ligand_dock_mol2 = output_dir / "ligand_dock.mol2"

    receptor_ms = output_dir / "rec.ms"

    print("[prepare] Gerando receptor_clean.pdb")
    write_clean_receptor_pdb(raw_structure, chain_id, receptor_clean_pdb, pdb_id)

    print("[prepare] Extraindo ligante cristalográfico para ligand_ref.pdb")
    extract_ligand_ref(raw_structure, ligand_id, ligand_ref_pdb)

    print("[prepare] Rodando pdb2pqr para gerar receptor_pH74.pqr")
    run_pdb2pqr(receptor_clean_pdb, receptor_pqr, pdb2pqr_exe, ph)

    print("[prepare] Convertendo receptor_pH74.pqr para receptor_fixed.pdb")
    convert_file(receptor_pqr, receptor_fixed_pdb, obabel_exe, input_format="pqr", output_format="pdb")

    print("[prepare] Convertendo receptor_pH74.pqr para receptor_fixed.mol2")
    convert_file(receptor_pqr, receptor_fixed_mol2, obabel_exe, input_format="pqr", output_format="mol2")

    print("[prepare] Normalizando receptor_fixed.mol2 para formato compatível com DOCK6")
    normalize_receptor_mol2_for_dock6(receptor_pqr, receptor_fixed_mol2)

    print("[prepare] Convertendo ligand_ref.pdb para ligand_ref.mol2")
    convert_file(ligand_ref_pdb, ligand_ref_mol2, obabel_exe, input_format="pdb", output_format="mol2")

    print("[prepare] Gerando ligand_dock.pdb com hidrogênios")
    convert_file(ligand_ref_pdb, ligand_dock_pdb, obabel_exe, input_format="pdb", output_format="pdb", extra_args=["-h"])

    print("[prepare] Convertendo ligand_dock.pdb para ligand_dock.mol2")
    convert_file(
        ligand_dock_pdb,
        ligand_dock_mol2,
        obabel_exe,
        input_format="pdb",
        output_format="mol2",
        extra_args=["--gen3d", "-p", str(ph), "--partialcharge", "gasteiger"],
    )

    print("[prepare] Gerando rec.ms com dms")
    run_dms(receptor_clean_pdb, receptor_ms, dms_exe)

    result = Dock6PreparedInputs(
        output_dir=str(output_dir),
        raw_structure=str(raw_structure),
        receptor_clean_pdb=str(receptor_clean_pdb),
        receptor_pqr=str(receptor_pqr),
        receptor_fixed_pdb=str(receptor_fixed_pdb),
        receptor_fixed_mol2=str(receptor_fixed_mol2),
        ligand_ref_pdb=str(ligand_ref_pdb),
        ligand_ref_mol2=str(ligand_ref_mol2),
        ligand_dock_pdb=str(ligand_dock_pdb),
        ligand_dock_mol2=str(ligand_dock_mol2),
        receptor_ms=str(receptor_ms),
    )

    (output_dir / "prepared_inputs.json").write_text(json.dumps(asdict(result), indent=2), encoding="utf-8")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepara entradas para o pipeline DOCK6.")
    parser.add_argument("pdb_id")
    parser.add_argument("chain")
    parser.add_argument("ligand_id")
    parser.add_argument("output_dir")
    parser.add_argument("--pdb2pqr-exe", default="pdb2pqr")
    parser.add_argument("--obabel-exe", default="obabel")
    parser.add_argument("--dms-exe", default="dms")
    parser.add_argument("--ph", type=float, default=7.4)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        result = prepare_dock6_inputs(
            pdb_id=args.pdb_id,
            chain_id=args.chain,
            ligand_id=args.ligand_id,
            output_dir=Path(args.output_dir),
            pdb2pqr_exe=args.pdb2pqr_exe,
            obabel_exe=args.obabel_exe,
            dms_exe=args.dms_exe,
            ph=args.ph,
        )
        print()
        print("Entradas DOCK6 preparadas:")
        print(json.dumps(asdict(result), indent=2))
        return 0
    except Exception as exc:
        print(f"Erro: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())