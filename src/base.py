# src/base.py: download e preparação de estruturas (proteínas e ligantes)

from __future__ import annotations
import argparse
import shlex
import shutil
import sys
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

WATER_AND_COMMON_SMALL_MOLECULES = {
    "HOH",
    "WAT",
    "DOD",
    "H2O",
    "NA",
    "K",
    "CL",
    "CA",
    "MG",
    "MN",
    "ZN",
    "FE",
    "CU",
    "CO",
    "NI",
    "SO4",
    "PO4",
    "GOL",
    "PEG",
    "EDO",
    "ACT",
    "FMT",
    "BME",
    "MES",
    "MPD",
    "TRS",
    "NO3",
    "NAG",
    "MAN",
    "FUC",
    "BGC",
}
STANDARD_AMINO_ACIDS = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "SEC",
    "PYL",
}


def download_binary(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(request, timeout=60) as response:
        destination.write_bytes(response.read())


def find_local_structure_file(pdb_id: str) -> Path | None:
    workspace_root = Path(__file__).resolve().parents[1]
    search_names = [f"{pdb_id.lower()}.cif", f"{pdb_id.lower()}.pdb"]
    search_roots = [
        workspace_root / "homologia_9thj" / "data",
        workspace_root / "Backup",
        workspace_root / "Documentos",
        workspace_root / "data",
        workspace_root,
    ]
    for root in search_roots:
        if not root.exists():
            continue
        for search_name in search_names:
            candidate = root / search_name
            if candidate.exists():
                return candidate
    for search_name in search_names:
        for candidate in workspace_root.rglob(search_name):
            if candidate.is_file():
                return candidate
    return None


def download_structure(pdb_id: str, output_dir: Path) -> Path:
    pdb_id = pdb_id.lower()
    remote_errors: list[Exception] = []
    for suffix in ("cif", "pdb"):
        destination = output_dir / "raw" / f"{pdb_id}.{suffix}"
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.{suffix}"
        try:
            download_binary(url, destination)
            return destination
        except HTTPError as exc:
            if exc.code != 404:
                remote_errors.append(exc)
        except URLError as exc:
            remote_errors.append(exc)
    local_source = find_local_structure_file(pdb_id)
    if local_source is not None:
        local_destination = output_dir / "raw" / local_source.name
        local_destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(local_source, local_destination)
        return local_destination
    if remote_errors:
        raise remote_errors[-1]
    raise FileNotFoundError(f"Não foi possível localizar a estrutura {pdb_id.upper()}")


def extract_chain_from_pdb(source_file: Path, chain_id: str, output_dir: Path, pdb_id: str) -> Path:
    chain_id = chain_id.strip().upper()
    destination = output_dir / "receptor" / f"{pdb_id.lower()}_{chain_id}.pdb"
    destination.parent.mkdir(parents=True, exist_ok=True)
    filtered_lines: list[str] = []
    with source_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21 and line[21].strip().upper() == chain_id:
                filtered_lines.append(line)
    filtered_lines.append("END\n")
    destination.write_text("".join(filtered_lines), encoding="utf-8")
    return destination


def iter_atom_site_rows(cif_file: Path):
    lines = cif_file.read_text(encoding="utf-8", errors="ignore").splitlines()
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        if line != "loop_":
            index += 1
            continue
        index += 1
        headers: list[str] = []
        while index < len(lines):
            header = lines[index].strip()
            if header.startswith("_atom_site."):
                headers.append(header)
                index += 1
                continue
            break
        if not headers:
            continue
        buffer = ""
        while index < len(lines):
            current = lines[index].rstrip("\n")
            stripped = current.strip()
            if not stripped or stripped.startswith("#"):
                index += 1
                continue
            if stripped == "loop_" or stripped.startswith(("data_", "_")):
                break
            buffer = f"{buffer} {stripped}".strip()
            tokens = shlex.split(buffer, posix=True)
            if len(tokens) < len(headers):
                index += 1
                continue
            if len(tokens) != len(headers):
                buffer = ""
                index += 1
                continue
            yield dict(zip(headers, tokens))
            buffer = ""
            index += 1


def pick_value(row: dict[str, str], *suffixes: str, default: str = "") -> str:
    for suffix in suffixes:
        value = row.get(f"_atom_site.{suffix}", "")
        if value and value not in {"?", "."}:
            return value
    return default


def parse_float(value: str, default: float = 0.0) -> float:
    if not value or value in {"?", "."}:
        return default
    try:
        return float(value)
    except ValueError:
        return default


def format_pdb_atom_name(atom_name: str, element: str) -> str:
    atom_name = atom_name[:4].strip()
    if not atom_name:
        return "    "
    if len(atom_name) >= 4:
        return atom_name[:4]
    if len(element.strip()) == 1:
        return atom_name.rjust(4)
    return atom_name.ljust(4)


def format_pdb_line(record: str, serial: int, atom_name: str, alt_loc: str, res_name: str, chain_id: str, res_seq: str, ins_code: str, x: float, y: float, z: float, occupancy: float, temp_factor: float, element: str,) -> str:
    formatted_atom_name = format_pdb_atom_name(atom_name, element)
    # Build line character by character to ensure correct column positions
    # PDB format (1-based indexing in comments, 0-based in code):
    # 1-6: record name (ATOM  or HETATM)
    # 7-11: atom serial number (right-justified, 5 chars)
    # 12: blank
    # 13-16: atom name (4 chars)
    # 17: alternate location indicator
    # 18-20: residue name (3 chars, right-justified)
    # 21: blank
    # 22: chain identifier
    # 23-26: residue sequence number (4 chars, right-justified)
    # 27: insertion code
    # 28-30: blank (3 spaces)
    # 31-38: X coordinate (8 chars, right-justified)
    # 39-46: Y coordinate (8 chars, right-justified)
    # 47-54: Z coordinate (8 chars, right-justified)
    # 55-60: occupancy (6 chars, right-justified)
    # 61-66: temperature factor (6 chars, right-justified)
    # 67-76: blank (10 spaces)
    # 77-78: element symbol (2 chars, right-justified)
    # 79-80: charge
    parts = [
        f"{record:<6}",           # 1-6
        f"{serial:>5d}",          # 7-11
        " ",                      # 12
        f"{formatted_atom_name:<4}",  # 13-16
        f"{alt_loc:1}",           # 17
        f"{res_name:>3}",         # 18-20
        " ",                      # 21
        f"{chain_id:1}",          # 22
        f"{res_seq:>4}",          # 23-26
        f"{ins_code:1}",          # 27
        "   ",                    # 28-30
        f"{x:>8.3f}",             # 31-38
        f"{y:>8.3f}",             # 39-46
        f"{z:>8.3f}",             # 47-54
        f"{occupancy:>6.2f}",     # 55-60
        f"{temp_factor:>6.2f}",   # 61-66
        "          ",             # 67-76
        f"{element[:2]:>2s}",     # 77-78
        "  \n",                   # 79-82
    ]
    return "".join(parts)


def extract_chain_from_mmcif(source_file: Path, chain_id: str, output_dir: Path, pdb_id: str) -> Path:
    chain_id = chain_id.strip().upper()
    destination = output_dir / "receptor" / f"{pdb_id.lower()}_{chain_id}.pdb"
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    serial = 1
    selected_model: str | None = None
    for row in iter_atom_site_rows(source_file):
        group_pdb = pick_value(row, "group_PDB").upper()
        if group_pdb != "ATOM":
            continue
        model_number = pick_value(row, "pdbx_PDB_model_num", default="1")
        if selected_model is None:
            selected_model = model_number
        if model_number != selected_model:
            continue
        row_chain_id = pick_value(row, "auth_asym_id", "label_asym_id").upper()
        if row_chain_id != chain_id:
            continue
        res_name = pick_value(row, "auth_comp_id", "label_comp_id")
        atom_name = pick_value(row, "auth_atom_id", "label_atom_id", "id")
        alt_loc = pick_value(row, "label_alt_id")
        if alt_loc in {"?", "."}:
            alt_loc = ""
        res_seq = pick_value(row, "auth_seq_id", "label_seq_id")
        ins_code = pick_value(row, "pdbx_PDB_ins_code")
        if ins_code in {"?", "."}:
            ins_code = ""
        element = pick_value(row, "type_symbol", default=(atom_name[:1] if atom_name else ""))
        x = parse_float(pick_value(row, "Cartn_x"))
        y = parse_float(pick_value(row, "Cartn_y"))
        z = parse_float(pick_value(row, "Cartn_z"))
        occupancy = parse_float(pick_value(row, "occupancy"), default=1.0)
        temp_factor = parse_float(pick_value(row, "B_iso_or_equiv"))
        lines.append(
            format_pdb_line(
                "ATOM",
                serial,
                atom_name,
                alt_loc,
                res_name,
                chain_id[:1] if chain_id else " ",
                res_seq,
                ins_code,
                x,
                y,
                z,
                occupancy,
                temp_factor,
                element,
            )
        )
        serial += 1
    lines.append("END\n")
    destination.write_text("".join(lines), encoding="utf-8")
    return destination


def detect_ligands(source_file: Path) -> list[str]:
    ligands: set[str] = set()
    with source_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("HETATM"):
                continue
            residue_name = line[17:20].strip().upper()
            if not residue_name:
                continue
            if residue_name in WATER_AND_COMMON_SMALL_MOLECULES:
                continue
            if residue_name in STANDARD_AMINO_ACIDS:
                continue
            ligands.add(residue_name)
    return sorted(ligands)


def detect_ligands_from_mmcif(source_file: Path) -> list[str]:
    ligands: set[str] = set()
    for row in iter_atom_site_rows(source_file):
        if pick_value(row, "group_PDB").upper() != "HETATM":
            continue
        residue_name = pick_value(row, "auth_comp_id", "label_comp_id").upper()
        if not residue_name:
            continue
        if residue_name in WATER_AND_COMMON_SMALL_MOLECULES:
            continue
        if residue_name in STANDARD_AMINO_ACIDS:
            continue
        ligands.add(residue_name)
    return sorted(ligands)


def download_ligand(ligand_id: str, output_dir: Path) -> Path:
    ligand_id = ligand_id.upper()
    destination = output_dir / "ligands" / f"{ligand_id}.sdf"
    url = f"https://files.rcsb.org/ligands/download/{ligand_id}_ideal.sdf"
    download_binary(url, destination)
    return destination


def parse_ligand_ids(raw_ids: Iterable[str] | None) -> list[str]:
    if not raw_ids:
        return []
    parsed: list[str] = []
    for item in raw_ids:
        for ligand_id in item.split(","):
            ligand_id = ligand_id.strip().upper()
            if ligand_id:
                parsed.append(ligand_id)
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Baixa a estrutura da proteína e os ligantes do PDB.")
    parser.add_argument("--pdb-id", required=True, help="Identificador PDB, por exemplo 9THJ")
    parser.add_argument("--chain", required=True, help="Cadeia da proteína a extrair, por exemplo A")
    parser.add_argument("--output-dir", required=True, help="Diretório de saída")
    parser.add_argument(
        "--ligand-id",
        action="append",
        dest="ligand_ids",
        help="Código do ligante no PDB. Pode ser repetido ou separado por vírgula.",
    )
    parser.add_argument(
        "--no-auto-ligands",
        action="store_true",
        help="Desativa a detecção automática de ligantes a partir da estrutura baixada.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    pdb_id = args.pdb_id.strip().upper()
    chain_id = args.chain.strip().upper()
    output_dir = Path(args.output_dir)
    try:
        print(f"Baixando estrutura {pdb_id}...")
        raw_structure = download_structure(pdb_id, output_dir)
        print(f"Extraindo cadeia {chain_id}...")
        if raw_structure.suffix.lower() == ".cif":
            receptor_file = extract_chain_from_mmcif(raw_structure, chain_id, output_dir, pdb_id)
            ligand_detector = detect_ligands_from_mmcif
        else:
            receptor_file = extract_chain_from_pdb(raw_structure, chain_id, output_dir, pdb_id)
            ligand_detector = detect_ligands
        print(f"Proteína salva em: {receptor_file}")
        ligand_ids = parse_ligand_ids(args.ligand_ids)
        if not ligand_ids and not args.no_auto_ligands:
            ligand_ids = ligand_detector(raw_structure)
        if ligand_ids:
            print("Baixando ligantes: " + ", ".join(ligand_ids))
            for ligand_id in ligand_ids:
                try:
                    ligand_file = download_ligand(ligand_id, output_dir)
                    print(f"Ligante salvo em: {ligand_file}")
                except (HTTPError, URLError, TimeoutError) as exc:
                    print(f"Aviso: não foi possível baixar {ligand_id}: {exc}", file=sys.stderr)
        else:
            print("Nenhum ligante detectado/solicitado.")
    except (HTTPError, URLError, TimeoutError, OSError) as exc:
        print(f"Erro: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
