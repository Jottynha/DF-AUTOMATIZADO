# Pipeline de docking apoiado em `src/base.py`. Exemplo:
# python3 src/docking.py --pdb-id 9THJ --chain A --output-dir resultados/9thj

from __future__ import annotations
import argparse
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
try:
	import base as base_module
except ImportError:  # pragma: no cover - fallback para execução como módulo
	from . import base as base_module  # type: ignore

@dataclass(frozen=True)
class Box:
	center_x: float
	center_y: float
	center_z: float
	size_x: float
	size_y: float
	size_z: float

def parse_float_or_none(value: str | None) -> float | None:
	if value is None:
		return None
	return float(value)

def ensure_executable(command: str) -> None:
	if shutil.which(command) is None:
		raise FileNotFoundError(f"Não encontrei o executável `{command}` no PATH")

def extract_pdb_coordinates(pdb_file: Path) -> list[tuple[float, float, float]]:
	coordinates: list[tuple[float, float, float]] = []
	with pdb_file.open("r", encoding="utf-8", errors="ignore") as handle:
		for line in handle:
			if not line.startswith(("ATOM  ", "HETATM")):
				continue
			try:
				x = float(line[30:38])
				y = float(line[38:46])
				z = float(line[46:54])
			except ValueError:
				continue
			coordinates.append((x, y, z))
	return coordinates

def compute_box(receptor_file: Path, padding: float) -> Box:
	coordinates = extract_pdb_coordinates(receptor_file)
	if not coordinates:
		raise ValueError(f"Não encontrei átomos válidos em `{receptor_file}`")
	x_values = [point[0] for point in coordinates]
	y_values = [point[1] for point in coordinates]
	z_values = [point[2] for point in coordinates]
	center_x = (min(x_values) + max(x_values)) / 2.0
	center_y = (min(y_values) + max(y_values)) / 2.0
	center_z = (min(z_values) + max(z_values)) / 2.0
	size_x = max(max(x_values) - min(x_values) + padding, 20.0)
	size_y = max(max(y_values) - min(y_values) + padding, 20.0)
	size_z = max(max(z_values) - min(z_values) + padding, 20.0)
	return Box(center_x, center_y, center_z, size_x, size_y, size_z)

def write_vina_config(config_file: Path, receptor_pdbqt: Path, ligand_pdbqt: Path, box: Box, exhaustiveness: int, num_modes: int, energy_range: int, output_pdbqt: Path) -> None:
	config_file.parent.mkdir(parents=True, exist_ok=True)
	config_file.write_text(
		"\n".join(
			[
				f"receptor = {receptor_pdbqt}",
				f"ligand = {ligand_pdbqt}",
				"",
				f"center_x = {box.center_x:.3f}",
				f"center_y = {box.center_y:.3f}",
				f"center_z = {box.center_z:.3f}",
				"",
				f"size_x = {box.size_x:.3f}",
				f"size_y = {box.size_y:.3f}",
				f"size_z = {box.size_z:.3f}",
				"",
				f"exhaustiveness = {exhaustiveness}",
				f"num_modes = {num_modes}",
				f"energy_range = {energy_range}",
				"",
				f"out = {output_pdbqt}",
				"",
			]
		),
		encoding="utf-8",
	)

def run_command(command: list[str], cwd: Path | None = None) -> None:
	completed = subprocess.run(command, cwd=cwd, check=False, text=True, capture_output=True)
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

def convert_to_pdbqt(input_file: Path, output_file: Path, obabel_exe: str, input_format: str, output_format: str, extra_args: Iterable[str] = ()) -> None:
	output_file.parent.mkdir(parents=True, exist_ok=True)
	command = [obabel_exe, f"-{input_format}", str(input_file), f"-{output_format}", "-O", str(output_file)]
	command.extend(extra_args)
	run_command(command)

def prepare_receptor_pdbqt(receptor_file: Path, output_file: Path, obabel_exe: str) -> None:
	convert_to_pdbqt(receptor_file, output_file, obabel_exe, "ipdb", "opdbqt", extra_args=("-xr",))

def prepare_ligand_pdbqt(ligand_file: Path, output_file: Path, obabel_exe: str) -> None:
	input_format = f"i{ligand_file.suffix.lstrip('.').lower()}"
	convert_to_pdbqt(ligand_file, output_file, obabel_exe, input_format, "opdbqt", extra_args=("--partialcharge", "gasteiger"))

def resolve_ligands(raw_structure: Path, ligand_ids: Iterable[str] | None, no_auto_ligands: bool) -> list[str]:
	manual_ligands = base_module.parse_ligand_ids(ligand_ids)
	if manual_ligands:
		return manual_ligands
	if no_auto_ligands:
		return []
	if raw_structure.suffix.lower() == ".cif":
		return base_module.detect_ligands_from_mmcif(raw_structure)
	return base_module.detect_ligands(raw_structure)

def resolve_structure(pdb_id: str, chain_id: str, output_dir: Path) -> tuple[Path, Path]:
	raw_structure = base_module.download_structure(pdb_id, output_dir)
	if raw_structure.suffix.lower() == ".cif":
		receptor_file = base_module.extract_chain_from_mmcif(raw_structure, chain_id, output_dir, pdb_id)
	else:
		receptor_file = base_module.extract_chain_from_pdb(raw_structure, chain_id, output_dir, pdb_id)
	return raw_structure, receptor_file

def build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(description="Executa docking do ligante contra a proteína baixada por `base.py`.")
	parser.add_argument("--pdb-id", required=True, help="Identificador PDB, por exemplo 9THJ")
	parser.add_argument("--chain", required=True, help="Cadeia do receptor, por exemplo A")
	parser.add_argument("--output-dir", required=True, help="Diretório raiz de saída")
	parser.add_argument("--ligand-id", action="append", dest="ligand_ids", help="Ligante do PDB. Pode repetir ou usar vírgula.")
	parser.add_argument("--ligand-file", action="append", dest="ligand_files", help="Arquivo de ligante local em SDF/MOL2/PDB/PDBQT. Pode repetir.")
	parser.add_argument("--no-auto-ligands", action="store_true", help="Desliga a detecção automática de ligantes.")
	parser.add_argument("--center-x", help="Centro X da caixa de docking")
	parser.add_argument("--center-y", help="Centro Y da caixa de docking")
	parser.add_argument("--center-z", help="Centro Z da caixa de docking")
	parser.add_argument("--size-x", help="Tamanho X da caixa de docking")
	parser.add_argument("--size-y", help="Tamanho Y da caixa de docking")
	parser.add_argument("--size-z", help="Tamanho Z da caixa de docking")
	parser.add_argument("--padding", type=float, default=8.0, help="Padding aplicado ao box automático")
	parser.add_argument("--exhaustiveness", type=int, default=8, help="Parâmetro de busca do Vina")
	parser.add_argument("--num-modes", type=int, default=9, help="Número de modos retornados pelo Vina")
	parser.add_argument("--energy-range", type=int, default=3, help="Faixa energética do Vina")
	parser.add_argument("--vina-exe", default="vina", help="Executável do AutoDock Vina")
	parser.add_argument("--obabel-exe", default="obabel", help="Executável do Open Babel")
	parser.add_argument("--prepare-only", action="store_true", help="Só prepara os arquivos PDBQT e a configuração")
	return parser

def prepare_ligand_sources(output_dir: Path, ligand_ids: list[str], ligand_files: list[str] | None) -> list[tuple[str, Path]]:
	prepared: list[tuple[str, Path]] = []
	ligand_dir = output_dir / "ligands"
	ligand_dir.mkdir(parents=True, exist_ok=True)
	for ligand_id in ligand_ids:
		if ligand_id.endswith(".pdbqt") or ligand_id.endswith(".pdb") or ligand_id.endswith(".sdf") or ligand_id.endswith(".mol2"):
			label = Path(ligand_id).stem
			prepared.append((label, Path(ligand_id)))
			continue
		ligand_path = base_module.download_ligand(ligand_id, output_dir)
		prepared.append((ligand_id, ligand_path))
	for ligand_file in ligand_files or []:
		path = Path(ligand_file).expanduser().resolve()
		if not path.exists():
			raise FileNotFoundError(f"Ligante não encontrado: {path}")
		prepared.append((path.stem, path))
	return prepared

def run_docking_for_ligand(
	label: str,
	ligand_source: Path,
	receptor_file: Path,
	box: Box,
	args: argparse.Namespace,
	base_output_dir: Path,
) -> Path:
	ligand_run_dir = base_output_dir / "docking" / label
	ligand_run_dir.mkdir(parents=True, exist_ok=True)
	receptor_pdbqt = ligand_run_dir / "receptor.pdbqt"
	ligand_pdbqt = ligand_run_dir / "ligand.pdbqt"
	config_file = ligand_run_dir / "config.txt"
	output_pdbqt = ligand_run_dir / "out.pdbqt"
	log_file = ligand_run_dir / "log.txt"
	prepare_receptor_pdbqt(receptor_file, receptor_pdbqt, args.obabel_exe)
	prepare_ligand_pdbqt(ligand_source, ligand_pdbqt, args.obabel_exe)
	write_vina_config(config_file, receptor_pdbqt, ligand_pdbqt, box, args.exhaustiveness, args.num_modes, args.energy_range, output_pdbqt)
	if args.prepare_only:
		return ligand_run_dir
	command = [
		args.vina_exe,
		"--config",
		str(config_file),
	]
	completed = subprocess.run(command, check=False, text=True, capture_output=True)
	if completed.returncode != 0:
		raise RuntimeError(
			"\n".join(
				[
					f"Vina falhou para `{label}`",
					completed.stdout.strip(),
					completed.stderr.strip(),
				]
			)
		)
	if completed.stdout.strip():
		print(completed.stdout.strip())
	if completed.stderr.strip():
		print(completed.stderr.strip(), file=sys.stderr)
	log_file.write_text(
		"\n".join(filter(None, [completed.stdout.strip(), completed.stderr.strip()])),
		encoding="utf-8",
	)
	if not output_pdbqt.exists():
		raise FileNotFoundError(f"Vina terminou sem gerar `{output_pdbqt}`")
	return ligand_run_dir

def main(argv: list[str] | None = None) -> int:
	parser = build_parser()
	args = parser.parse_args(argv)
	pdb_id = args.pdb_id.strip().upper()
	chain_id = args.chain.strip().upper()
	output_dir = Path(args.output_dir)
	try:
		ensure_executable(args.obabel_exe)
		if not args.prepare_only:
			ensure_executable(args.vina_exe)
		print(f"Preparando receptor a partir de {pdb_id} cadeia {chain_id}...")
		raw_structure, receptor_file = resolve_structure(pdb_id, chain_id, output_dir)
		ligand_ids = resolve_ligands(raw_structure, args.ligand_ids, args.no_auto_ligands)
		ligand_sources = prepare_ligand_sources(output_dir, ligand_ids, args.ligand_files)
		if not ligand_sources:
			raise ValueError("Nenhum ligante foi informado ou detectado")
		if all(value is None for value in (parse_float_or_none(args.center_x), parse_float_or_none(args.center_y), parse_float_or_none(args.center_z), parse_float_or_none(args.size_x), parse_float_or_none(args.size_y), parse_float_or_none(args.size_z))):
			box = compute_box(receptor_file, args.padding)
		else:
			missing = [name for name, value in {
				"center_x": args.center_x,
				"center_y": args.center_y,
				"center_z": args.center_z,
				"size_x": args.size_x,
				"size_y": args.size_y,
				"size_z": args.size_z,
			}.items() if value is None]
			if missing:
				raise ValueError("Para caixa manual, informe todos os parâmetros: " + ", ".join(missing))
			box = Box(
				parse_float_or_none(args.center_x) or 0.0,
				parse_float_or_none(args.center_y) or 0.0,
				parse_float_or_none(args.center_z) or 0.0,
				parse_float_or_none(args.size_x) or 20.0,
				parse_float_or_none(args.size_y) or 20.0,
				parse_float_or_none(args.size_z) or 20.0,
			)
		print(f"Box: center=({box.center_x:.3f}, {box.center_y:.3f}, {box.center_z:.3f}) size=({box.size_x:.3f}, {box.size_y:.3f}, {box.size_z:.3f})")
		for label, ligand_source in ligand_sources:
			print(f"Executando docking para {label}...")
			run_dir = run_docking_for_ligand(label, ligand_source, receptor_file, box, args, output_dir)
			print(f"Docking concluído em: {run_dir}")
		return 0
	except Exception as exc:
		print(f"Erro: {exc}", file=sys.stderr)
		return 1

if __name__ == "__main__":
	raise SystemExit(main())