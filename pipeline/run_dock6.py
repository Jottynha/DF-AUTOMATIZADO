from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))

import base
from dock6 import analysis as dock6_analysis
from dock6 import docking as dock6_docking
from dock6.prepare import prepare_dock6_inputs


def add_common_executable_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--pdb2pqr-exe", default="pdb2pqr")
    parser.add_argument("--obabel-exe", default="obabel")
    parser.add_argument("--dms-exe", default="dms")
    parser.add_argument("--sphgen-exe", default="sphgen")
    parser.add_argument("--sphere-selector-exe", default="sphere_selector")
    parser.add_argument("--showbox-exe", default="showbox")
    parser.add_argument("--grid-exe", default="grid")
    parser.add_argument("--dock6-exe", default="dock6")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Entrada única para preparar, executar e analisar o fluxo DOCK6."
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare",
        help="Prepara as entradas químicas/geométricas do DOCK6.",
    )
    prepare_parser.add_argument("pdb_id")
    prepare_parser.add_argument("chain")
    prepare_parser.add_argument("ligand_id")
    prepare_parser.add_argument("output_dir")
    prepare_parser.add_argument("--ph", type=float, default=7.4)
    add_common_executable_args(prepare_parser)

    pipeline_parser = subparsers.add_parser(
        "pipeline",
        help="Executa o pipeline completo de preparação, docking e análise.",
    )
    pipeline_parser.add_argument("pdb_id")
    pipeline_parser.add_argument("chain")
    pipeline_parser.add_argument("ligand_id")
    pipeline_parser.add_argument("output_dir")
    pipeline_parser.add_argument("--dock6-param-dir", required=True)
    pipeline_parser.add_argument("--ph", type=float, default=7.4)
    pipeline_parser.add_argument("--sphere-radius", type=float, default=10.0)
    pipeline_parser.add_argument("--prepare-only", action="store_true")
    add_common_executable_args(pipeline_parser)

    redocking_parser = subparsers.add_parser(
        "redocking",
        help="Executa o redocking com DOCK6 a partir de entradas opcionais já preparadas.",
    )
    redocking_parser.add_argument("pdb_id")
    redocking_parser.add_argument("chain")
    redocking_parser.add_argument("ligand_id")
    redocking_parser.add_argument("output_dir")
    redocking_parser.add_argument("--receptor-ms", required=True)
    redocking_parser.add_argument("--receptor-mol2")
    redocking_parser.add_argument("--ligand-ref-mol2")
    redocking_parser.add_argument("--ligand-dock-mol2")
    redocking_parser.add_argument("--dock6-param-dir", required=True)
    redocking_parser.add_argument("--sphere-radius", type=float, default=10.0)
    redocking_parser.add_argument("--prepare-only", action="store_true")
    add_common_executable_args(redocking_parser)

    return parser


def run_prepare(args: argparse.Namespace) -> int:
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


def run_pipeline(args: argparse.Namespace) -> int:
    output_dir = Path(args.output_dir)
    prepared_dir = output_dir / "prepared_inputs"
    docking_output_dir = output_dir / "run"

    print()
    print("===" * 20)
    print("EXECUTANDO O PIPELINE DOCK6 COMPLETO")
    print("===" * 20)
    print()

    prepared = prepare_dock6_inputs(
        pdb_id=args.pdb_id,
        chain_id=args.chain,
        ligand_id=args.ligand_id,
        output_dir=prepared_dir,
        pdb2pqr_exe=args.pdb2pqr_exe,
        obabel_exe=args.obabel_exe,
        dms_exe=args.dms_exe,
        ph=args.ph,
    )

    print()
    print("[2/4] Preparando/rodando DOCK6")
    raw_structure = base.download_structure(args.pdb_id, docking_output_dir)
    if raw_structure.suffix.lower() == ".cif":
        receptor_file = base.extract_chain_from_mmcif(raw_structure, args.chain, docking_output_dir, args.pdb_id)
    else:
        receptor_file = base.extract_chain_from_pdb(raw_structure, args.chain, docking_output_dir, args.pdb_id)

    run_dir = dock6_docking.run_dock6_pipeline_for_ligand(
        ligand_label=args.ligand_id,
        ligand_source=Path(prepared.ligand_dock_mol2),
        receptor_pdb=receptor_file,
        receptor_ms=Path(prepared.receptor_ms),
        receptor_mol2_source=Path(prepared.receptor_fixed_mol2),
        ligand_ref_mol2_source=Path(prepared.ligand_ref_mol2),
        ligand_dock_mol2_source=Path(prepared.ligand_dock_mol2),
        base_output_dir=docking_output_dir,
        dock6_param_dir=Path(args.dock6_param_dir),
        obabel_exe=args.obabel_exe,
        sphgen_exe=args.sphgen_exe,
        sphere_selector_exe=args.sphere_selector_exe,
        showbox_exe=args.showbox_exe,
        grid_exe=args.grid_exe,
        dock6_exe=args.dock6_exe,
        sphere_radius=args.sphere_radius,
        prepare_only=args.prepare_only,
    )
    print(f"Diretório da rodada DOCK6: {run_dir}")
    print()

    if args.prepare_only:
        print("Modo prepare-only ativado. Encerrando antes de grid/dock6.")
        return 0

    print("[3/4] Analisando resultado")
    result = dock6_analysis.analyze_dock6_result(
        dock6_run_dir=run_dir,
        crystal_ligand_mol2=Path(prepared.ligand_ref_mol2),
    )

    score = result.get("grid_score")
    rmsd = result.get("rmsd")

    print(f"Grid Score DOCK6: {score:.4f}" if score is not None else "Grid Score DOCK6: N/A")
    print(f"RMSD: {rmsd:.2f} Å" if rmsd is not None else "RMSD: N/A")

    if rmsd is not None and rmsd <= 2.0:
        print("Redocking bem-sucedido: RMSD <= 2.0 Å")
    elif rmsd is not None:
        print("RMSD elevado: verificar pose manualmente.")
    print()

    print("[4/4] Salvando JSON")
    final_json = output_dir / "resultado_pipeline_dock6_completo.json"
    final_json.parent.mkdir(parents=True, exist_ok=True)
    final_json.write_text(
        json.dumps(
            {
                "prepared_inputs": asdict(prepared),
                "docking_result": result,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Resultado salvo em: {final_json}")
    return 0


def run_redocking(args: argparse.Namespace) -> int:
    output_dir_path = Path(args.output_dir)

    print()
    print("Pipeline da 2ª Entrega: Redocking com DOCK6 + RMSD")
    print()

    print(f"[1/8] Baixando estrutura {args.pdb_id}")
    raw_structure = base.download_structure(args.pdb_id, output_dir_path)
    print(f"Estrutura salva em: {raw_structure}")
    print()

    print(f"[2/8] Preparando proteína cadeia {args.chain}")
    if raw_structure.suffix.lower() == ".cif":
        receptor_file = base.extract_chain_from_mmcif(raw_structure, args.chain, output_dir_path, args.pdb_id)
    else:
        receptor_file = base.extract_chain_from_pdb(raw_structure, args.chain, output_dir_path, args.pdb_id)

    print(f"Proteína salva em: {receptor_file}")
    print()

    print(f"[3/8] Baixando ligante cristalográfico {args.ligand_id}")
    ligand_file = base.download_ligand(args.ligand_id, output_dir_path)
    print(f"Ligante salvo em: {ligand_file}")
    print()

    print("[4/8] Preparando e executando fluxo DOCK6")
    run_dir = dock6_docking.run_dock6_pipeline_for_ligand(
        ligand_label=args.ligand_id,
        ligand_source=ligand_file,
        receptor_pdb=receptor_file,
        receptor_ms=Path(args.receptor_ms),
        receptor_mol2_source=Path(args.receptor_mol2) if args.receptor_mol2 else None,
        ligand_ref_mol2_source=Path(args.ligand_ref_mol2) if args.ligand_ref_mol2 else None,
        ligand_dock_mol2_source=Path(args.ligand_dock_mol2) if args.ligand_dock_mol2 else None,
        base_output_dir=output_dir_path,
        dock6_param_dir=Path(args.dock6_param_dir),
        obabel_exe=args.obabel_exe,
        sphgen_exe=args.sphgen_exe,
        sphere_selector_exe=args.sphere_selector_exe,
        showbox_exe=args.showbox_exe,
        grid_exe=args.grid_exe,
        dock6_exe=args.dock6_exe,
        sphere_radius=args.sphere_radius,
        prepare_only=args.prepare_only,
    )

    print(f"Diretório DOCK6: {run_dir}")
    print()

    if args.prepare_only:
        print("[5/8] Modo prepare-only ativado")
        print("Arquivos de entrada preparados, mas grid/dock6 não foram executados.")
        print()
        return 0

    print("[5/8] Preparando ligante cristalográfico para RMSD")
    crystal_ligand_mol2 = run_dir / "ligand_ref.mol2"

    if not crystal_ligand_mol2.exists():
        dock6_docking.prepare_ligand_mol2(
            ligand_file,
            crystal_ligand_mol2,
            obabel_exe=args.obabel_exe,
        )

    print(f"Ligante cristalográfico MOL2 salvo em: {crystal_ligand_mol2}")
    print()

    print("[6/8] Analisando score e RMSD")
    result = dock6_analysis.analyze_dock6_result(
        dock6_run_dir=run_dir,
        crystal_ligand_mol2=crystal_ligand_mol2,
    )

    score = result.get("grid_score")
    rmsd = result.get("rmsd")

    if score is not None:
        print(f"Grid Score DOCK6: {score:.4f}")
    else:
        print("Grid Score DOCK6: N/A")

    if rmsd is not None:
        print(f"RMSD: {rmsd:.2f} Å")
        if rmsd <= 2.0:
            print("Redocking bem-sucedido: RMSD <= 2.0 Å")
        else:
            print("RMSD elevado: verificar pose manualmente.")
    else:
        print("RMSD: N/A")

    print()

    print("[7/8] Salvando resultado JSON")
    results_json = output_dir_path / "resultados_redocking_dock6.json"
    results_json.write_text(json.dumps([result], indent=2), encoding="utf-8")
    print(f"Resultados salvos em: {results_json}")
    print()

    print("[8/8] Resumo")
    print(f"-> PDB ID: {args.pdb_id}")
    print(f"-> Cadeia: {args.chain}")
    print(f"-> Ligante cristalográfico: {args.ligand_id}")
    print(f"-> Grid Score: {score:.4f}" if score is not None else "-> Grid Score: N/A")
    print(f"-> RMSD: {rmsd:.2f} Å" if rmsd is not None else "-> RMSD: N/A")
    print(f"-> Diretório de saída: {output_dir_path}")
    print()

    print("Pipeline DOCK6 completo.")
    print()
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command == "prepare":
            return run_prepare(args)
        if args.command == "pipeline":
            return run_pipeline(args)
        if args.command == "redocking":
            return run_redocking(args)
        raise ValueError(f"Comando desconhecido: {args.command}")
    except Exception as exc:
        print(f"Erro: {exc}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())