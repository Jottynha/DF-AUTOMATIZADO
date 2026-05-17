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
from dock6.prepare import prepare_dock6_inputs
from dock6 import docking as dock6_docking
from dock6 import analysis as dock6_analysis


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Executa pipeline DOCK6 completo.")
    parser.add_argument("pdb_id")
    parser.add_argument("chain")
    parser.add_argument("ligand_id")
    parser.add_argument("output_dir")
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
    parser.add_argument("--prepare-only", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
    prepared_dir = output_dir / "prepared_inputs"
    docking_output_dir = output_dir / "run"

    print()
    print("===" * 20)
    print("EXECUTANDO O PIPELINE DOCK6 COMPLETO")
    print("===" * 20)
    print()

    try:
        print("[1/4] Preparando entradas químicas/geométricas do DOCK6")
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

    except Exception as exc:
        print(f"Erro: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())