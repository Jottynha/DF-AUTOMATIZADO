from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from . import analysis as dock6_analysis
from . import docking as dock6_docking
from .prepare import prepare_dock6_inputs


@dataclass(frozen=True)
class Dock6RedockingResult:
    pdb_id: str
    chain: str
    ligand_id: str
    output_dir: str
    prepared_inputs: dict[str, Any]
    run_dir: str
    scored_mol2: str | None
    dock_out: str | None
    grid_score: float | None
    rmsd: float | None
    success: bool
    result_json: str


def redock(
    pdb_id: str,
    chain: str,
    ligand_id: str,
    output_dir: str | Path,
    dock6_param_dir: str | Path,
    *,
    pdb2pqr_exe: str = "pdb2pqr",
    obabel_exe: str = "obabel",
    dms_exe: str = "dms",
    sphgen_exe: str = "sphgen",
    sphere_selector_exe: str = "sphere_selector",
    showbox_exe: str = "showbox",
    grid_exe: str = "grid",
    dock6_exe: str = "dock6",
    ph: float = 7.4,
    sphere_radius: float = 10.0,
    prepare_only: bool = False,
) -> Dock6RedockingResult:
    """
    Executa redocking completo com DOCK6.

    Entrada mínima: PDB ID, cadeia, ligante cristalográfico e diretório de saída.
    O fluxo gera automaticamente receptor, ligante de referência, ligante de docking,
    superfície molecular, esferas, caixa, grid, docking e análise de score/RMSD.
    """
    out = Path(output_dir)
    prepared_dir = out / "prepared_inputs"

    prepared = prepare_dock6_inputs(
        pdb_id=pdb_id,
        chain_id=chain,
        ligand_id=ligand_id,
        output_dir=prepared_dir,
        pdb2pqr_exe=pdb2pqr_exe,
        obabel_exe=obabel_exe,
        dms_exe=dms_exe,
        ph=ph,
    )

    run_dir = dock6_docking.run_dock6_pipeline_for_ligand(
        ligand_label=ligand_id,
        ligand_source=Path(prepared.ligand_dock_mol2),
        receptor_pdb=Path(prepared.receptor_clean_pdb),
        receptor_ms=Path(prepared.receptor_ms),
        receptor_mol2_source=Path(prepared.receptor_fixed_mol2),
        ligand_ref_mol2_source=Path(prepared.ligand_ref_mol2),
        ligand_dock_mol2_source=Path(prepared.ligand_dock_mol2),
        base_output_dir=out,
        dock6_param_dir=Path(dock6_param_dir),
        obabel_exe=obabel_exe,
        sphgen_exe=sphgen_exe,
        sphere_selector_exe=sphere_selector_exe,
        showbox_exe=showbox_exe,
        grid_exe=grid_exe,
        dock6_exe=dock6_exe,
        sphere_radius=sphere_radius,
        prepare_only=prepare_only,
    )

    if prepare_only:
        analysis_result = {
            "scored_mol2": None,
            "dock_out": None,
            "grid_score": None,
            "rmsd": None,
        }
    else:
        analysis_result = dock6_analysis.analyze_dock6_result(
            dock6_run_dir=run_dir,
            crystal_ligand_mol2=Path(prepared.ligand_ref_mol2),
        )

    result_json = out / "dock6_redocking_result.json"
    result_json.parent.mkdir(parents=True, exist_ok=True)

    result = Dock6RedockingResult(
        pdb_id=pdb_id.upper(),
        chain=chain,
        ligand_id=ligand_id.upper(),
        output_dir=str(out),
        prepared_inputs=asdict(prepared),
        run_dir=str(run_dir),
        scored_mol2=analysis_result.get("scored_mol2"),
        dock_out=analysis_result.get("dock_out"),
        grid_score=analysis_result.get("grid_score"),
        rmsd=analysis_result.get("rmsd"),
        success=bool(analysis_result.get("rmsd") is not None and analysis_result.get("rmsd") <= 2.0),
        result_json=str(result_json),
    )

    result_json.write_text(json.dumps(asdict(result), indent=2), encoding="utf-8")
    return result


def preprocess_and_dock(
    pdb_id: str,
    chain: str,
    ligand_id: str,
    output_dir: str | Path,
    dock6_param_dir: str | Path,
    **kwargs: Any,
) -> tuple[Path | None, float | None]:
    """
    Alias compatível com o estilo do módulo Vina.
    Retorna: caminho do MOL2 pontuado e RMSD.
    """
    result = redock(
        pdb_id=pdb_id,
        chain=chain,
        ligand_id=ligand_id,
        output_dir=output_dir,
        dock6_param_dir=dock6_param_dir,
        **kwargs,
    )
    scored = Path(result.scored_mol2) if result.scored_mol2 else None
    return scored, result.rmsd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Pipeline de redocking DOCK6 organizado.")
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

    try:
        result = redock(
            pdb_id=args.pdb_id,
            chain=args.chain,
            ligand_id=args.ligand_id,
            output_dir=args.output_dir,
            dock6_param_dir=args.dock6_param_dir,
            pdb2pqr_exe=args.pdb2pqr_exe,
            obabel_exe=args.obabel_exe,
            dms_exe=args.dms_exe,
            sphgen_exe=args.sphgen_exe,
            sphere_selector_exe=args.sphere_selector_exe,
            showbox_exe=args.showbox_exe,
            grid_exe=args.grid_exe,
            dock6_exe=args.dock6_exe,
            ph=args.ph,
            sphere_radius=args.sphere_radius,
            prepare_only=args.prepare_only,
        )

        print(json.dumps(asdict(result), indent=2))
        if result.grid_score is not None:
            print(f"Grid Score DOCK6: {result.grid_score:.4f}")
        if result.rmsd is not None:
            print(f"RMSD: {result.rmsd:.2f} Å")
        print(f"Resultado salvo em: {result.result_json}")
        return 0
    except Exception as exc:
        print(f"Erro durante redocking DOCK6: {exc}", file=sys.stderr)
        import traceback

        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
