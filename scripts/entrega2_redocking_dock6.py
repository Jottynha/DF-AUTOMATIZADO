"""

"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))

import base
from dock6 import analysis as dock6_analysis
from dock6 import docking as dock6_docking


def main(
    pdb_id: str,
    chain_id: str,
    ligand_id: str,
    output_dir: str,
    receptor_ms: str,
    receptor_mol2: str | None,
    ligand_ref_mol2: str | None,
    ligand_dock_mol2: str | None,
    dock6_param_dir: str,
    sphere_radius: float,
    prepare_only: bool,
) -> int:
    output_dir_path = Path(output_dir)

    print()
    print("Pipeline da 2ª Entrega: Redocking com DOCK6 + RMSD")
    print()

    try:
        print(f"[1/8] Baixando estrutura {pdb_id}")
        raw_structure = base.download_structure(pdb_id, output_dir_path)
        print(f"Estrutura salva em: {raw_structure}")
        print()

        print(f"[2/8] Preparando proteína cadeia {chain_id}")
        if raw_structure.suffix.lower() == ".cif":
            receptor_file = base.extract_chain_from_mmcif(
                raw_structure,
                chain_id,
                output_dir_path,
                pdb_id,
            )
        else:
            receptor_file = base.extract_chain_from_pdb(
                raw_structure,
                chain_id,
                output_dir_path,
                pdb_id,
            )

        print(f"Proteína salva em: {receptor_file}")
        print()

        print(f"[3/8] Baixando ligante cristalográfico {ligand_id}")
        ligand_file = base.download_ligand(ligand_id, output_dir_path)
        print(f"Ligante salvo em: {ligand_file}")
        print()

        print("[4/8] Preparando e executando fluxo DOCK6")
        run_dir = dock6_docking.run_dock6_pipeline_for_ligand(
            ligand_label=ligand_id,
            ligand_source=ligand_file,
            receptor_pdb=receptor_file,
            receptor_ms=Path(receptor_ms),
            receptor_mol2_source=Path(receptor_mol2) if receptor_mol2 else None,
            ligand_ref_mol2_source=Path(ligand_ref_mol2) if ligand_ref_mol2 else None,
            ligand_dock_mol2_source=Path(ligand_dock_mol2) if ligand_dock_mol2 else None,
            base_output_dir=output_dir_path,
            dock6_param_dir=Path(dock6_param_dir),
            sphere_radius=sphere_radius,
            prepare_only=prepare_only,
        )

        print(f"Diretório DOCK6: {run_dir}")
        print()

        if prepare_only:
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
                obabel_exe="obabel",
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
        results_json.write_text(
            json.dumps([result], indent=2),
            encoding="utf-8",
        )
        print(f"Resultados salvos em: {results_json}")
        print()

        print("[8/8] Resumo")
        print(f"-> PDB ID: {pdb_id}")
        print(f"-> Cadeia: {chain_id}")
        print(f"-> Ligante cristalográfico: {ligand_id}")
        print(f"-> Grid Score: {score:.4f}" if score is not None else "-> Grid Score: N/A")
        print(f"-> RMSD: {rmsd:.2f} Å" if rmsd is not None else "-> RMSD: N/A")
        print(f"-> Diretório de saída: {output_dir_path}")
        print()

        print("Pipeline DOCK6 completo.")
        print()

        return 0

    except Exception as exc:
        print()
        print(f"X - Erro: {exc}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Executa pipeline da 2ª entrega: redocking com DOCK6 + análise."
    )

    parser.add_argument("pdb_id", help="Identificador PDB. Exemplo: 9THJ")
    parser.add_argument("chain", help="Cadeia da proteína. Exemplo: A")
    parser.add_argument("ligand_id", help="ID do ligante cristalográfico. Exemplo: A1JV9")
    parser.add_argument("output_dir", help="Diretório de saída")

    parser.add_argument(
        "--receptor-ms",
        required=True,
        help="Arquivo .ms do receptor, usado pelo sphgen. Exemplo: Backup/dock6/rec.ms",
    )

    parser.add_argument(
        "--receptor-mol2",
        help=(
            "Arquivo receptor MOL2 já preparado para DOCK6. "
            "Recomendado. Exemplo: Backup/dock6/receptor_fixed.mol2"
        ),
    )

    parser.add_argument(
        "--ligand-ref-mol2",
        help=(
            "Ligante cristalográfico em MOL2, usado para selecionar as spheres. "
            "Exemplo: Backup/dock6/ligand_ref.mol2"
        ),
    )

    parser.add_argument(
        "--ligand-dock-mol2",
        help=(
            "Ligante preparado para docking em MOL2. "
            "Exemplo: Backup/dock6/ligand_dock.mol2"
        ),
    )

    parser.add_argument(
        "--dock6-param-dir",
        required=True,
        help="Pasta de parâmetros do DOCK6 contendo vdw_AMBER_parm99.defn, flex.defn e flex_drive.tbl.",
    )

    parser.add_argument(
        "--sphere-radius",
        type=float,
        default=10.0,
        help="Raio para selecionar esferas próximas ao ligante cristalográfico.",
    )

    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepara arquivos, mas não executa grid/dock6.",
    )

    args = parser.parse_args()

    raise SystemExit(
        main(
            pdb_id=args.pdb_id,
            chain_id=args.chain,
            ligand_id=args.ligand_id,
            output_dir=args.output_dir,
            receptor_ms=args.receptor_ms,
            receptor_mol2=args.receptor_mol2,
            ligand_ref_mol2=args.ligand_ref_mol2,
            ligand_dock_mol2=args.ligand_dock_mol2,
            dock6_param_dir=args.dock6_param_dir,
            sphere_radius=args.sphere_radius,
            prepare_only=args.prepare_only,
        )
    )