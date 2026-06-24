from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path


def _add_repo_to_path() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root))


def _add_common_executable_args(parser: argparse.ArgumentParser) -> None:
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Pipeline organizado para DOCK6: redocking e triagem de biblioteca. Saídas brutas em <saida>/dock6_out/."
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    redocking_parser = subparsers.add_parser(
        "redocking",
        help="Executa redocking completo com DOCK6.",
    )
    redocking_parser.add_argument("pdb_id")
    redocking_parser.add_argument("chain")
    redocking_parser.add_argument("ligand_id")
    redocking_parser.add_argument("output_dir")
    redocking_parser.add_argument("--prepare-only", action="store_true")
    _add_common_executable_args(redocking_parser)

    library_parser = subparsers.add_parser(
        "library",
        help="Faz docking de biblioteca de ligantes com DOCK6.",
    )
    library_parser.add_argument("pdb_id")
    library_parser.add_argument("chain")
    library_parser.add_argument("ligand_ref", help="Ligante cristalográfico usado para definir o sítio ativo.")
    library_parser.add_argument("ligand_library", help="Diretório/arquivo, self/ref/crystal, spec pubchem:, chembl:, zinc: ou pubchem-file:, chembl-file:, zinc-file:.")
    library_parser.add_argument("output_dir")
    library_parser.add_argument("--calculate-rmsd-against-ref", action="store_true")
    library_parser.add_argument("--stop-on-error", action="store_true")
    _add_common_executable_args(library_parser)

    download_parser = subparsers.add_parser(
        "download-library",
        help="Baixa/converte ligantes sem executar docking.",
    )
    download_parser.add_argument("ligand_library", help="Arquivo/diretório, spec pubchem:, chembl:, zinc: ou pubchem-file:, chembl-file:, zinc-file:.")
    download_parser.add_argument("output_dir")
    download_parser.add_argument("--obabel-exe", default="obabel")

    similar_parser = subparsers.add_parser(
        "download-similar",
        help="Busca similares no PubChem por SMILES e prepara biblioteca sem docking.",
    )
    similar_parser.add_argument("query_smiles")
    similar_parser.add_argument("output_dir")
    similar_parser.add_argument("--threshold", type=int, default=90)
    similar_parser.add_argument("--max-ligands", type=int, default=1000)
    similar_parser.add_argument("--obabel-exe", default="obabel")

    screen_library_parser = subparsers.add_parser(
        "screen-library",
        help="Baixa/prepara uma biblioteca e executa triagem virtual DOCK6.",
    )
    screen_library_parser.add_argument("pdb_id")
    screen_library_parser.add_argument("chain")
    screen_library_parser.add_argument("ligand_ref")
    screen_library_parser.add_argument("ligand_library")
    screen_library_parser.add_argument("output_dir")
    screen_library_parser.add_argument("--calculate-rmsd-against-ref", action="store_true")
    screen_library_parser.add_argument("--stop-on-error", action="store_true")
    _add_common_executable_args(screen_library_parser)

    screen_similar_parser = subparsers.add_parser(
        "screen-similar",
        help="Busca similares no PubChem por SMILES, prepara ligantes e executa triagem virtual DOCK6.",
    )
    screen_similar_parser.add_argument("pdb_id")
    screen_similar_parser.add_argument("chain")
    screen_similar_parser.add_argument("ligand_ref")
    screen_similar_parser.add_argument("query_smiles")
    screen_similar_parser.add_argument("output_dir")
    screen_similar_parser.add_argument("--threshold", type=int, default=70)
    screen_similar_parser.add_argument("--max-ligands", type=int, default=1000)
    screen_similar_parser.add_argument("--calculate-rmsd-against-ref", action="store_true")
    screen_similar_parser.add_argument("--stop-on-error", action="store_true")
    _add_common_executable_args(screen_similar_parser)

    search_parser = subparsers.add_parser(
        "search",
        help="Busca similares no PubChem por SMILES e faz docking com DOCK6.",
    )
    search_parser.add_argument("pdb_id")
    search_parser.add_argument("chain")
    search_parser.add_argument("ligand_ref", help="Ligante cristalográfico usado para definir o sítio ativo.")
    search_parser.add_argument("query_smiles")
    search_parser.add_argument("output_dir")
    search_parser.add_argument("--threshold", type=int, default=90)
    search_parser.add_argument("--max-ligands", type=int, default=20)
    search_parser.add_argument("--calculate-rmsd-against-ref", action="store_true")
    search_parser.add_argument("--stop-on-error", action="store_true")
    _add_common_executable_args(search_parser)

    return parser


def _print_top_results(results: list, max_rows: int = 10) -> None:
    ranked = [result for result in results if getattr(result, "grid_score", None) is not None]
    ranked.sort(key=lambda item: item.grid_score)

    if not ranked:
        print("Nenhum score DOCK6 válido foi encontrado.")
        return

    print("\nTop resultados:")
    for result in ranked[:max_rows]:
        score = f"{result.grid_score:.4f}" if result.grid_score is not None else "N/A"
        rmsd = f"{result.rmsd:.2f} Å" if result.rmsd is not None else "N/A"
        print(f"#{result.rank or '-'} {result.ligand_label} | Grid Score: {score} | RMSD: {rmsd}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _add_repo_to_path()

    try:
        if args.command == "download-library":
            from src.dock6 import dock_library

            ligand_files = dock_library.prepare_library_only(
                ligand_library=args.ligand_library,
                output_dir=Path(args.output_dir),
                obabel_exe=args.obabel_exe,
            )
            print(f"Biblioteca preparada com {len(ligand_files)} ligantes.")
            print(f"Diretório: {Path(args.output_dir) / 'ligands'}")
            print(f"Manifesto CSV: {Path(args.output_dir) / 'ligands_manifest.csv'}")
            print(f"Manifesto JSON: {Path(args.output_dir) / 'ligands_manifest.json'}")
            return 0

        if args.command == "download-similar":
            from src.dock6 import dock_library

            ligand_files = dock_library.prepare_similar_library_only(
                query_smiles=args.query_smiles,
                output_dir=Path(args.output_dir),
                threshold=args.threshold,
                max_ligands=args.max_ligands,
                obabel_exe=args.obabel_exe,
            )
            print(f"Biblioteca similar preparada com {len(ligand_files)} ligantes.")
            print(f"CIDs: {Path(args.output_dir) / 'pubchem_similar_cids.txt'}")
            print(f"Diretório: {Path(args.output_dir) / 'ligands'}")
            print(f"Manifesto CSV: {Path(args.output_dir) / 'ligands_manifest.csv'}")
            print(f"Manifesto JSON: {Path(args.output_dir) / 'ligands_manifest.json'}")
            return 0

        if args.command == "screen-library":
            from src.dock6 import dock_library

            output_dir = Path(args.output_dir)
            library_dir = output_dir / "downloaded_library"

            print("[1/2] Baixando/preparando biblioteca de ligantes")
            ligand_files = dock_library.prepare_library_only(
                ligand_library=args.ligand_library,
                output_dir=library_dir,
                obabel_exe=args.obabel_exe,
            )
            print(f"Biblioteca preparada com {len(ligand_files)} ligantes.")
            print(f"Manifesto CSV: {library_dir / 'ligands_manifest.csv'}")
            print()

            print("[2/2] Executando triagem virtual com DOCK6")
            results = dock_library.dock_library(
                pdb_id=args.pdb_id,
                chain=args.chain,
                ligand_ref=args.ligand_ref,
                ligand_library=str(library_dir / "ligands"),
                output_dir=output_dir,
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
                calculate_rmsd_against_ref=args.calculate_rmsd_against_ref,
                continue_on_error=not args.stop_on_error,
            )
            _print_top_results(results)
            print(f"\nResumo CSV: {output_dir / 'dock6_results_summary.csv'}")
            print(f"Resumo JSON: {output_dir / 'dock6_results_summary.json'}")
            print(f"Biblioteca baixada/preparada: {library_dir / 'ligands'}")
            return 0

        if args.command == "screen-similar":
            from src.dock6 import dock_library

            output_dir = Path(args.output_dir)
            library_dir = output_dir / "downloaded_library"

            print("[1/2] Buscando e preparando biblioteca similar no PubChem")
            ligand_files = dock_library.prepare_similar_library_only(
                query_smiles=args.query_smiles,
                output_dir=library_dir,
                threshold=args.threshold,
                max_ligands=args.max_ligands,
                obabel_exe=args.obabel_exe,
            )
            print(f"Biblioteca preparada com {len(ligand_files)} ligantes.")
            print(f"CIDs: {library_dir / 'pubchem_similar_cids.txt'}")
            print(f"Manifesto CSV: {library_dir / 'ligands_manifest.csv'}")
            print()

            print("[2/2] Executando triagem virtual com DOCK6")
            results = dock_library.dock_library(
                pdb_id=args.pdb_id,
                chain=args.chain,
                ligand_ref=args.ligand_ref,
                ligand_library=str(library_dir / "ligands"),
                output_dir=output_dir,
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
                calculate_rmsd_against_ref=args.calculate_rmsd_against_ref,
                continue_on_error=not args.stop_on_error,
            )
            _print_top_results(results)
            print(f"\nResumo CSV: {output_dir / 'dock6_results_summary.csv'}")
            print(f"Resumo JSON: {output_dir / 'dock6_results_summary.json'}")
            print(f"Biblioteca baixada: {library_dir / 'ligands'}")
            return 0

        if args.command == "redocking":
            from src.dock6 import redocking

            result = redocking.redock(
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

            print(f"Diretório DOCK6: {result.run_dir}")
            print(f"Grid Score DOCK6: {result.grid_score:.4f}" if result.grid_score is not None else "Grid Score DOCK6: N/A")
            print(f"RMSD: {result.rmsd:.2f} Å" if result.rmsd is not None else "RMSD: N/A")
            print(f"Resultado JSON: {result.result_json}")
            return 0

        if args.command == "library":
            from src.dock6 import dock_library

            results = dock_library.dock_library(
                pdb_id=args.pdb_id,
                chain=args.chain,
                ligand_ref=args.ligand_ref,
                ligand_library=args.ligand_library,
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
                calculate_rmsd_against_ref=args.calculate_rmsd_against_ref,
                continue_on_error=not args.stop_on_error,
            )
            _print_top_results(results)
            print(f"\nResumo CSV: {Path(args.output_dir) / 'dock6_results_summary.csv'}")
            print(f"Resumo JSON: {Path(args.output_dir) / 'dock6_results_summary.json'}")
            return 0

        if args.command == "search":
            from src.dock6 import dock_library

            results = dock_library.search_and_dock(
                pdb_id=args.pdb_id,
                chain=args.chain,
                ligand_ref=args.ligand_ref,
                query_smiles=args.query_smiles,
                output_dir=args.output_dir,
                dock6_param_dir=args.dock6_param_dir,
                threshold=args.threshold,
                max_ligands=args.max_ligands,
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
                calculate_rmsd_against_ref=args.calculate_rmsd_against_ref,
                continue_on_error=not args.stop_on_error,
            )
            _print_top_results(results)
            print(f"\nResumo CSV: {Path(args.output_dir) / 'dock6_results_summary.csv'}")
            print(f"Resumo JSON: {Path(args.output_dir) / 'dock6_results_summary.json'}")
            return 0

        raise ValueError(f"Comando desconhecido: {args.command}")

    except Exception as exc:
        print(f"Erro durante pipeline DOCK6: {exc}", file=sys.stderr)
        import traceback

        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
