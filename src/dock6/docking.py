from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Dock6Paths:
    run_dir: Path
    receptor_mol2: Path
    ligand_mol2: Path
    receptor_ms: Path
    sphgen_in: Path
    all_spheres: Path
    selected_spheres: Path
    showbox_in: Path
    box_pdb: Path
    grid_in: Path
    grid_prefix: Path
    dock_in: Path
    dock_out: Path
    scored_mol2: Path


def ensure_executable(command: str) -> None:
    if shutil.which(command) is None:
        raise FileNotFoundError(f"Não encontrei o executável `{command}` no PATH")


def run_command(
    command: list[str],
    cwd: Path | None = None,
    stdin_text: str | None = None,
    log_file: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
        input=stdin_text,
        text=True,
        capture_output=True,
        check=False,
    )

    log_text = "\n".join(
        part for part in [completed.stdout.strip(), completed.stderr.strip()] if part
    )

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


def convert_to_mol2(
    input_file: Path,
    output_file: Path,
    obabel_exe: str,
    extra_args: Iterable[str] = (),
) -> None:
    if not input_file.exists():
        raise FileNotFoundError(f"Arquivo de entrada não encontrado: {input_file}")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    input_format = input_file.suffix.lstrip(".").lower()
    if not input_format:
        raise ValueError(f"Não consegui identificar o formato de {input_file}")

    command = [
        obabel_exe,
        f"-i{input_format}",
        str(input_file),
        "-omol2",
        "-O",
        str(output_file),
    ]
    command.extend(extra_args)

    run_command(command)


def prepare_receptor_mol2(receptor_pdb: Path, output_file: Path, obabel_exe: str) -> None:
    # Mantido como fallback. Para DOCK6, prefira passar --receptor-mol2 já preparado.
    convert_to_mol2(
        receptor_pdb,
        output_file,
        obabel_exe,
        extra_args=("-xr",),
    )


def prepare_ligand_mol2(ligand_file: Path, output_file: Path, obabel_exe: str) -> None:
    convert_to_mol2(
        ligand_file,
        output_file,
        obabel_exe,
        extra_args=("--partialcharge", "gasteiger"),
    )


def make_paths(base_output_dir: Path, ligand_label: str) -> Dock6Paths:
    run_dir = base_output_dir / "dock6_out" / ligand_label

    return Dock6Paths(
        run_dir=run_dir,
        receptor_mol2=run_dir / "receptor_fixed.mol2",
        ligand_mol2=run_dir / "ligand_dock.mol2",
        receptor_ms=run_dir / "rec.ms",
        sphgen_in=run_dir / "INSPH",
        all_spheres=run_dir / "rec.sph",
        selected_spheres=run_dir / "selected_spheres.sph",
        showbox_in=run_dir / "box.in",
        box_pdb=run_dir / "rec_box.pdb",
        grid_in=run_dir / "grid.in",
        grid_prefix=run_dir / "grid",
        dock_in=run_dir / "dock6.in",
        dock_out=run_dir / "dock6.out",
        scored_mol2=run_dir / "scored.mol2",
    )


def clean_previous_run(paths: Dock6Paths) -> None:
    """
    Remove arquivos que programas do DOCK6 não costumam sobrescrever.
    Isso evita erros como: Cannot open file 'OUTSPH': Arquivo existe.
    """
    old_files = [
        paths.run_dir / "OUTSPH",
        paths.all_spheres,
        paths.selected_spheres,
        paths.box_pdb,
        paths.run_dir / "selected_spheres.pdb",
        paths.run_dir / "grid.bmp",
        paths.run_dir / "grid.nrg",
        paths.run_dir / "grid.cnt",
        paths.run_dir / "grid.out",
        paths.run_dir / "grid_command.log",
        paths.dock_out,
        paths.run_dir / "dock6_command.log",
        paths.run_dir / "dock6_run_scored.mol2",
        paths.run_dir / "dock6_pose_.mol2",
        paths.scored_mol2,
        paths.run_dir / "ligand_ref.mol2",
        paths.run_dir / "sphgen.log",
        paths.run_dir / "sphere_selector.log",
        paths.run_dir / "showbox.log",
    ]

    for old_file in old_files:
        if old_file.exists():
            old_file.unlink()


def copy_or_use_ms(receptor_ms: Path, destination: Path) -> None:
    if not receptor_ms.exists():
        raise FileNotFoundError(
            "Arquivo .ms não encontrado. O sphgen precisa do arquivo de superfície molecular, "
            "como o rec.ms."
        )

    destination.parent.mkdir(parents=True, exist_ok=True)

    if receptor_ms.resolve() != destination.resolve():
        shutil.copyfile(receptor_ms, destination)


def copy_or_prepare_receptor_mol2(
    receptor_pdb: Path,
    destination: Path,
    obabel_exe: str,
    receptor_mol2_source: Path | None,
) -> None:
    if receptor_mol2_source is not None:
        if not receptor_mol2_source.exists():
            raise FileNotFoundError(f"Receptor MOL2 não encontrado: {receptor_mol2_source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        if receptor_mol2_source.resolve() != destination.resolve():
            shutil.copyfile(receptor_mol2_source, destination)
        return

    print(
        "Aviso: --receptor-mol2 não foi informado. "
        "Vou tentar gerar receptor_fixed.mol2 com Open Babel, mas isso pode falhar no DOCK6.",
        file=sys.stderr,
    )
    prepare_receptor_mol2(receptor_pdb, destination, obabel_exe)


def copy_or_prepare_ligands(
    ligand_source: Path,
    ligand_dock_destination: Path,
    ligand_ref_destination: Path,
    obabel_exe: str,
    ligand_dock_source: Path | None,
    ligand_ref_source: Path | None,
) -> None:
    """
    Para DOCK6, o ideal é separar:
    - ligand_ref.mol2: ligante cristalográfico usado para selecionar spheres.
    - ligand_dock.mol2: ligante preparado para o docking.

    Se os arquivos prontos forem informados, o script copia esses arquivos.
    Caso contrário, usa o ligante baixado e gera ambos por Open Babel.
    """
    ligand_dock_destination.parent.mkdir(parents=True, exist_ok=True)

    if ligand_dock_source is not None:
        if not ligand_dock_source.exists():
            raise FileNotFoundError(f"Ligante de docking MOL2 não encontrado: {ligand_dock_source}")
        if ligand_dock_source.resolve() != ligand_dock_destination.resolve():
            shutil.copyfile(ligand_dock_source, ligand_dock_destination)
    else:
        prepare_ligand_mol2(ligand_source, ligand_dock_destination, obabel_exe)

    if ligand_ref_source is not None:
        if not ligand_ref_source.exists():
            raise FileNotFoundError(f"Ligante de referência MOL2 não encontrado: {ligand_ref_source}")
        if ligand_ref_source.resolve() != ligand_ref_destination.resolve():
            shutil.copyfile(ligand_ref_source, ligand_ref_destination)
    else:
        # Fallback: usa o mesmo ligante preparado como referência.
        # Isso pode falhar quando o ligante baixado não está no sistema de coordenadas do cristal.
        if ligand_dock_destination.resolve() != ligand_ref_destination.resolve():
            shutil.copyfile(ligand_dock_destination, ligand_ref_destination)


def write_sphgen_input(paths: Dock6Paths) -> None:
    paths.sphgen_in.write_text(
        "\n".join(
            [
                paths.receptor_ms.name,
                "R",
                "X",
                "0.0",
                "4.0",
                "1.4",
                paths.all_spheres.name,
                "",
            ]
        ),
        encoding="utf-8",
    )


def run_sphgen(paths: Dock6Paths, sphgen_exe: str) -> None:
    write_sphgen_input(paths)

    run_command(
        [sphgen_exe],
        cwd=paths.run_dir,
        log_file=paths.run_dir / "sphgen.log",
    )

    if not paths.all_spheres.exists():
        candidates = sorted(paths.run_dir.glob("*.sph"))
        if candidates:
            shutil.copyfile(candidates[0], paths.all_spheres)

    if not paths.all_spheres.exists():
        raise FileNotFoundError(f"sphgen não gerou {paths.all_spheres}")


def run_sphere_selector(
    paths: Dock6Paths,
    ligand_ref_mol2: Path,
    sphere_selector_exe: str,
    radius: float,
) -> None:
    run_command(
        [
            sphere_selector_exe,
            paths.all_spheres.name,
            ligand_ref_mol2.name,
            str(radius),
        ],
        cwd=paths.run_dir,
        log_file=paths.run_dir / "sphere_selector.log",
    )

    if paths.selected_spheres.exists():
        return

    candidates = sorted(paths.run_dir.glob("*selected*.sph"))
    if candidates:
        shutil.copyfile(candidates[0], paths.selected_spheres)
        return

    sph_candidates = sorted(paths.run_dir.glob("*.sph"))
    for candidate in sph_candidates:
        if candidate.name != paths.all_spheres.name:
            shutil.copyfile(candidate, paths.selected_spheres)
            return

    raise FileNotFoundError(
        "sphere_selector não gerou selected_spheres.sph. "
        "Envie sphere_selector.log para ajuste."
    )


def write_showbox_input(paths: Dock6Paths) -> None:
    paths.showbox_in.write_text(
        "\n".join(
            [
                "Y",
                "5",
                paths.selected_spheres.name,
                "1",
                paths.box_pdb.name,
                "",
            ]
        ),
        encoding="utf-8",
    )


def run_showbox(paths: Dock6Paths, showbox_exe: str) -> None:
    write_showbox_input(paths)

    stdin_text = paths.showbox_in.read_text(encoding="utf-8")

    run_command(
        [showbox_exe],
        cwd=paths.run_dir,
        stdin_text=stdin_text,
        log_file=paths.run_dir / "showbox.log",
    )

    if not paths.box_pdb.exists():
        candidates = sorted(paths.run_dir.glob("*box*.pdb"))
        if candidates:
            shutil.copyfile(candidates[0], paths.box_pdb)

    if not paths.box_pdb.exists():
        raise FileNotFoundError(
            "showbox não gerou rec_box.pdb. Envie showbox.log para ajuste."
        )


def write_grid_input(
    paths: Dock6Paths,
    dock6_param_dir: Path,
    grid_spacing: float = 0.3,
) -> None:
    vdw_defn_file = dock6_param_dir / "vdw_AMBER_parm99.defn"

    if not vdw_defn_file.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {vdw_defn_file}")

    paths.grid_in.write_text(
        "\n".join(
            [
                "compute_grids yes",
                f"grid_spacing {grid_spacing}",
                "output_molecule no",
                "contact_score no",
                "energy_score yes",
                "energy_cutoff_distance 10",
                "atom_model a",
                "attractive_exponent 6",
                "repulsive_exponent 12",
                "distance_dielectric yes",
                "dielectric_factor 4",
                "allow_non_integral_charges yes",
                "bump_filter yes",
                "bump_overlap 0.75",
                "",
                f"receptor_file {paths.receptor_mol2.name}",
                f"box_file {paths.box_pdb.name}",
                f"vdw_definition_file {vdw_defn_file}",
                f"score_grid_prefix {paths.grid_prefix.name}",
                "",
            ]
        ),
        encoding="utf-8",
    )


def run_grid(paths: Dock6Paths, grid_exe: str) -> None:
    run_command(
        [grid_exe, "-i", paths.grid_in.name, "-o", "grid.out"],
        cwd=paths.run_dir,
        log_file=paths.run_dir / "grid_command.log",
    )

    # No fluxo usado no material da disciplina, as saídas essenciais do grid são
    # grid.bmp e grid.nrg. Algumas versões/configurações também podem gerar
    # grid.cnt, mas ele não é obrigatório quando contact_score está desligado.
    required_extensions = [".bmp", ".nrg"]
    missing_required = [
        ext for ext in required_extensions if not Path(str(paths.grid_prefix) + ext).exists()
    ]

    if missing_required:
        raise FileNotFoundError(
            f"grid terminou, mas faltam arquivos essenciais do grid: {missing_required}. "
            f"Verifique {paths.run_dir / 'grid.out'}"
        )

    optional_cnt = Path(str(paths.grid_prefix) + ".cnt")
    if not optional_cnt.exists():
        print(
            "Aviso: grid.cnt não foi gerado. Isso pode ser normal nesta configuração "
            "quando contact_score está desligado. Continuando com grid.bmp e grid.nrg.",
            file=sys.stderr,
        )


def write_dock6_input(
    paths: Dock6Paths,
    dock6_param_dir: Path,
    max_orientations: int = 500,
    ligand_outfile_prefix: str = "dock6_run",
) -> None:
    vdw_defn_file = dock6_param_dir / "vdw_AMBER_parm99.defn"
    flex_defn_file = dock6_param_dir / "flex.defn"
    flex_drive_file = dock6_param_dir / "flex_drive.tbl"

    for required_file in [vdw_defn_file, flex_defn_file, flex_drive_file]:
        if not required_file.exists():
            raise FileNotFoundError(f"Arquivo não encontrado: {required_file}")

    paths.dock_in.write_text(
        "\n".join(
            [
                "conformer_search_type                                        flex",
                "write_fragment_libraries                                     no",
                "user_specified_anchor                                        no",
                "limit_max_anchors                                            no",
                "min_anchor_size                                              5",
                "pruning_use_clustering                                       yes",
                "pruning_max_orients                                          1000",
                "pruning_clustering_cutoff                                    100",
                "pruning_orient_score_cutoff                                  1000.0",
                "pruning_conformer_score_cutoff                               100.0",
                "pruning_conformer_score_scaling_factor                       1.0",
                "use_clash_overlap                                            no",
                "write_growth_tree                                            no",
                "use_internal_energy                                          yes",
                "internal_energy_rep_exp                                      12",
                "internal_energy_cutoff                                       100.0",
                f"ligand_atom_file                                             {paths.ligand_mol2.name}",
                "limit_max_ligands                                            no",
                "skip_molecule                                                no",
                "read_mol_solvation                                           no",
                "calculate_rmsd                                               no",
                "use_database_filter                                          no",
                "orient_ligand                                                yes",
                "automated_matching                                           yes",
                f"receptor_site_file                                           {paths.selected_spheres.name}",
                f"max_orientations                                             {max_orientations}",
                "critical_points                                              no",
                "chemical_matching                                            no",
                "use_ligand_spheres                                           no",
                "bump_filter                                                  no",
                "score_molecules                                              yes",
                "contact_score_primary                                        no",
                "grid_score_primary                                           yes",
                "grid_score_rep_rad_scale                                     1",
                "grid_score_vdw_scale                                         1",
                "grid_score_es_scale                                          1",
                "grid_lig_efficiency                                          no",
                f"grid_score_grid_prefix                                       {paths.grid_prefix.name}",
                "minimize_ligand                                              yes",
                "minimize_anchor                                              yes",
                "minimize_flexible_growth                                     yes",
                "use_advanced_simplex_parameters                              no",
                "minimize_flexible_growth_ramp                                yes",
                "simplex_max_cycles                                           1",
                "simplex_score_converge                                       0.1",
                "simplex_initial_score_coverge                                5",
                "simplex_cycle_converge                                       1.0",
                "simplex_trans_step                                           1.0",
                "simplex_rot_step                                             0.1",
                "simplex_tors_step                                            10.0",
                "simplex_anchor_max_iterations                                500",
                "simplex_grow_max_iterations                                  250",
                "simplex_grow_tors_premin_iterations                          0",
                "simplex_final_min                                            no",
                "simplex_random_seed                                          0",
                "simplex_restraint_min                                        no",
                "atom_model                                                   a",
                f"vdw_defn_file                                                {vdw_defn_file}",
                f"flex_defn_file                                               {flex_defn_file}",
                f"flex_drive_file                                              {flex_drive_file}",
                f"ligand_outfile_prefix                                        {ligand_outfile_prefix}",
                "write_mol_solvation                                          no",
                "write_orientations                                           no",
                "num_final_scored_poses                                       5",
                "num_preclustered_conformers                                  5",
                "write_conformations                                          no",
                "cluster_conformations                                        yes",
                "cluster_rmsd_threshold                                       2.0",
                "score_threshold                                              100.0",
                "rank_ligands                                                 no",
                "",
            ]
        ),
        encoding="utf-8",
    )


def run_dock6(paths: Dock6Paths, dock6_exe: str) -> None:
    run_command(
        [dock6_exe, "-i", paths.dock_in.name, "-o", paths.dock_out.name],
        cwd=paths.run_dir,
        log_file=paths.run_dir / "dock6_command.log",
    )

    candidates = [
        paths.run_dir / "dock6_run_scored.mol2",
        paths.run_dir / "scored.mol2",
        paths.run_dir / "scored_scored.mol2",
    ]

    for candidate in candidates:
        if candidate.exists():
            if candidate != paths.scored_mol2:
                shutil.copyfile(candidate, paths.scored_mol2)
            return

    mol2_candidates = sorted(paths.run_dir.glob("*.mol2"))
    for candidate in mol2_candidates:
        if "score" in candidate.name.lower():
            shutil.copyfile(candidate, paths.scored_mol2)
            return

    raise FileNotFoundError(
        "DOCK6 terminou, mas não encontrei o arquivo MOL2 de saída. "
        "Envie dock6.out e a lista de arquivos da pasta."
    )


def run_dock6_pipeline_for_ligand(
    ligand_label: str,
    ligand_source: Path,
    receptor_pdb: Path,
    receptor_ms: Path,
    base_output_dir: Path,
    dock6_param_dir: Path,
    receptor_mol2_source: Path | None = None,
    ligand_ref_mol2_source: Path | None = None,
    ligand_dock_mol2_source: Path | None = None,
    obabel_exe: str = "obabel",
    sphgen_exe: str = "sphgen",
    sphere_selector_exe: str = "sphere_selector",
    showbox_exe: str = "showbox",
    grid_exe: str = "grid",
    dock6_exe: str = "dock6",
    sphere_radius: float = 10.0,
    prepare_only: bool = False,
) -> Path:
    paths = make_paths(base_output_dir, ligand_label)
    paths.run_dir.mkdir(parents=True, exist_ok=True)

    clean_previous_run(paths)

    copy_or_prepare_receptor_mol2(
        receptor_pdb=receptor_pdb,
        destination=paths.receptor_mol2,
        obabel_exe=obabel_exe,
        receptor_mol2_source=receptor_mol2_source,
    )
    ligand_ref_mol2 = paths.run_dir / "ligand_ref.mol2"
    copy_or_prepare_ligands(
        ligand_source=ligand_source,
        ligand_dock_destination=paths.ligand_mol2,
        ligand_ref_destination=ligand_ref_mol2,
        obabel_exe=obabel_exe,
        ligand_dock_source=ligand_dock_mol2_source,
        ligand_ref_source=ligand_ref_mol2_source,
    )
    copy_or_use_ms(receptor_ms, paths.receptor_ms)

    run_sphgen(paths, sphgen_exe)
    run_sphere_selector(paths, ligand_ref_mol2, sphere_selector_exe, sphere_radius)
    run_showbox(paths, showbox_exe)

    write_grid_input(paths, dock6_param_dir=dock6_param_dir)
    write_dock6_input(paths, dock6_param_dir=dock6_param_dir)

    if prepare_only:
        return paths.run_dir

    run_grid(paths, grid_exe)
    run_dock6(paths, dock6_exe)

    return paths.run_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Executa docking com DOCK6.")

    parser.add_argument("--receptor-pdb", required=True)
    parser.add_argument("--receptor-ms", required=True)
    parser.add_argument("--receptor-mol2")
    parser.add_argument("--ligand-ref-mol2")
    parser.add_argument("--ligand-dock-mol2")
    parser.add_argument("--ligand-file", required=True)
    parser.add_argument("--ligand-label", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--dock6-param-dir", required=True)

    parser.add_argument("--obabel-exe", default="obabel")
    parser.add_argument("--sphgen-exe", default="sphgen")
    parser.add_argument("--sphere-selector-exe", default="sphere_selector")
    parser.add_argument("--showbox-exe", default="showbox")
    parser.add_argument("--grid-exe", default="grid")
    parser.add_argument("--dock6-exe", default="dock6")

    parser.add_argument("--sphere-radius", type=float, default=10.0)
    parser.add_argument("--prepare-only", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        for exe in [
            args.obabel_exe,
            args.sphgen_exe,
            args.sphere_selector_exe,
            args.showbox_exe,
        ]:
            ensure_executable(exe)

        if not args.prepare_only:
            ensure_executable(args.grid_exe)
            ensure_executable(args.dock6_exe)

        run_dir = run_dock6_pipeline_for_ligand(
            ligand_label=args.ligand_label,
            ligand_source=Path(args.ligand_file),
            receptor_pdb=Path(args.receptor_pdb),
            receptor_ms=Path(args.receptor_ms),
            receptor_mol2_source=Path(args.receptor_mol2) if args.receptor_mol2 else None,
            ligand_ref_mol2_source=Path(args.ligand_ref_mol2) if args.ligand_ref_mol2 else None,
            ligand_dock_mol2_source=Path(args.ligand_dock_mol2) if args.ligand_dock_mol2 else None,
            base_output_dir=Path(args.output_dir),
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

        print(f"Pipeline DOCK6 concluído em: {run_dir}")
        return 0

    except Exception as exc:
        print(f"Erro: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())