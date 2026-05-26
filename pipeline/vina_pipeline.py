from __future__ import annotations
import argparse
import sys
from pathlib import Path

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Pipeline que roda pré-processamento (src.base) e Vina para redocking")
    parser.add_argument("pdb_id")
    parser.add_argument("chain")
    parser.add_argument("ligand_id")
    parser.add_argument("output_dir")
    parser.add_argument("--vina-exe", default="vina")
    args = parser.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root))
    
    from src.vina import redocking as redock

    try:
        best_pose, best_rmsd = redock.preprocess_and_dock(
            args.pdb_id, args.chain, args.ligand_id, args.output_dir, args.vina_exe
        )
        print(f"Melhor pose: {best_pose}")
        print(f"RMSD da melhor pose: {best_rmsd:.3f} Å")
    except Exception as exc:
        import traceback
        print(f"Erro durante o pipeline: {exc}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
