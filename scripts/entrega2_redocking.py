"""
Script automatizado para redocking + RMSD + análise.
Executa:
1) Download da estrutura
2) Preparação da proteína
3) Download do ligante cristalográfico
4) Docking do ligante
5) Cálculo de RMSD
6) Resumo dos resultados
Exemplo: python3 scripts/entrega2_redocking.py 9THJ A A1JV9 resultados/entrega2
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
# Import dos módulos locais
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
import base
import docking
import analysis

def main(
    pdb_id: str,
    chain_id: str,
    ligand_id: str,
    output_dir: str,
    tight_box: bool = False,
    exhaustiveness: int = 12,
    num_modes: int = 15,
    energy_range: int = 4,
) -> int:
    output_dir = Path(output_dir)
    print(f"\n{'-x-'*20}")
    print(f"Pipeline da 2ª Entrega: Redocking + RMSD")
    print(f"{'-x-'*20}\n")
    try:
        print(f"[1/6] Baixando estrutura {pdb_id}")
        raw_structure = base.download_structure(pdb_id, output_dir)
        print(f"Estrutura salva em: {raw_structure}\n")
        print(f"[2/6] Preparando proteína (cadeia {chain_id})")
        if raw_structure.suffix.lower() == ".cif":
            receptor_file = base.extract_chain_from_mmcif(
                raw_structure, chain_id, output_dir, pdb_id
            )
        else:
            receptor_file = base.extract_chain_from_pdb(
                raw_structure, chain_id, output_dir, pdb_id
            )
        print(f"Proteína salva em: {receptor_file}\n")
        print(f"[3/6] Baixando ligante cristalográfico ({ligand_id})...")
        ligand_file = base.download_ligand(ligand_id, output_dir)
        print(f"Ligante salvo em: {ligand_file}\n")
        print(f"[4/6] Executando docking")
        docking_dir = output_dir / "docking" / ligand_id
        ligand_pdbqt = docking_dir / "ligand.pdbqt"
        receptor_pdbqt = docking_dir / "receptor.pdbqt"
        docking.prepare_receptor_pdbqt(receptor_file, receptor_pdbqt, "obabel")
        docking.prepare_ligand_pdbqt(ligand_file, ligand_pdbqt, "obabel")
        print(f"Arquivos PDBQT preparados")
        padding = 4.0 if tight_box else 8.0
        box = docking.compute_box(receptor_file, padding=padding)
        print(f"Box calculado: center=({box.center_x:.2f}, {box.center_y:.2f}, {box.center_z:.2f})")
        print(f"Modo tight-box: {'SIM (padding=4.0)' if tight_box else 'NÃO (padding=8.0)'}")
        print(f"Exaustão: {exhaustiveness}, Modos: {num_modes}, Intervalo energético: {energy_range}")
        config_file = docking_dir / "config.txt"
        output_pdbqt = docking_dir / "out.pdbqt"
        docking.write_vina_config(
            config_file, receptor_pdbqt, ligand_pdbqt, box,
            exhaustiveness=exhaustiveness, num_modes=num_modes, energy_range=energy_range, output_pdbqt=output_pdbqt
        )
        import subprocess
        completed = subprocess.run(
            ["vina", "--config", str(config_file)],
            check=False, text=True, capture_output=True
        )
        if completed.returncode != 0:
            print(f"X - Vina falhou: {completed.stderr}")
            return 1
        print(f"Docking concluído")
        log_file = docking_dir / "log.txt"
        log_file.write_text(
            "\n".join(filter(None, [completed.stdout.strip(), completed.stderr.strip()])),
            encoding="utf-8",
        )
        print(f"Log salvo em: {log_file}\n")
        print(f"[5/6] Calculando RMSD")
        ligand_pdb = output_dir / "ligands" / f"{ligand_id}_crystal.pdb"
        ligand_pdb.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["obabel", "-isdf", str(ligand_file), "-opdb", "-O", str(ligand_pdb)],
            check=False, capture_output=True
        )
        result = analysis.analyze_docking_result(docking_dir, ligand_pdb)
        score = result.get("best_score")
        rmsd = result.get("rmsd")
        if score is not None:
            print(f"Score (Vina): {score:.2f} kcal/mol")
        else:
            print(f"X - Score não encontrado no log")
        if rmsd is not None:
            print(f"RMSD: {rmsd:.2f} Å")
            if rmsd <= 2.0:
                print(f"Redocking bem-sucedido! (RMSD <= 2.0 Å)")
            else:
                print(f"RMSD elevado. Verificar pose manualmente.")
        else:
            print(f"RMSD não pôde ser calculado")
        print()
        print(f"[6/6] Gerando resumo")
        results_json = output_dir / "resultados_redocking.json"
        import json
        json.dump([result], results_json.open("w"), indent=2)
        print(f"Resultados salvos em: {results_json}\n")
        print(f"{'-x-'*20}")
        print(f"Pipeline completo! Arquivos gerados em: {output_dir}")
        print(f"{'-x-'*20}\n")
        # Resumo visual
        print("\nRESUMO:")
        print(f"-> PDB ID: {pdb_id}")
        print(f"-> Cadeia: {chain_id}")
        print(f"-> Ligante cristalográfico: {ligand_id}")
        print(f"-> Score: {score:.2f} kcal/mol" if score else "-> Score: N/A")
        print(f"-> RMSD: {rmsd:.2f} Å" if rmsd else "-> RMSD: N/A")
        print(f"-> Diretório de saída: {output_dir}")
        return 0
    except Exception as exc:
        print(f"\nX - Erro: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Executa pipeline da 2ª entrega: redocking + análise."
    )
    parser.add_argument("pdb_id", help="Identificador PDB (ex: 9THJ)")
    parser.add_argument("chain", help="Cadeia da proteína (ex: A)")
    parser.add_argument("ligand_id", help="ID do ligante cristalográfico (ex: A1JV9)")
    parser.add_argument("output_dir", help="Diretório de saída")
    parser.add_argument("--tight-box", action="store_true", help="Usar caixa de docking mais restrita (padding=4.0)")
    parser.add_argument("--exhaustiveness", type=int, default=12, help="Exaustão de busca do Vina (padrão: 12)")
    parser.add_argument("--num-modes", type=int, default=15, help="Número de modos retornados (padrão: 15)")
    parser.add_argument("--energy-range", type=int, default=4, help="Intervalo energético em kcal/mol (padrão: 4)")
    args = parser.parse_args()
    raise SystemExit(main(
        args.pdb_id, args.chain, args.ligand_id, args.output_dir,
        tight_box=args.tight_box,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
        energy_range=args.energy_range
    ))
