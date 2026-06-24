import pandas as pd
from pathlib import Path

dock_path = Path("resultados/dock6_validation_pb2/dock6_results_summary.csv")
manifest_path = Path("data/validation_pb2_mol2/validation_pb2_manifest.csv")
out = Path("resultados/dock6_validation_pb2/dock6_validation_ranked_with_activity.csv")

dock = pd.read_csv(dock_path)
manifest = pd.read_csv(manifest_path)

dock["grid_score"] = pd.to_numeric(dock["grid_score"], errors="coerce")
dock["rmsd"] = pd.to_numeric(dock["rmsd"], errors="coerce")

valid = dock[dock["grid_score"].notna()].copy()

merged = valid.merge(
    manifest[["ligand_label", "smiles", "is_active"]],
    on="ligand_label",
    how="left"
)

merged = merged.sort_values("grid_score", ascending=True)
merged["rank_dock6"] = range(1, len(merged) + 1)

cols = [
    "rank_dock6",
    "ligand_label",
    "is_active",
    "smiles",
    "grid_score",
    "rmsd",
    "status",
    "scored_mol2",
    "dock_out",
    "run_dir",
    "ligand_file",
    "error",
]

merged = merged[[c for c in cols if c in merged.columns]]
merged.to_csv(out, index=False)

print(f"Arquivo salvo em: {out}")
print(f"Total original: {len(dock)}")
print(f"Resultados válidos: {len(merged)}")
print(f"Ativos: {merged['is_active'].sum() if 'is_active' in merged.columns else 'NA'}")
print()
print(merged.head(20).to_string(index=False))

path = "resultados/dock6_validation_pb2/dock6_validation_ranked_with_activity.csv"
df = pd.read_csv(path)

df["is_active"] = pd.to_numeric(df["is_active"], errors="coerce").fillna(0).astype(int)

total = len(df)
total_actives = df["is_active"].sum()
active_rate = total_actives / total if total > 0 else 0

print(f"Total de ligantes válidos: {total}")
print(f"Total de ativos: {total_actives}")
print(f"Taxa global de ativos: {active_rate:.4f}")
print()

for n in [10, 20, 50, 100]:
    top = df.head(n)
    actives_top = top["is_active"].sum()
    top_rate = actives_top / n
    ef = top_rate / active_rate if active_rate > 0 else 0

    print(f"Top {n}:")
    print(f"  Ativos no top {n}: {actives_top}")
    print(f"  Taxa de ativos no top {n}: {top_rate:.4f}")
    print(f"  Enrichment Factor EF{n}: {ef:.2f}")
    print()
    
path = "resultados/dock6_validation_pb2/dock6_validation_ranked_with_activity.csv"
df = pd.read_csv(path)

actives = df[df["is_active"] == 1].copy()

cols = ["rank_dock6", "ligand_label", "is_active", "grid_score", "rmsd"]
print(actives[cols].to_string(index=False))