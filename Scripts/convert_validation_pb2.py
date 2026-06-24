from pathlib import Path
import subprocess
import pandas as pd


csv_path = Path("data/validation_set_pb2_dude_simple.csv")
out_dir = Path("data/validation_pb2_mol2")
out_dir.mkdir(parents=True, exist_ok=True)

if not csv_path.exists():
    raise FileNotFoundError(f"Arquivo não encontrado: {csv_path}")

df = pd.read_csv(csv_path)

required_cols = {"ligand_label", "smiles", "is_active"}
missing = required_cols - set(df.columns)

if missing:
    raise ValueError(f"Colunas ausentes no CSV: {missing}")

manifest_rows = []
errors = []

for _, row in df.iterrows():
    raw_label = str(row["ligand_label"])

    label = (
        raw_label
        .replace("/", "_")
        .replace("\\", "_")
        .replace(" ", "_")
        .replace(":", "_")
    )

    smiles = str(row["smiles"]).strip()
    is_active = row["is_active"]

    smi_file = out_dir / f"{label}.smi"
    mol2_file = out_dir / f"{label}.mol2"

    smi_file.write_text(f"{smiles} {label}\n")

    cmd = [
        "obabel",
        str(smi_file),
        "-O",
        str(mol2_file),
        "--gen3d",
        "--partialcharge",
        "gasteiger",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    ok = (
        result.returncode == 0
        and mol2_file.exists()
        and mol2_file.stat().st_size > 0
    )

    if not ok:
        errors.append(
            {
                "ligand_label": label,
                "smiles": smiles,
                "is_active": is_active,
                "error": result.stderr.strip(),
            }
        )
        continue

    manifest_rows.append(
        {
            "ligand_label": label,
            "smiles": smiles,
            "is_active": is_active,
            "mol2_file": str(mol2_file),
        }
    )

manifest = pd.DataFrame(manifest_rows)
manifest.to_csv(out_dir / "validation_pb2_manifest.csv", index=False)

if errors:
    pd.DataFrame(errors).to_csv(
        out_dir / "validation_pb2_conversion_errors.csv",
        index=False,
    )

print("Conversão finalizada.")
print(f"CSV de entrada: {csv_path}")
print(f"Pasta de saída: {out_dir}")
print(f"Ligantes no CSV: {len(df)}")
print(f"MOL2 gerados com sucesso: {len(manifest_rows)}")
print(f"Erros de conversão: {len(errors)}")
print(f"Manifesto: {out_dir / 'validation_pb2_manifest.csv'}")

if errors:
    print(f"Erros salvos em: {out_dir / 'validation_pb2_conversion_errors.csv'}")