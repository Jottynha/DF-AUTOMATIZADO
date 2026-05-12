input_fasta = "data/2ynd_A.fasta"
output_pir = "data/target.ali"

sequence = ""
with open(input_fasta) as f:
    for line in f:
        line = line.strip()
        if not line.startswith(">"):
            sequence += line

with open(output_pir, "w") as f:
    f.write(">P1;target\n")
    f.write("sequence:target:::::::0.00: 0.00\n")
    for i in range(0, len(sequence), 60):
        f.write(sequence[i:i+60] + "\n")
    f.write("*\n")

print("Arquivo PIR criado:", output_pir)
