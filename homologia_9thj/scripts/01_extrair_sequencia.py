import os
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

# Get the directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)

pdb_file = os.path.join(parent_dir, "data/9thj.pdb")
chain_id = "A"
output_file = os.path.join(parent_dir, "data/9thj_A.fasta")

parser = PDBParser(QUIET=True)
structure = parser.get_structure("9thj", pdb_file)

sequence = ""

for model in structure:
    chain = model[chain_id]

    for residue in chain:
        if residue.id[0] == " ":
            resname = residue.get_resname()

            try:
                sequence += seq1(resname)
            except Exception:
                sequence += "X"

    break

with open(output_file, "w") as f:
    f.write(">9thj_A\n")

    for i in range(0, len(sequence), 60):
        f.write(sequence[i:i+60] + "\n")

print("Sequncia salva em:", output_file)
print("Tamanho:", len(sequence), "aminocidos")
