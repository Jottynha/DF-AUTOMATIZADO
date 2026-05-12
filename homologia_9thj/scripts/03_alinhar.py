from modeller import *

env = Environ()
env.io.atom_files_directory = ["data"]

aln = Alignment(env)

mdl = Model(
    env,
    file="template_A_clean",
    model_segment=("FIRST:A", "LAST:A")
)

aln.append_model(
    mdl,
    align_codes="template",
    atom_files="template_A_clean.pdb"
)

aln.append(
    file="data/target.ali",
    align_codes="target",
    alignment_format="PIR"
)

aln.align2d(max_gap_length=50)

aln.write(
    file="data/target-template.ali",
    alignment_format="PIR"
)

aln.write(
    file="data/target-template.pap",
    alignment_format="PAP"
)

print("Alinhamento gerado.")
