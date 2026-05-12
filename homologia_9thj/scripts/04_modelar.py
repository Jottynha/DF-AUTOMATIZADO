from modeller import *
from modeller.automodel import *

log.verbose()
env = Environ()
env.io.atom_files_directory = ["data"]

a = AutoModel(
    env,
    alnfile="data/target-template.ali",
    knowns="template",
    sequence="target",
    assess_methods=(assess.DOPE, assess.GA341)
)
a.starting_model = 1
a.ending_model = 10
a.make()

