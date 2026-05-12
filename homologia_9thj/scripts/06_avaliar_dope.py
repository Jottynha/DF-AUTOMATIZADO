from modeller import *
from modeller.scripts import complete_pdb

log.verbose()
env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

modelo = "resultados/modelos/target.B99990004.pdb"
mdl = complete_pdb(env, modelo)
s = Selection(mdl)
s.assess_dope(
    output='ENERGY_PROFILE NO_REPORT',
    file='resultados/modelo_dope.profile',
    normalize_profile=True,
    smoothing_window=15
)
print("Perfil DOPE salvo em resultados/modelo_dope.profile")

