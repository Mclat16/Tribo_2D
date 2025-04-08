from tribo_2D.afm import *
from tribo_2D.tools import *


# run=afm('afm_config.ini')
# run.system()
# run.load()
# run.slide()
# run.pbs()





from tribo_2D.sheet import *

run=sheetvsheet('sheet_config.ini')
run.system()
run.pbs()

for file in os.listdir():
    if file.endswith(".cif") or file.endswith(".lmp") or file.endswith(".json"):
        os.remove(file)