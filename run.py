from tribo_2D.afm import *

run=afm('afm_config.ini')
run.system()
run.load()
run.slide()
run.pbs()

# for file in os.listdir():
#     if file.endswith(".cif") or file.endswith(".lmp") or file.endswith(".json"):
#         os.remove(file)



# from tribo_2D.sheet import *

# run=sheet('sheet_config.ini')
# run.system()