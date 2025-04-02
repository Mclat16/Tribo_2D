from tribo_2D_2.afm import *

run=afm('afm_config.ini')
run.system()
# run.load()
# run.slide()

# for file in os.listdir():
#     if file.endswith(".cif") or file.endswith(".lmp") or file.endswith(".json"):
#         os.remove(file)



from tribo_2D_2.sheet import *

run=sheet('sheet_config.ini')
run.system()