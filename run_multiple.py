from tribo_2D.afm import *
from tribo_2D.tools import *
from tribo_2D.sheet import *
import os

with open("material_list.txt", "r") as f:
    materials = [line.strip() for line in f]

# files = ["list_system", "list_load", "list_slide"]

# for f in files:
#     with open(Path("scripts/afm2/scripts/", f), "w"):
#         pass 

with open(f"afm_config_temp.ini", "r") as file:
    template_afm = file.read()

with open(f"sheet_config_temp.ini", "r") as file:
    template_sheet = file.read()

# materials = ["h-CrO2"]
for m in materials:
    updated_afm = template_afm.replace("{mat}", m)

    with open(f"afm_config.ini", "w") as file:
        file.write(updated_afm)

    updated_sheet = template_sheet.replace("{mat}", m)
    with open(f"sheet_config.ini", "w") as file:
        file.write(updated_sheet)

    run=afm('afm_config.ini')
    run.system()
    # run.slide()
    # run.pbs()

    # run=sheetvsheet('sheet_config.ini')
    # run.sheet_system()
    # run.sheet_pbs()

for file in os.listdir():
    if file.endswith(".cif") or file.endswith(".lmp") or file.endswith(".json"):
        os.remove(file)