from .settings import *
from mp_api.client import MPRester
import subprocess
import os
import json
from .Potentials.lj_params import lj_params
import numpy as np
import configparser
import re

def matsearch(material_id):
    material_id = str(material_id)
    with MPRester("ASGI0EvO83K5vj5GFmdJCYOpd7qgVTAL") as mpr:
        data = mpr.materials.summary.search(material_ids=[material_id], fields=["structure", "formula_pretty","composition_reduced"])
        formula = data[0].formula_pretty
        comp = {str(element): int(count) for element, count in data[0].composition_reduced.items()}
        elements = []
        for element, count in comp.items():
            elements.extend([element] * count)
        nelements = len(elements)
        structure = data[0].structure
        filename = formula + ".cif"
        structure.to(filename)

    cif={}

    cif= {
        'lat_a': structure.lattice.a,
        'lat_b': structure.lattice.b,
        'lat_c': structure.lattice.c,
        'formula': formula,
        'ang_a': structure.lattice.alpha,
        'ang_b': structure.lattice.beta,
        'ang_g': structure.lattice.gamma,
        "elem_comp": comp,
        "nelements": nelements,
        "filename": filename,
        "elements": elements
    }
    
    return cif
    
def cifread(cif):

    filename = os.path.basename(cif).split('.')[0]

    with open(cif, 'r') as f:
        lines = f.readlines()
    
    elements = []
    cif = {}
    reading_elements = False 
    header_skipped =  False

    keys = {
        '_cell_length_a': 'lat_a',
        '_cell_length_b': 'lat_b',
        '_cell_length_c': 'lat_c',
        '_chemical_formula_structural': 'formula',
        '_cell_angle_alpha': 'ang_a',
        '_cell_angle_beta': 'ang_b',
        '_cell_angle_gamma': 'ang_g'
    }

    for line in lines:
        for key, var in keys.items():
            if line.startswith(key):
                value = line.split(maxsplit=1)[1].strip()
                cif[var] = float(value) if 'formula' not in var else value
        print(cif)
        # Detect when element section starts
        if line.strip().startswith('loop_'):
            reading_elements = True  # Start reading element symbols
            continue
            # If reading elements, first skip headers
        elif reading_elements == True and not header_skipped== True:
            if line.strip().startswith('_'):  # These are column headers
                continue  # Skip header lines
            header_skipped = True  # Once we hit real data, stop skipping

        # Now read element names
        if reading_elements and header_skipped:
            element = line.strip().split()[0]  # First column is the element name
            elements.append(element)
            # Convert formula (e.g., "Al2 O2") to dictionary {"Al": 2, "O": 2}
    
    
    elem_count = {}
    matches = re.findall(r'([A-Z][a-z]*)(\d*)', cif['formula'])  # Match elements and their count
    for element, count in matches:
        elem_count[element] = int(count) if count else 1  # Default count to 1 if missing
    nelements = len(elements)

    cif.update = ({
        "elements" : elements,
        "elem_comp": elem_count,
        "nelements": nelements,
        "filename" : filename
    })

    return cif

def count_elemtypes(file):

    elem_type = {}

    matches = re.compile(r'([A-Za-z]+)(\d*)')  # Matches "Al1", "O2", or just "Al", "O"

    # Open the potential file
    with open(file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        stripped_line = line.strip()

        # Skip comments and empty lines
        if stripped_line.startswith('#') or not stripped_line:
            continue
        
        # Split the line into components (by spaces or commas)
        parts = stripped_line.split()

        # We want to count the first three elements (which are atomic species)
        if len(parts) >= 3:
            for element in parts[:3]:  # Check first three items in the line
                match = matches.match(element)
                if match:
                    element_name = match.group(1)  # Extract the element (e.g., "Al", "O")
                    element_number = match.group(2)  # Extract the number if it exists

                    if element_number:  # If there's a number (e.g., "Al2")
                        element_number = int(element_number)
                                                
                    else:  # If there's no number (e.g., "Al")
                        element_number = 1  

                    if element_name not in elem_type:
                        elem_type[element_name] = 0  
                    
                    
                    elem_type[element_name] = max(elem_type[element_name], element_number)
                    

    return elem_type

def slab_generator(file,x,y,z):

    atomsk_command = f"atomsk {file}.cif -duplicate 2 2 1 -orthogonal-cell cif -ow -v 0"
    subprocess.run(atomsk_command, shell=True, check=True)
    cif = cifread(file+".cif")
    x2 = round(x/cif["lattice"][0])
    y2 = round(y/cif["lattice"][1])
    z2 = round(z/cif["lattice"][2])
    atomsk_command2 = f"atomsk {file}.cif -duplicate {x2} {y2} {z2} lammps -ow -v 0"
    subprocess.run(atomsk_command2, shell=True, check=True)

def get_model_dimensions(lmp):
    xlo , xhi, ylo, yhi, zlo, zhi = [None] *6
    with open(lmp, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "xlo xhi" in line:
                xlo, xhi = map(float, line.split()[0:2])
            elif "ylo yhi" in line:
                ylo, yhi = map(float, line.split()[0:2])
            elif "zlo zhi" in line:
                zlo, zhi = map(float, line.split()[0:2])
    dim = {
        "xlo": xlo,
        "xhi": xhi, 
        "ylo": ylo, 
        "yhi": yhi, 
        "zlo": zlo,
        "zhi": zhi
    }            
    return dim

def LJparams(X,Y):
    # Lorentz-Bertholt Mixing rules for obtaining the LJ parameters between two atoms using the Universal Force Field parameters
    e1 = lj_params[X][1]
    e2 = lj_params[Y][1]
    s1 = lj_params[X][0]
    s2 = lj_params[Y][0]
    epsilon = np.sqrt(e1*e2)
    sigma = (s1 + s2)/2
    return epsilon, sigma

def removeInlineComments(config):

    for section in config.sections():
        for item in config.items(section):
            config.set(section, item[0], item[1].split("#")[0].strip())
    return config

def read_config(input):
    config = configparser.ConfigParser()
    config.read(input)

    config = removeInlineComments(config)

    var = {}
    for section in config.sections():
        var[section] = {}
        for key in config[section]:
            # Assign variable
            value = config.get(section,key)

            if value.endswith(']'):
                var[section][key]= json.loads(value)
            elif value.isdigit():
                var[section][key] = int(value)
            elif '.' in value and value.replace('.', '', 1).isdigit():
                var[section][key] = float(value)
            elif value.replace('.', '', 1).replace('e', '', 1).replace('-', '', 1).isdigit():
                var[section][key] = float(value)
            else:
                var[section][key] = value
    return var