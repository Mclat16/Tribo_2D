from .settings import *
from mp_api.client import MPRester
import subprocess
import os
import json
from .Potentials.lj_params import lj_params
from mpi4py import MPI
from lammps import lammps
import numpy as np
from ase import io, data
from pathlib import Path
from CifFile import ReadCif
import configparser
import argparse
import re
from collections import defaultdict

def matsearch(material_id):
    material_id = str(material_id)
    with MPRester("ASGI0EvO83K5vj5GFmdJCYOpd7qgVTAL") as mpr:
        data = mpr.materials.summary.search(material_ids=[material_id], fields=["structure", "formula_pretty","composition_reduced",'nsites','symmetry'])
        formula = data[0].formula_pretty
        comp = {str(element): int(count) for element, count in data[0].composition_reduced.items()}
        elements = []
        for element, count in comp.items():
            elements.extend([element] * count)
        nelements = len(elements)
        structure = data[0].structure
        nsites = data[0].nsites
        symmetry = data[0].symmetry
        filename = formula + ".cif"
        structure.to(filename)

        return structure,formula,nelements,elements
    
def cifread(cif):

    with open(cif, 'r') as f:
        lines = f.readlines()
    
    lat_a,lat_b,lat_c = [None]*3
    elements = []
    reading_elements = False  # Flag to start reading elements
    header_skipped =  False
    formula = None
    filename = os.path.basename(cif).split('.')[0]

    for line in lines:
        if line.startswith('_cell_length_a'):
            lat_a = float(line.split()[1])

        elif line.startswith('_cell_length_b'):
            lat_b = float(line.split()[1])

        elif line.startswith('_cell_length_c'):
            lat_c = float(line.split()[1])

        elif line.startswith('_chemical_formula_structural'):
            formula = line.split(maxsplit=1)[1].strip()

        elif line.startswith('_cell_angle_alpha'):
            ang_a = float(line.split()[1])

        elif line.startswith('_cell_angle_beta'):
            ang_b = float(line.split()[1])
        
        elif line.startswith('_cell_angle_gamma'):
            ang_g = float(line.split()[1])

            # Detect when element section starts
        elif line.strip().startswith('loop_'):
            reading_elements = True  # Start reading element symbols
            continue
            # If reading elements, first skip headers
        elif reading_elements == True and not header_skipped== True:
            if line.strip().startswith('_'):  # These are column headers
                continue  # Skip header lines
            header_skipped = True  # Once we hit real data, stop skipping

        # Now read element names
        if reading_elements and header_skipped:
            if line.strip().startswith('_') or line.strip().startswith('loop_') or line.strip() == "":
                break  # Stop reading when another section starts
            element = line.strip().split()[0]  # First column is the element name
            elements.append(element)

    # Use filename if formula is missing
        if not formula:
            formula = filename

        # Convert formula (e.g., "Al2 O2") to dictionary {"Al": 2, "O": 2}
        elem_count = defaultdict(int)
        matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula)  # Match elements and their counts

        for element, count in matches:
            elem_count[element] = int(count) if count else 1  # Default count to 1 if missing


    nelements = len(elements)
    cif = {
        "lattice": (lat_a, lat_b, lat_c),
        "angles": (ang_a, ang_b, ang_g),
        "elem_count": dict(elem_count),
        "elements": elements,
        "formula": formula,
        "nelements": nelements,
        "filename": filename
    }

    return cif

def count_elemtypes(file):
    # Dictionary to store the highest number for each element
    elem_type = defaultdict(int)

    # Regular expression to match element symbols (with or without numbers)
    pattern = re.compile(r'([A-Za-z]+)(\d*)')  # Matches "Al1", "O2", or just "Al", "O"

    # Open the potential file
    with open(file, 'r') as f:
        lines = f.readlines()

    # Loop through the lines in the file
    for line in lines:
        clean_line = line.strip()

        # Skip comments and empty lines
        if clean_line.startswith('#') or not clean_line:
            continue
        
        # Split the line into components (by spaces or commas)
        parts = clean_line.split()

        # We want to count the first three elements (which are atomic species)
        if len(parts) >= 3:
            for element in parts[:3]:  # Check first three items in the line
                match = pattern.match(element)
                if match:
                    element_name = match.group(1)  # Extract the element (e.g., "Al", "O")
                    element_number = match.group(2)  # Extract the number if it exists

                    if element_number:  # If there's a number (e.g., "Al2")
                        element_number = int(element_number)
                    else:  # If there's no number (e.g., "Al")
                        element_number = 1  

                    # Update the highest number for this element
                    elem_type[element_name] = max(elem_type[element_name], element_number)

    return dict(elem_type)


def slab_generator(file,x,y,z):

    atomsk_command = f"atomsk {file}.cif -orthogonal-cell cif -ow -v 0"
    subprocess.run(atomsk_command, shell=True, check=True)
    cif = readcif(file+".cif")
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
                
    return xlo, xhi, ylo, yhi, zlo, zhi

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

def center(system,filename,var,natype):

        lmp = lammps(cmdargs=["-log", "none", "-screen", "none",  "-nocite"])
        # lmp.file("lammps/in.init")
        lmp.commands_list([
        "units           metal\n",
        "atom_style      atomic\n",
        "neighbor        0.3 bin\n",
        "boundary        p p p",
        "neigh_modify    every 1 delay 0 check yes #every 5\n\n",
        f"region box block {var['dim']['xlo']} {var['dim']['xhi']} {var['dim']['ylo']} {var['dim']['yhi']} -50 50\n",
        f"create_box      {var['data'][system][2]} box\n\n",
        # "read_data       %s/system_build/%s_%d.lmp" % (self.directory_l,self.data["2D"][1],self.layers),
        "read_data       %s add append" % filename,
        ])
        
        elem = var['data'][system][3].copy()

        count = {}
        result=[]
        mass = []
        for i in range(len(elem)):
            mass=data.atomic_masses[data.atomic_numbers[elem[i]]]
            lmp.command("mass %d %d" % (i+1, mass))
            element = elem[i]
            count[element] = count.get(element, 0) + 1
            result.append(element + str(count[element]))
        for i in range(len(elem)):
            element = elem[i]
            if count[element] > 1:
                elem[i]=result[i]
        
        lmp.commands_list([
        "pair_style sw",
        "pair_coeff * * Potentials/%s.sw %s" % (var['data'][system][1], ' '.join(elem)),
        "compute zmin all reduce min z",
        "variable disp equal -c_zmin",
        "thermo_style  custom step temp pe ke etotal c_zmin v_disp",
        "run 0",
        "displace_atoms all move 0.0 0.0 v_disp units box",
        f"change_box all z final {var['dim']['zlo']} {var['dim']['zhi']}",
        "run 0",
        "write_data  %s" % filename
        ])
        lmp.close

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