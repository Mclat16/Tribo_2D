from .tools import *
from .settings import *
from .materials import *

from mp_api.client import MPRester
import sys
import subprocess
import os
import json
from mpi4py import MPI
from lammps import lammps
import numpy as np
from ase import io, data
from pathlib import Path
from CifFile import ReadCif
import configparser
import argparse

def tip(var):
        
    filename = var['data']["tip"][1] + ".lmp"
    
    x = 2*var['tip']['r']
    y = x
    z = round(var['tip']['r']/3)
    side = round(2*var['tip']['r']/2.5)
    
    if var['tip']['amorph'] == 'a':
        
        am_filename = Path(__file__).parent / f"materials/amor_{filename}"
            
        if os.path.exists(am_filename):
            print("File exists")
        else:
            slab_generator(var['data']['tip'][1], 200, 200, 20)
            amorph(filename, 2500,var)
            os.remove(filename)
            print("File Removed")
        
        filename = am_filename
        
    else:
        slab_generator(var['data']['tip'][1],x,y,z)
    
    lmp = lammps(cmdargs=["-log", "none", "-screen", "none",  "-nocite"])
    lmp.commands_list([
    "boundary p p p",
    "units metal",
    "atom_style      atomic",
    "region          box block -%f %f -%f %f -3 %f units box" % (var['tip']['r'],var['tip']['r'],var['tip']['r'],var['tip']['r'],z),
    "create_box      %d box" % var['data']['tip'][2],
    "read_data       %s add append shift -%f -%f 0" % (filename,var['tip']['r'],var['tip']['r']),
    "region          afm_tip sphere 0 0 %f %f side in units box" % (var['tip']['r'],var['tip']['r']),
    "region tip intersect 2 afm_tip box",
    "group           tip region tip",
    "group           box subtract all tip",
    "delete_atoms    group box",
    "change_box all x final -%f %f y final -%f %f z final -3 %f " % (side,side,side,side,z),
    #------------------Save the final configuration to a data file
    "reset_atoms     id",
    "write_data      %s/system_build/tip.lmp" % var['dir']
            ])
    lmp.close

def sub(var):
    
    filename = var['data']['sub'][1] + ".lmp"
    
    if var['sub']['amorph'] == 'a':
        
        am_filename = Path(__file__).parent / f"materials/amor_{filename}"
            
        if os.path.exists(am_filename):
            print("File exists")
        else:
            slab_generator(var['data']["sub"][1], 200, 200, 20)
            amorph(filename, 2500,var)
            
            os.remove(filename)
            print("File Removed")
        
        filename = am_filename
    
    else:
        slab_generator(var['data']['sub'][1], var['dim']['xhi'], var['dim']['yhi'], 10)
    
    lmp = lammps(cmdargs=["-log", "none", "-screen", "none",  "-nocite"])
    lmp.commands_list([
        "boundary p p p",
        "units metal",
        "atom_style      atomic",
        "read_data %s" % filename,
        "region sub block 0 %f 0 %f 0 10 units box" % (var['dim']['xhi'], var['dim']['yhi']),
        "group sub region sub",
        "group box subtract all sub",
        "delete_atoms group box",
        "change_box all x final 0 %f y final 0 %f z final 0 10" % (var['dim']['xhi'], var['dim']['yhi']),
        "write_data      %s/system_build/sub.lmp" % var['dir']
    ])
    
    lmp.close

def amorph(filename,tempmelt,var):
    
    quench = (tempmelt-var['general']['temproom'])*100
    
    # lmp = lammps(cmdargs=["-log", "none", "-screen", "none",  "-nocite"])
    lmp = lammps()
    lmp.commands_list([
    "boundary p p p",
    "units metal",
    "atom_style      atomic",

    "read_data %s" % filename,
    #------------------Apply potentials-----------------------
    "pair_style sw",
    "pair_coeff * * tribo_2D/Potentials/Si.sw Si",

    ##########################################################
    #-------------------Energy Minimization------------------#
    ##########################################################
    "min_style       cg",
    "minimize        1.0e-4 1.0e-8 100 1000",
    ##########################################################
    #-------------------Equilibrate System-------------------#
    ##########################################################
    "timestep        0.001",
    "thermo          100",
    "thermo_style    custom step temp pe ke etotal press",
    # Specify melting temperature and timestep
    "velocity        all create %f 1234579 rot yes dist gaussian" % tempmelt,
    "run             0",
    # Equilibration at temperature
    "fix             melt all nvt temp %f %f $(100.0*dt)" % (tempmelt,tempmelt),
    "run             5000",
    "unfix           melt",
    
    ##########################################################
    #----------------------Quench System---------------------#
    ##########################################################
    "fix             quench all nvt temp %f %f $(100.0*dt) " % (tempmelt, var['general']['temproom']),
    "run             %d " % quench,
    "unfix           quench ",
    "write_data materials/amor_%s" % filename
    ])
    
    lmp.close

def sheet(var):

    x=var['2D'][1]
    y=var['2D'][2]

    filename = var['2D'][0]

    cif = cifread(filename +".cif")
    pot = count_elemtypes(f"{filename.lower()}/{filename.lower()}.sw")

    multiples = {}

    for element, cif_count in cif["elem_count"].items():
        potential_count = pot.get(element, 0)  # Get count from potential file or default to 0
    if cif_count > 0 and potential_count != 1:  # Avoid division by zero
        multiples[element] = potential_count / cif_count
    else:
        multiples[element] = None  # If the CIF count is zero (which shouldn't happen), mark as None
    
    first_multiple = next(iter(multiples.values()))  # Get the first multiple
    for multiple in multiples.values():
        if multiple != first_multiple:  
            raise ValueError("multiples must be the same")

    typecount = 0
    for atomcount in pot.values():
        typecount += atomcount

    print(first_multiple)     
    if os.path.exists("a.cif"):
        os.remove("a.cif")  # Delete the existing file
    if os.path.exists(f"a.lmp"):
        os.remove(f"a.lmp")  # Delete the existing file
    if os.path.exists(f"{filename}.lmp"):
        os.remove(f"{filename}.lmp")  # Delete the existing file
    if first_multiple == None:
        atomsk_command = f"atomsk {filename}.cif -duplicate 2 2 1 -ow a.lmp"
        subprocess.run(atomsk_command, shell=True, check=True)
    else:
        if first_multiple ==1:
            atomsk_command = f"atomsk {filename}.cif -ow a.lmp"
            subprocess.run(atomsk_command, shell=True, check=True)
        else:
            m = np.sqrt(first_multiple)
            atomsk_command = f"atomsk {filename}.cif -duplicate {int(m)} {int(m)} 1 a.lmp"
            subprocess.run(atomsk_command, shell=True, check=True)
            print("a.lmp created")

        with open("a.lmp", 'r') as file:
            lines = file.readlines()

        modified_lines = []
        atom_types_updated = False
        masses_section = False
        atoms_section = False
        atom_type_counter = 1  # To track incremental atom type assignment
        # Step 2: Create a dictionary of atom types from the 'Masses' section
        atom_types = {}
        elem2D = {}
        for i, line in enumerate(lines):
            stripped_line = line.strip()

            # Update the "atom types" value
            if re.match(r'^\s*\d+\s+atom types\s*$', stripped_line) and not atom_types_updated:
                lines[i]= f"   {typecount}  atom types\n"
                atom_types_updated = True
                continue

            if line.strip() == "Masses":
                masses_section = True
                continue  # Skip the "Masses" header line
            
            if masses_section:
                if "Atoms" in line:
                     # End of the Masses section
                    break
                
                parts = line.split()
                if len(parts) < 2:
                    continue  # Skip empty lines
                
                try:
                    atom_type_id = int(parts[0])  # First column is the atom type ID (1, 2, ...)
                    mass = float(parts[1])  # Second column is the mass
                    # Extract the atom type name from the comment section (after the `#`)
                    if "#" in line:
                        atom_type_name = line.split("#")[-1].strip()
                        lines[i]= ""
                    else:
                        atom_type_name = f"Unknown_{atom_type_id}"  # Fallback in case there's no comment
                    atom_types[atom_type_id] = (atom_type_name, mass)
                except ValueError:
                    continue  # Skip lines that cannot be converted properly

        modified_lines = set() 

        for i in range(1,len(atom_types)+1):
            atoms_section = False
            for l, line in enumerate(lines):
                stripped_line = line.strip()
                if "Atoms" in line:
                    atoms_section = True
                    continue

                if atoms_section and stripped_line and l not in modified_lines:
                        parts = stripped_line.split()

                        if parts[1] == str(i):
                            parts[1] = str(atom_type_counter)  # Update atom type
                            lines[l]= "  ".join(parts) + "\n"
                            modified_lines.add(l)
                            print(parts[1])
                            elem2D[atom_type_counter] = atom_types[i]
                            atom_type_counter += 1  # Increment for next line
                            continue

                        

        masses_section = False
        for i, line in enumerate(lines):
            print(elem2D)
            print("rewrite masses")
            stripped_line = line.strip()
            if line.strip() == "Masses":
                masses_section = True
                continue
            
            if masses_section:
                for l in range(1,len(elem2D)+1):
                    lines[i] += f"{l} {elem2D[l][1]}  #{elem2D[l][0]}\n"
                break

        # Save modified file
        with open("a.lmp", 'w') as file:
            file.writelines(lines)

        if first_multiple ==1:
            atomsk_command = f"atomsk a.lmp -duplicate 2 2 1 -ow lmp"
            subprocess.run(atomsk_command, shell=True, check=True)

        atomsk_command = f"atomsk a.lmp -orthogonal-cell -ow lmp"
        subprocess.run(atomsk_command, shell=True, check=True)


        dim = get_model_dimensions('a.lmp')
        print(dim)
        duplicate_a = round(x / dim['xhi'])
        duplicate_b = round(y / dim['yhi'])

        if os.path.exists(f"{filename}.lmp"):
            os.remove(f"{filename}.lmp")  # Delete the existing file
        atomsk_command = f"atomsk a.lmp -duplicate {duplicate_a} {duplicate_b} 1 {filename}.lmp"
        subprocess.run(atomsk_command, shell=True, check=True)