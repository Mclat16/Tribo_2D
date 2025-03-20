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
    for atomcount in pot.values:
        typecount += atomcount
    
          
    if os.path.exists("a.cif"):
        os.remove("a.cif")  # Delete the existing file

    if os.path.exists(f"{filename}.lmp"):
        os.remove(f"{filename}.lmp")  # Delete the existing file

    if first_multiple == None:
        atomsk_command = f"atomsk {filename}.cif -duplicate 2 2 1 -orthogonal-cell  -ow a.lmp"
        subprocess.run(atomsk_command, shell=True, check=True)

    else:
        if first_multiple ==1:
            atomsk_command = f"atomsk {filename}.cif -ow a.lmp"
            subprocess.run(atomsk_command, shell=True, check=True)
        else:
            m = np.sqrt(first_multiple)
            atomsk_command = f"atomsk a.lmp -duplicate {m} {m} 1 -orthogonal-cell -ow a.lmp"
            subprocess.run(atomsk_command, shell=True, check=True)

        with open("a.lmp", 'r') as file:
            lines = file.readlines()
    
            # Step 1: Update the atom types line to match pot_count
            atom_types_line_index = None
            for i, line in enumerate(lines):
                if line.strip().startswith("atom types"):
                    atom_types_line_index = i
                    lines[i] = f"  {typecount}  atom types\n"  # Update atom types line

            # Step 2: Create a dictionary of atom types from the 'Masses' section
            atom_types = {}
            masses_section_started = False
            for i, line in enumerate(lines):
                if line.strip() == "Masses":
                    masses_section_started = True
                    continue  # Skip the "Masses" header line
                if masses_section_started:
                    if line.strip() == "Atoms # atomic":  # End of the Masses section
                        break
                    parts = line.split()
                    atom_type_id = parts[0]  # First column is the atom type ID (1, 2, ...)
                    mass = float(parts[1])  # Second column is the mass
                    # Extract the atom type name from the comment section (after the `#`)
                    atom_type_name = line.split("#")[-1].strip()
                    atom_types[int(atom_type_id)] = (atom_type_name, mass)
            
            
            # Step 3: Modify the second column in the Atoms section
            atoms_section_started = False
            for i, line in enumerate(lines):
                if line.strip() == "Atoms # atomic":
                    atoms_section_started = True
                    continue  # Skip the "Atoms # atomic" header line
                if atoms_section_started:
                    if line.strip() == "Masses":  # End of the Atoms section
                        break
                    parts = line.split()
                    parts[1] = str(int(parts[1]))  # Increment the atom type index in the second column
                    lines[i] = '  '.join(parts) + "\n"  # Join and update the line

            # Save the modified file
            with open("a.lmp", 'w') as file:
                file.writelines(lines)




                atomsk_command = f"atomsk a.lmp -duplicate 2 2 1 -orthogonal-cell -ow a.lmp"
                subprocess.run(atomsk_command, shell=True, check=True)




    xlo, xhi, ylo, yhi, zlo, zhi = get_model_dimensions(filename)
    duplicate_a = round(x / xhi)
    duplicate_b = round(y / yhi)

    atomsk_command2 = f"atomsk a.lmp -duplicate {duplicate_a} {duplicate_b} 1  -ow {filename}.lmp"
    os.remove("a.cif")

    # atomsk_command = f"atomsk {filename} -duplicate 2 2 1 -cut above 0.5*BOX z -cut below 0.2*BOX z -orthogonal-cell cif -ow"
    # atomsk_command = f"atomsk {filename} -duplicate 2 2 1 -orthogonal-cell cif -ow"


    subprocess.run(atomsk_command2, shell=True, check=True)