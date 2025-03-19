from .tools import *
from .settings import *
from .materials import *

from mp_api.client import MPRester
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