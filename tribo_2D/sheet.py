from .tools import *
from .build import *
from .settings.file import *
from .Potentials import *
from .materials.Build2D import *

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
import shutil


class sheet:
    def __init__(self,input):

        var = read_config(input) # Read input file settings

        var['data'] = {'2D': cifread(var['2D']['cif_path'])}  # Read materials  
        var['pot'] = {'2D': count_elemtypes(var['2D']['pot_path'])}  # Count potentials  
        
        var['data']['2D'].update({
            'natype':sum(var['pot']['2D'].values())
            })
        
        var['dir'] = f"scripts/sheetvsheet/{var['2D']['mat']}/size_{var['2D']['x']}x_{var['2D']['y']}y/K{self.var['general']['temproom']}"
        
        self.scripts = var['dir'] +"/scripts"
        dirs = ["data","lammps","visuals", "results", "system_build", "potentials", "scripts"]
        for d in dirs:
            Path(var['dir'], d).mkdir(parents=True, exist_ok=True)

        Path(self.scripts).mkdir(parents=True, exist_ok=True)
        with open(f"{self.scripts}/sheetvsheet", 'w'):
            pass

        self.var = sheet(var)
        self.var['dim'] = {}

        self.file = f"{self.directory}/system_build/{self.var['data']['2D'][1]}_4.lmp"
        self.var['ngroups'] = self.natype*4

    def system(self):
        lat_c = self.var['data']['2D'][0].lattice.c/2
        for force in self.var['general']['force']:
            if force == self.var['general']['scan_angle'][3]:
                scan_angle = np.arange(self.var['general']['scan_angle'][0],self.var['tip']['scan_angle'][1]+1,self.var['tip']['scan_angle'][2])
            else:
                scan_angle = [0]
                
            for a in scan_angle:
                filename = f"{self.directory}/lammps/{self.var['data']['2D'][1]}.lmp"
                with open(filename, 'w') as f:
                    init(f)
                    f.writelines([
                        "#------------------Create Geometry------------------------\n",
                        "#----------------- Define the simulation box -------------\n",
                        f"region          box block {self.var['dim']['xlo']} {self.var['dim']['xhi']} {self.var['dim']['ylo']} {self.var['dim']['yhi']} -40.0 40.0 units box\n",
                        f"create_box      {self.var['ngroups']} box bond/types 1 extra/bond/per/atom 100\n\n",

                        f"read_data       {self.file} add append\n\n",

                        "variable Cspring equal 5 # in eV/A^2\n\n ",


                    "#----------------- Create visualisation files ------------\n\n "
                    ])

                    if dump == True:
                        f.write(f"dump            sys all atom 100 ./{self.directory}/visuals/load_{force}N_l{layer}.lammpstrj\n\n",)

                    for t in range(self.natype):
                        t+=1
                        f.writelines([
                        f"group 2D_{t} type {t}\n",
                        ])

                    i = 0  
                     
                    for l in range(4):
                    
                        zlo= l*lat_c -1
                        zhi= zlo + lat_c
                        l+=1
                        f.writelines([
                        f"region layer_{l} block INF INF INF INF {zlo} {zhi} units box\n",
                        f"group layer_{l} region layer_{l} \n",
                        f"region layer_{l} delete\n",
                        ])

                        for t in range(self.natype):   
                            i+=1
                            f.writelines([
                            f"group layer intersect 2D_{t} layer_{l}\n",
                            f"set group layer type {i}\n",
                            "group layer delete\n\n"
                            ])
                            if l == 3:
                                f.write(f"group 2D_{t} delete\n\n")

                    settings_filename = self.directory+"/lammps/system.in.settings"

                    self.settings(settings_filename) 
                    
                    f.writelines([
                        "# Create bonds\n",
                        "bond_style harmonic\n",
                        f"bond_coeff 1 {self.var['general']['Cspring']} {lat_c} \n",
                        f"create_bonds many layer_1 layer_2 1 {lat_c-0.25} {lat_c+0.25}\n",
                        f"create_bonds many layer_3 layer_4 1 {lat_c-0.25} {lat_c+0.25}\n\n ",
                        "##########################################################\n",
                        "#-------------------Energy Minimization------------------#\n",
                        "##########################################################\n\n ",

                        "min_style       cg\n",
                        "minimize        1.0e-4 1.0e-8 100000 100000\n\n ",

                        "##########################################################\n",
                        "#------------------- Apply Constraints ------------------#\n",
                        "##########################################################\n\n ",


                        "#----------------- Apply Langevin thermostat -------------\n",
                        "group center union layer_2 layer_3\n ",
                        f"velocity        center create {self.var['general']['temproom']} 492847948\n",
                        f"fix             lang center langevin {self.var['general']['temproom']} {self.var['general']['temproom']} $(100.0*dt) 2847563 zero yes\n\n ",

                        "fix             nve_all all nve\n\n ",

                        "timestep        0.001\n",
                        "thermo          100\n\n ",

                        "compute COM_bot stage com\n",
                        "variable comx_bot equal c_COM_bot[1] \n",
                        "variable comy_bot equal c_COM_bot[2] \n",
                        "variable comz_bot equal c_COM_bot[3] \n\n ",

                        "fix             fstage2 layer_4 rigid single force * on on off torque * off off off\n",

                        "run 1000\n\n ",

                        f"variable omega equal {a}/1000\n",
                        "fix rot stage move rotate ${comx_bot} ${comy_bot} ${comz_bot} 0 0 1 ${omega}\n\n ",

                        "run             1000\n\n ",

                        "unfix fstage2\n",
                        "unfix rot\n\n",

                        "fix             fstage2 layer_4 rigid single force * off off on torque * off off off\n\n",

                        "fix             fsbot stage_bot setforce 0.0 0.0 0.0 \n ",
                        "velocity        stage_bot set 0.0 0.0 0.0 units box\n\n  ",
                        
                        f"variable Fatom equal -{force}/(count(layer_4)*1.602176565)\n",
                        "fix force tip_fix aveforce 0.0 0.0 ${Fatom}\n\n",

                        "run             10000\n\n ",
                        
                        "variable        fx   equal  f_fstage2[1]*1.602176565\n",
                        "variable        fy   equal  f_fstage2[2]*1.602176565\n",
                        "variable        fz   equal  f_force[3]*1.602176565\n\n",



                        "----------------- Output values -------------------------\n",
                        "fix             fc_ave all ave/time 1 1000 1000 v_fx v_fy v_fz file ../friction_measurements/fc_graphvgraph\n\n",

                        f"velocity        stage_top set {self.var['general']['scan_s']} 0.0 0.0 ",
                        "run             100000\n\n ",

                        "##########################################################\n",
                        "#-----------------------Write Data-----------------------#\n",
                        "##########################################################\n\n ",

                        "#----------------- Save final configuration in data file -\n",
                        "write_data     graphvgraph.data\n"
                    ])