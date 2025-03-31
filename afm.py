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

class afm:
    def __init__(self,input):
        
        # Generate data and variables

        var = read_config(input) # Read input file settings
        self.group = ['2D' ,'sub', 'tip']

        var['data'] = {mat: cifread(var[mat]['mat']) for mat in self.group}  # Read materials  
        var['pot'] = {mat: count_elemtypes(var['data'][mat]['pot_path']) for mat in self.group}  # Count potentials  

        var['data'] = {
            mat['natype']: sum(var['pot'][mat].values()) 
            for mat in self.group
        }

        # Build one layer of 2D material to generate important variables

        self.var = sheet(var)

        # Create file locations
        self.var['dir'] = f"scripts/{self.var['data']['2D'][1]}/size_{self.var['2D']['x']}x_{self.var['2D']['y']}y/sub_{self.var['sub']['amorph']}{self.var['data']['sub'][1]}/tip_{self.var['tip']['amorph']}{self.var['data']['tip'][1]}_r{self.var['tip']['r']}/K{self.var['general']['temproom']}"

        self.scripts = self.var['dir'] +"/scripts"

        dirs = ["visuals", "results", "system_build", "potentials", "scripts"]

        for d in dirs:
            Path(self.var['dir'], d).mkdir(parents=True, exist_ok=True)

        files = ["list_system", "list_load", "list_slide"]

        for f in files:
            with open(Path(self.var['dir'], "scripts", f), "w"):
                pass 
        

        #Expand to multiple layers if required

        self.var['ngroups'] = {}
        for l in self.var['2D']['layers']:
            self.directory[l] = Path(self.var['dir']) / f"l_{l}"

            for sub in ["data", "lammps"]:
                (self.directory[l] / sub).mkdir(parents=True, exist_ok=True)
            _ = stacking(self.var,l)
            self.var['ngroups'][l] = self.var['2D']['natype']*l + self.var['data']['sub']['nelements']*3 + self.var['data']['tip']['nelements']*3
            


        self.scan_angle = np.arange(self.var['tip']['scan_angle'][0],self.var['tip']['scan_angle'][1]+1,self.var['tip']['scan_angle'][2])

        # Limit the number of visual files generated
        self.dump_load = [self.var['general']['force'][i] for i in range(4, len(self.var['general']['force']), 5)]
        self.dump_slide = [self.scan_angle[i] for i in range(4, len(self.scan_angle), 5)]
        
        # Generate substrate and tip

        tip(self.var)
        sub(self.var)

    def system(self):
        for layer in self.var['2D']['layers']:

            [tip_x,tip_y,tip_z] = [self.var['dim']['xhi']/2, self.var['dim']['yhi']/2, 55+self.var['2D']['lat_c']*(layer-1)/2]# Tip placement
            tip_h       = round(self.var['tip']['r']*0.52)
            tipt        = tip_z+tip_h-3
            h_2D        = 10+self.var['2D']['lat_c']/3
            tip_c_bot   = tipt-2
            tip_c_max   = tipt   

            filename = self.directory[layer] +"/lammps/" + "system.lmp"
            with open(self.scripts + "/list_system", 'a') as f:
                f.write(f"{filename}\n")

            with open(filename, 'w') as f:
                init(f)
                f.writelines([
                    f"region box block {self.var['dim']['xlo']} {self.var['dim']['xhi']} {self.var['dim']['ylo']} {self.var['dim']['yhi']} -5 100\n",
                    f"create_box      {self.var['ngroups'][layer]} box\n\n",
                    "#----------------- Read data files -----------------------\n\n",
                    f"read_data       {self.var['dir']}/system_build/sub.lmp add append group sub\n",
                    f"read_data       {self.var['dir']}/system_build/tip.lmp add append shift {tip_x} {tip_y} {tip_z}  group tip\n",
                    f"read_data       {self.var['dir']}/system_build/{self.var['data']['2D']['mat']}_{layer}.lmp add append shift 0.0 0.0 {h_2D} group 2D\n\n"
                ])

                for t in range(max(self.var['data'][mat]['natype'] for mat in self.group)):
                    t+=1
                    f.write(f"group type_{t} type {t}\n") 
                
                # Identify atom regions
                f.writelines([
                    "\n#Identify the top atoms of AFM tip\n\n",
                    f"region          tip_fix block INF INF INF INF {tipt} INF units box\n",
                    "group           tip_fix region tip_fix\n\n",
                    "#Identify thermostat region of AFM tip\n\n",
                    f"region          tip_thermo block INF INF INF INF {tip_c_bot} {tip_c_max} units box\n",
                    "group           tip_thermo region tip_thermo\n\n",
                    "#Identify bottom atoms of Amorphous Silicon substrate\n\n",
                    "region          sub_fix block INF INF INF INF INF 2 units box\n",
                    "group           sub_fix region sub_fix\n\n",
                    "#Identify thermostat region of Amorphous Silicon substrate\n\n",
                    "region          sub_thermo block INF INF INF INF 2 5 units box\n",
                    "group           sub_thermo region sub_thermo\n\n",
                    "# Define sub groups and atom types\n\n"
                ])

                i = 1
                # Define sub groups and atom types
                for g in self.group:
                    for t in range(self.var['data'][g]['natype']):
                        if g == '2D':
                            t+=1
                            f.writelines([
                            f"group 2D_{t} intersect 2D type_{t}\n",
                            f"set group 2D_{t} type {i}\n",
                            ])
                            i+=1
                            for l in range(layer):
                                l+=1
                                zlo= h_2D + l*self.var['2D']['lat_c']/2 -1
                                zhi= zlo + self.var['2D']['lat_c']/2
                                f.writelines([
                                f"region layer_{l} block INF INF INF INF {zlo} {zhi} units box\n",
                                f"group layer_{l} region layer_{l} \n",
                                f"region layer_{l} delete\n",
                                f"group layer intersect 2D_{t} layer_{l}\n",
                                f"set group layer type {i}\n",
                                f"group layer_{l} delete\n"
                                f"group layer delete\n\n"
                                ])
                                i+=1
                                f.write(f"group 2D_{t} delete\n\n")
                        else:
                            t+=1
                            f.writelines([
                            f"group {g}_{t} intersect {g} type_{t}\n",
                            f"set group {g}_{t} type {i}\n",

                            f"group {g}_fix_{t} intersect {g}_fix type_{t}\n",
                            f"set group {g}_fix_{t} type {i+1}\n",
                            f"group {g}_fix_{t} delete\n\n",

                            f"group {g}_thermo_{t} intersect {g}_thermo type_{t}\n",
                            f"set group {g}_thermo_{t} type {i+2}\n",
                            f"group {g}_thermo_{t} delete\n\n"
                            ])
                            i+=3
                            f.write(f"group {g}_{t} delete\n\n")
            
                #generate potentials
                settings_filename = f"{self.directory[layer]}/lammps/system.in.settings"
                self.settings(settings_filename,layer) 

                f.writelines([
                "# Apply potentials\n\n",
                f"include        {self.directory[layer]}/lammps/system.in.settings\n\n",
                "#----------------- Create visualisation files ------------\n\n",
                f"dump            sys all atom 100 ./{self.var['dir']}/visuals/system_{layer}.lammpstrj\n\n",
                "#----------------- Minimize the system -------------------\n\n",
                "min_style       cg\n",
                "minimize        1.0e-4 1.0e-8 1000000 1000000\n\n",
                "timestep        0.001\n",
                "thermo          100\n\n",
                #----------------- Apply Nose-Hoover thermostat ----------
                "group           fixset union sub_fix tip_fix\n",
                "group           system subtract all fixset\n\n",
                f"velocity        system create {self.var['general']['temproom']} 492847948\n\n",
                "compute         temp_tip tip_thermo temp/partial 0 1 0\n",
                f"fix             lang_tip tip_thermo langevin {self.var['general']['temproom']} {self.var['general']['temproom']} $(100.0*dt) 699483 zero yes\n",
                "fix_modify      lang_tip temp temp_tip\n\n",
                "compute         temp_sub sub_thermo temp/partial 0 1 0\n",
                f"fix             lang_sub sub_thermo langevin {self.var['general']['temproom']} {self.var['general']['temproom']} $(100.0*dt) 2847563 zero yes\n",
                "fix_modify      lang_sub temp temp_sub\n\n",
                "fix             nve_all all nve\n\n",
                "fix             sub_fix sub_fix setforce 0.0 0.0 0.0 \n",
                "velocity        sub_fix set 0.0 0.0 0.0\n\n",
                "fix             tip_f tip_fix rigid/nve single force * off off off torque * off off off\n\n",
                "run             10000\n\n",
                "unfix           tip_f \n\n",
                "##########################################################\n",
                "#--------------------Tip Indentation---------------------#\n",
                "##########################################################\n",
                "#----------------- Displace tip closer -------------------\n\n",
                "displace_atoms  tip_all move 0.0 0.0 -20.0 units box\n\n",
                "#----------------- Apply constraints ---------------------\n\n",
                "#Fix the bottom layer of the base and the edges of the graphene\n\n",
                "fix             tip_f tip_fix rigid/nve single force * off off on torque * off off off\n\n",
                "#----------------- Set up initial parameters -------------\n\n",
                f"                variable find equal {self.var['general']['find']}\n",
                "variable        num_floads equal 500\n",
                "variable        r equal 0.0\n",
                "variable        f equal 0.0\n"
                "variable        fincr equal ${find}/(${num_floads})\n",
                "thermo_modify   lost ignore flush yes\n\n",
                "#----------------- Apply pressure to the tip -------------\n\n",
                "variable i loop ${num_floads}\n",
                "label loop_load\n\n",
                "variable f equal ${f}+${fincr} \n\n",
                "# Set force variable\n\n",
                "variable Fatom equal -v_f/(count(tip_fix)*1.602176565)\n",
                "fix forcetip tip_fix aveforce 0.0 0.0 ${Fatom}\n",
                "run 100 \n\n"
                "unfix forcetip\n\n",
                "next i\n",
                "jump SELF loop_load\n\n",
                "##########################################################\n",
                "#---------------------Equilibration----------------------#\n",
                "##########################################################\n\n",
                "fix forcetip tip_fix aveforce 0.0 0.0 ${Fatom}\n",
                "variable        dispz equal xcm(tip_fix,z)\n\n",
                "run 100 pre yes post no\n\n",
                "# Prepare to loop for displacement checks\n\n",
                "label check_r\n\n",
                "variable disp_l equal ${dispz}\n",
                "variable disp_h equal ${dispz}\n\n",
                "variable disploop loop 50\n",
                "label disp\n\n",
                "run 100 pre no post no\n\n",
                "if '${dispz}>${disp_h}' then 'variable disp_h equal ${dispz}'\n",
                "if '${dispz}<${disp_l}' then 'variable disp_l equal ${dispz}'\n\n",
                "next disploop\n",
                "jump SELF disp\n\n",
                "variable r equal ${disp_h}-${disp_l}\n\n",
                "# Check if r is less than 0.1\n\n",
                "if '${r} < 0.2' then 'jump SELF loop_end' else 'jump SELF check_r'\n\n",
                "# End of the loop\n\n",
                "label loop_end\n\n",
                f"write_data {self.directory[layer]}/data/load_{self.var['general']['find']}N.data"
                ])
     
    def load(self):
        for layer in self.var['2D']['layers']:

            for force in self.var['general']['force']:
                dump = False
                if force in self.dump_load:
                    dump = True

                filename = f"{self.directory[layer]}/lammps/load_{force}N.lmp"
                
                with open(self.scripts + "/list_load", 'a') as f:
                    f.write(f"{filename}\n")

                with open(filename, 'w') as f:
                    init(f)
                    f.writelines([
                    # f"include {self.directory[layer]}/lammps/in.init\n\n",
                    f"read_data       {self.directory[layer]}/data/load_{self.var['general']['find']}N.data # Read system data\n\n",
                    f"include         {self.directory[layer]}/lammps/system.in.settings\n\n",
                    "#----------------- Create visualisation files  ------------\n\n",
                    ])

                    if dump == True:
                        f.write(f"dump            sys all atom 100 ./{self.var['dir']}/visuals/load_{force}N_l{layer}.lammpstrj\n\n",)


                    f.writelines([
                    "##########################################################\n",
                    "#--------------------Tip Indentation---------------------#\n",
                    "##########################################################\n",
                    "#----------------- Apply constraints   ---------------------\n\n",
                    "#Fix the bottom layer of the base and the edges of the graphene\n\n",
                    "fix             sub_fix sub_fix setforce 0.0 0.0 0.0 \n",
                    "fix             tip_f tip_fix rigid/nve single force * off off on torque * off off off\n\n",
                    "#----------------- Apply Langevin thermostat   -------------\n\n",
                    "compute         temp_tip tip_thermo temp/partial 0 1 0\n",
                    f"fix             lang_tip tip_thermo langevin {self.var['general']['temproom']}   {self.var['general']['temproom']} $(100.0*dt) 699483 zero yes\n",
                    "fix_modify      lang_tip temp temp_tip\n\n",
                    "compute         temp_base sub_thermo temp/partial 0 1 0\n",
                    f"fix             lang_bot sub_thermo langevin {self.var['general']['temproom']}   {self.var['general']['temproom']} $(100.0*dt) 2847563 zero yes\n",
                    "fix_modify      lang_bot temp temp_base\n\n",
                    "fix             nve_all all nve\n",
                    "timestep        0.001\n",
                    "thermo          100\n\n",

                    "#----------------- Set up initial parameters   -------------\n\n",
                    f"variable find equal {force}\n",
                    "variable num_floads equal 500\n",
                    "variable r equal 0.0\n",
                    f"variable f equal {self.var['general']['find']}\n",
                    "variable fincr equal (${find}-${f})/${num_floads}\n",
                    "thermo_style    custom step temp v_f pe ke etotal press\n",
                    # "thermo_modify   lost ignore flush yes\n\n",
                    "run 0\n\n",
                    "#----------------- Apply pressure to the tip   -------------\n\n",
                    "variable i loop ${num_floads}\n",
                    "label loop_load\n\n",
                    "# Set force variable\n\n",
                    "variable Fatom equal -v_f/(count(tip_fix)*1.602176565)\n",
                    "fix forcetip tip_fix aveforce 0.0 0.0 ${Fatom}\n",
                    "run 100 \n\n",
                    "unfix forcetip\n\n",
                    "variable f equal ${f}+${fincr} \n\n",
                    "next i\n",
                    "jump SELF loop_load\n\n",
                    "##########################################################\n",
                    "#---------------------Equilibration----------------------#\n",
                    "##########################################################\n\n"    ,
                    "fix forcetip tip_fix aveforce 0.0 0.0 ${Fatom}\n",
                    "variable        dispz equal xcm(tip_fix,z)\n\n",
                    "run 100 pre yes post no\n\n",
                    "# Prepare to loop for displacement checks\n\n",
                    "label check_r\n\n",
                    "variable disp_l equal ${dispz}\n",
                    "variable disp_h equal ${dispz}\n\n",
                    "variable disploop loop 50\n",
                    "label disp\n\n",
                    "run 100 pre no post no\n\n",
                    "if '${dispz}>${disp_h}' then 'variable disp_h equal ${dispz}   '\n",
                    "if '${dispz}<${disp_l}' then 'variable disp_l equal ${dispz}   '\n\n",
                    "next disploop\n",
                    "jump SELF disp\n\n",
                    "variable r equal ${disp_h}-${disp_l}\n\n",
                    "# Check if r is less than 0.1\n\n",
                    "if '${r} < 0.2' then 'jump SELF loop_end' else 'jump SELF          check_r'\n\n",
                    "# End of the loop\n\n",
                    "label loop_end\n\n",
                    f"write_data {self.directory[layer]}/data/load_{force}N.data"
                    ])

    def slide(self):
        # 1.602176565x10^(-9) N = 1 eV/Angstrom and 1 Angstrom = 10^(-10) m
        springeV  = self.var['tip']['cspring']/16.02176565 # Spring constant in eV/A^2
        # 1.602176565 nN = 1 eV/Angstrom and 1 Angstrom = 10^(-10) m and 1 ps = 10^ (-12) s
        DspringeV  = self.var['tip']['dspring']/0.01602176565 # Spring damper in ev/(A^2/ps)
        tipps = self.var['tip']['s']/100 # in A/ps

        for layer in self.var['2D']['layers']:
            for force in self.var['general']['force']:
                
                if force == self.var['tip']['scan_angle'][3]:
                    dump = False
                    scan_angle = self.scan_angle
                else:
                    dump = False
                    scan_angle = [0]
                    if force in self.dump_load:
                        dump = True
                    
                for a in scan_angle:
                    if a in self.dump_slide:
                        dump = True
                        
                    spring_x = np.cos(np.deg2rad(a))
                    spring_y = np.sin((a))
                    filename = self.directory[layer] + "/lammps/" + "slide_" + str(force) +"N_"+ str(self.var['tip']['s']) + "ms_" + str(a) + "deg.lmp"
                    with open(self.scripts + "/list_slide", 'a') as f:
                        f.write(f"{filename}\n")
                    with open(filename, 'w') as f:
                        init(f)
                        f.writelines([
                        # f"include {self.directory[layer]}/lammps/in.init\n\n",
                        f"read_data       {self.directory[layer]}/data/load_{force}N.data # Read system data\n\n",
                        f"include         {self.directory[layer]}/lammps/system.in.settings\n\n",
                        "#----------------- Create visualisation files ------------\n\n",
                        ])
                        #CHANGE THIS
                        if dump == True:
                            f.write(f"dump            sys all atom 100 ./{self.var['dir']}/visuals/slide_{force}nN_{a}angle_{self.var['tip']['s']}ms_l{layer}.lammpstrj\n\n")
                        
                        f.writelines([
                        
                        "##########################################################\n",
                        "#--------------------Tip Indentation---------------------#\n",
                        "##########################################################\n",
                        "#----------------- Apply constraints ---------------------\n\n",
                        "#Fix the bottom layer of the base and the edges of the graphene\n\n",
                        "fix             sub_fix sub_fix setforce 0.0 0.0 0.0 \n",
                        "fix             tip_f tip_fix rigid/nve single force * off off on torque * off off off\n\n",
                        "#----------------- Apply Langevin thermostat -------------\n\n",
                        "compute         temp_tip tip_thermo temp/partial 0 1 0\n",
                        f"fix             lang_tip tip_thermo langevin {self.var['general']['temproom']} {self.var['general']['temproom']} $(100.0*dt) 699483 zero yes\n",
                        "fix_modify      lang_tip temp temp_tip\n\n",
                        "compute         temp_base sub_thermo temp/partial 0 1 0\n",
                        f"fix             lang_bot sub_thermo langevin {self.var['general']['temproom']} {self.var['general']['temproom']} $(100.0*dt) 2847563 zero yes\n",
                        "fix_modify      lang_bot temp temp_base\n\n",
                        "fix             nve_all all nve\n",
                        "timestep        0.001\n",
                        "thermo          100\n\n",
                        "#----------------- Apply pressure to the tip -------------\n\n",
                        f"variable        Ftotal          equal -{force}/1.602176565\n",
                        "variable        Fatom           equal v_Ftotal/count(tip_fix)\n",
                        "fix             forcetip tip_fix aveforce 0.0 0.0 ${Fatom}\n\n",
                        "##########################################################\n",
                        "#------------------------Compute-------------------------#\n",
                        "##########################################################\n\n",
                        "#----------------- Calculate total friction --------------\n\n",
                        "variable        fz_tip   equal  f_forcetip[3]*1.602176565\n\n",
                        "variable        fx_spr   equal  f_spr[1]*1.602176565\n\n",
                        "variable        fy_spr   equal f_spr[2]*1.602176565\n\n",
                        f"fix             fc_ave all ave/time 1 1000 1000 v_fz_tip v_fx_spr v_fy_spr file ./{self.var['dir']}/results/fc_ave_slide_{force}nN_{a}angle_{self.var['tip']['s']}ms_l{layer}\n\n",
                        "##########################################################\n",
                        "#---------------------Spring Loading---------------------#\n",
                        "##########################################################\n\n",
                        "#----------------- Add damping force ---------------------\n\n",
                        f"fix             damp tip_fix viscous {DspringeV}\n\n",
                        "#------------------Add lateral harmonic spring------------\n\n",
                        f"fix             spr tip_fix smd cvel {springeV} {tipps} tether {spring_x} {spring_y} NULL 0.0\n\n",
                        "run 200000\n\n",
                        "unfix spr\n\n",
                        "variable        fx_spr   equal  0\n",
                        "run 100000\n\n",
                        f"fix             spr tip_fix smd cvel  {springeV} {tipps} tether -{spring_x} -{spring_y} NULL 0.0\n\n",
                        "run 200000\n",
                        f"write_data {self.directory[layer]}/data/slide_{force}nN_{a}angle_{self.var['tip']['s']}ms.data"
                        ])
    
    def settings(self,filename,layer):
        """Writes the LAMMPS input file content to the specified filename.
        Args:
        filename (str): The name of the file to write to.
        """
        with open(filename, 'w') as f:
            group = ['sub','tip','2D']
            if self.var['flake']['flake'] == True:
                group.append("flake")

            
            sub_elem = self.var['data']['sub']['elements'].copy()
            tip_elem = self.var['data']['tip']['elements'].copy()
            twoD_elem = [str(sublist[0]) for sublist in self.elem2D]
            elems = [sub_elem,tip_elem,twoD_elem]

            # number elements

            for arr in elems:
                count = {}
                result=[]
                for i in range(len(arr)):
                    element = arr[i]
                    count[element] = count.get(element, 0) + 1
                    result.append(element + str(count[element]))
                for i in range(len(arr)):
                    element = arr[i]
                    if count[element] > 1:
                        arr[i]=result[i]
            
            # Set masses and create groups
            i = 0
            group_def = {}
            elemgroup = {}
            potentials= {}

            for g in group:

                if g == 'sub' or g == 'tip':
                    for t in range(self.var['data'][g]['nelements']):
                        m = self.var['data'][g]['elements'][t]
                        group_def.update({
                            i:   [f"{g}_b_t{t+1}", str(i+1), str(m), sub_elem[t]],
                            i+1: [f"{g}_fix_t{t+1}", str(i+2), str(m), sub_elem[t]],
                            i+2: [f"{g}_thermo_t{t+1}", str(i+3), str(m), sub_elem[t]]
                        })
                        elemgroup[g][m].extend([i, i+1, i+2])           
                        i+=3

                elif g == '2D':
                    for t in range(self.var['2D']['natype']):
                        m = self.elem2D[t]

                        for l in range(layer):
                            group_def.update({i: [str(g)+"_l" + str(l+1) + "_t"+str(t+1), str(i+1),str(m),twoD_elem[t],l+1]})
                            i+=1
                            elemgroup[g][l][m].append(i)

                for m in self.var['data'][g]['elements']: 
                    mass=data.atomic_masses[data.atomic_numbers[m]] 
                    if g =='2D':
                        f.write(f"mass {elemgroup[g][0][m][0]}*{elemgroup[g][-1][m][-1]} {mass} #{m} layer {l+1}\n")
                    else:
                        f.write(f"mass {elemgroup[g][m][0]}*{elemgroup[g][m][-1]} {mass} #{m} layer {l+1}\n")

                all = [group_def[i][1] for i in range(self.var['ngroups'][layer]) if g in group_def[i][0]]
                f.write(f"group {g}_all type {' '.join(all)}\n")
                
                if g == 'sub' or g == 'tip':
                    for n in ["_fix", "_thermo"]:
                        sub_group = [group_def[i][1] for i in range(self.var['ngroups'][layer]) if g+n in group_def[i][0]]
                        f.write(f"group {g}{n} type {' '.join(sub_group)}\n")

                if g == '2D':
                    for l in range(layer):
                        potentials[g][l] = [
                            group_def[i][3] if any("2D_l"+str(l+1) in group_def[i][0]) else "NULL"
                            for i in range(self.var['ngroups'][layer])
                        ]
                else:
                    potentials[g]=[group_def[i][2] if any(g in group_def[i][0]) else "NULL"
                                   for i in range(self.var['ngroups'][layer])]

            f.writelines(["group mobile union tip_thermo sub_thermo\n",
                          f"pair_style hybrid {'sw ' * (l + 2)}lj/cut 8.0\n"
                        ])
            
            i = 1 
            for key in group:
                if key == '2D':
                    for l in range(layer):  
                        f.write(f"pair_coeff * * sw {i} {self.var['dir']}/potentials/{self.var['data'][key]['formula']}.sw {potentials[key][l]} # interlayer {key.capitalize()} Layer {l+1}\n")
                        i += 1 
                else:
                    f.write(f"pair_coeff * * sw {i} {self.var['dir']}/potentials/{self.var['data'][key]['formula']}.sw {potentials[key]} # interlayer {key.capitalize()}\n")
                    i += 1 


            # This is for potentials that include VdW interactions
            # layer_potentials = [
            #     group_def[i][2] if any(group_def[i][2] == element and "2D" in group_def[i][0] for element in set(self.data["2D"][3])) else "NULL"
            #     for i in range(1, self.ngroups + 1)
            # ]
            # f.write(f"pair_coeff * * sw 3 Potentials/{self.data['2D'][1]}.sw {'  '.join(layer_potentials)} #interlayer substrate\n")
            
            ############ COPY POTENTIALS TO FILE
            
            # file = [Path(__file__).parent / f"potentials/{self.var['data']['2D'][1]}.sw",Path(__file__).parent / f"potentials/{self.var['data']['tip'][1]}.sw",Path(__file__).parent / f"potentials/{self.var['data']['sub'][1]}.sw"]
            # shutil.copy2(file, self.var['dir'] / 'potentials')

            #Consider LabelMaps to reduce amount of lines but doesn't really matter
            
            for t in self.var['data']['2D']['elements']:    
                for key in ('sub','tip'):
                    for s in self.var['data'][key]['elements']:
                        e,sigma = LJparams(t,s)
                        if len(elemgroup['2D'][t]) == 1 and layer == 1:
                            f.write(f"pair_coeff {elemgroup[key][t][0]}*{elemgroup[key][t][-1]} {elemgroup['2D'][0][t][0]} lj/cut {e} {sigma}\n")
                        else:  
                            f.write(f"pair_coeff {elemgroup[key][t][0]}*{elemgroup[key][t][-1]} {elemgroup['2D'][0][t][0]}*{elemgroup['2D'][-1][t][-1]} lj/cut {e} {sigma}\n")

                if layer>1:
                    for s in self.var['data']['2D']['elements']:
                        e,sigma = LJparams(s,t)

                        for l in range(1,layer):
                            t1 = f"{elemgroup['2D'][l][t][0]}*{elemgroup['2D'][l][t][-1]}"
                            t2 = f"{elemgroup['2D'][l+1][s][0]}*{elemgroup['2D'][-1][s][-1]}"

                            if elemgroup['2D'][l][t][0] == elemgroup['2D'][l][t][-1]:
                                t1 = f"{elemgroup['2D'][l][t][0]}"
                            if elemgroup['2D'][l+1][s][0] == elemgroup['2D'][-1][s][-1]:
                                t2 = f"{elemgroup['2D'][l+1][s][0]}"

                            if elemgroup['2D'][l][t][0]>elemgroup['2D'][l+1][s][0]:
                                t1, t2 = t2, t1

                            f.write(f"pair_coeff {t1} {t2} lj/cut {e} {sigma} \n")

            for s in self.var['data']['sub']['elements']:
                for t in self.var['data']['tip']['elements']:
                    e,sigma = LJparams(s,t)
                    f.write(f"pair_coeff {elemgroup['sub'][t][0]}*{elemgroup['sub'][t][-1]} {elemgroup['tip'][t][0]}*{elemgroup['tip'][t][-1]}  lj/cut {e} {sigma} \n")
                    
    def pbs(self):

        pbs_type = ['system','load','slide']
        for type in pbs_type:
            filename = self.scripts + self.data['2D'][1] + f"_{type}.pbs"
            PBS = '"${PBS_ARRAY_INDEX}p"'
            PBS_log = "{PBS_ARRAY_INDEX}"
            with open(self.scripts + f"list_{type}", 'r' ) as f:
                n = len(f.readlines())
            with open(filename,'w') as f: 
                f.writelines([
                    "#!/bin/bash\n",
                    "#PBS -l select=1:ncpus=32:mem=62gb:mpiprocs=32:cpu_type=rome\n",
                    "#PBS -l walltime=08:00:00\n",
                    f"#PBS -J 1-{n}\n",
                    f"#PBS -o /rds/general/user/mv923/home/{self.data['2D'][1]}/\n",
                    f"#PBS -e /rds/general/user/mv923/home/{self.data['2D'][1]}/\n\n",

                    "module purge\n",
                    "module load tools/dev\n",
                    "module load LAMMPS/23Jun2022-foss-2021b-kokkos\n",
                    "#module load OpenMPI/4.1.4-GCC-11.3.0\n\n",

                    "#Go to the temp directory (ephemeral) and create a new folder for this run\n",
                    "cd $EPHEMERAL\n\n",


                    "# $PBS_O_WORKDIR is the directory where the pbs script was sent from. Copy everything from the work directory to the temporary directory to prepare for the run\n\n",

                    f"mpiexec lmp -l {self.data['2D'][1]}/${PBS_log}.log -in $(sed -n {PBS} {self.var['dir']}/scripts/list_{type})\n\n",

                    # "#After the end of the run copy everything back to the parent directory\n",
                    # f"cp -r ./{self.var['dir']}/ $PBS_O_WORKDIR/{self.var['dir']}\n\n"
                ])

        filename = self.scripts + self.data['2D'][1] + "_transfer.pbs"
        with open(filename,'w') as f: 
            f.writelines([
                "#!/bin/bash\n",
                "#PBS -l select=1:ncpus=1:mem=62gb:cpu_type=rome\n",
                "#PBS -l walltime=00:30:00\n\n",
                f"#PBS -o /rds/general/user/mv923/home/{self.data['2D'][1]}/\n",
                f"#PBS -e /rds/general/user/mv923/home/{self.data['2D'][1]}/\n\n",

                "cd $HOME\n",
                f"mkdir -p logs_{self.data['2D'][1]}/\n\n",
                "cd $EPHEMERAL\n",
                f"mkdir -p {self.var['dir']}/\n\n",

                f"cp -r $PBS_O_WORKDIR/{self.var['dir']}/* {self.var['dir']}\n",
                "cp -r $PBS_O_WORKDIR/Potentials/ .\n"
            ])

        filename = self.scripts + self.data['2D'][1] + "_transfer2.pbs"
        with open(filename,'w') as f: 
            f.writelines([
                "#!/bin/bash\n",
                "#PBS -l select=1:ncpus=1:mem=62gb:cpu_type=rome\n",
                "#PBS -l walltime=00:30:00\n\n",
                f"#PBS -o /rds/general/user/mv923/home/logs_{self.data['2D'][1]}/\n",
                f"#PBS -e /rds/general/user/mv923/home/logs_{self.data['2D'][1]}/\n\n",

                "cd $EPHEMERAL\n",
                "#After the end of the run copy everything back to the parent directory\n",
                f"cp -r ./{self.var['dir']}/* $PBS_O_WORKDIR/{self.var['dir']}\n\n"
            ])

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description="Run the program with an input file.")
#     parser.add_argument("input_file", type=str, help="The input file to process")
#     args = parser. 