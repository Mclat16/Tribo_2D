from ase import data
from .tools import *


def settings_sheet(var,filename,layer):
    """Writes the LAMMPS input file content to the specified filename.
    Args:
    filename (str): The name of the file to write to.
    """    
    with open(filename, 'w') as f:
        
        group_def = {}
        elemgroup = {}
        potentials= {}
        
        arr = {}
        count = {}

        i=0    
        for element, count in var['pot']['2D'].items():
            if not count or count == 1:
                arr[i] = element
                i+=1
            else:
                for t in range(1,count+1):
                    arr[i] = element + str(t)
                    i+=1

        i = 0 
        c = 0
        g = 0
        for element,count in var['pot']['2D'].items():
            i += c
            for l in range(layer):
                for t in range(1,count+1):
                    n = i + t
                    m = element
                    g+=1
                    group_def.update({g: [f"2D_l{l+1}_t{n}", str(g),str(m),arr[n-1],l+1]})
                    elemgroup.setdefault(l, {}).setdefault(m, []).append(g)
                    c = count

        for m in var['data']['2D']['elem_comp']: 
                mass=data.atomic_masses[data.atomic_numbers[m]] 
                f.write(f"mass {elemgroup[0][m][0]}*{elemgroup[layer-1][m][-1]} {mass} #{m}\n")

        f.write(f"pair_style hybrid {var['data']['2D']['pot_type'] * layer} lj/cut 8.0\n")

        for l in range(layer):
            potentials[l] = [
                group_def[i][3] if group_def[i][4]==l+1 else "NULL"
                for i in range(1,var['2D']['natype']*layer+1)
            ]            
            f.write(f"pair_coeff * * {var['data']['2D']['pot_type']} {l+1} {var['data']['2D']['pot_path']} {'  '.join(potentials[l])} # interlayer '2D' Layer {l+1}\n")
        
        if var['data']['2D']['pot_type'] == 'sw':
            for t in var['data']['2D']['elem_comp']:
                for s in var['data']['2D']['elem_comp']:
                    e,sigma = LJparams(s,t)
                    lat_c = max(lat_c,sigma)
                    for l in range(layer-1):
                        t1 = f"{elemgroup[l][t][0]}*{elemgroup[l][t][-1]}"
                        t2 = f"{elemgroup[l+1][s][0]}*{elemgroup[layer-1][s][-1]}"
                        if elemgroup[l][t][0] == elemgroup[l][t][-1]:
                            t1 = f"{elemgroup[l][t][0]}"
                        if elemgroup[l+1][s][0] == elemgroup[layer-1][s][-1]:
                            t2 = f"{elemgroup[l+1][s][0]}"
                        if elemgroup[l][t][0]>elemgroup[l+1][s][0]:
                            t1, t2 = t2, t1
                        f.write(f"pair_coeff {t1} {t2} lj/cut {e} {sigma} \n")

def settings_afm(var,layer):

    """Writes the LAMMPS input file content to the specified filename.
    Args:
    filename (str): The name of the file to write to.
    """
    filename = f"{var['dir']}/l_{l}/lammps/system.in.settings"
    with open(filename, 'w') as f:
        group = ['sub','tip','2D']

        arr = {}
        count = {}
        elemgroup = {}
        group_def = {}
         
        for gr in group:
            i=0
            for element, count in var['pot'][gr].items():
                if not count or count == 1:
                    arr[gr][i] = element
                    i+=1
                else:
                    for t in range(1,count+1):
                        arr[gr][i] = element + str(t)
                        i+=1

            i = 0 
            c = 0
            g = 1
            for element,count in var['pot'][gr].items():
                i += c
                if gr == '2D':
                    for l in range(layer):
                        for t in range(1,count+1):
                            n = i + t
                            m = element
                            mass=data.atomic_masses[data.atomic_numbers[m]] 
                            group_def.update({g: [f"{gr}_l{l+1}_t{n}", str(g),str(m),arr[g][n-1],l+1]})
                            elemgroup.setdefault(gr, {}).setdefault(l, {}).setdefault(m, []).append(g)
                            c = count
                            g+=1

                else:
                    for t in range(1,count+1):
                        n = i + t
                        m = element
                        group_def.update({
                        g:   [f"{gr}_b_t{t+1}", str(i+1), str(m),      arr[g][n-1]],
                        g+1: [f"{gr}_fix_t{t+1}", str(i+2), str(m),    arr[g][n-1]],
                        g+2: [f"{gr}_thermo_t{t+1}", str(i+3), str(m), arr[g][n-1]]
                        })         
                        elemgroup.setdefault(gr, {}).setdefault(m, []).append([g, g+1, g+2]) 
                        g+=3
                        c = count
            
            for m in var['data'][g]['elem_comp']: 
                mass=data.atomic_masses[data.atomic_numbers[m]] 
                if gr =='2D':
                    f.write(f"mass {elemgroup[gr][0][m][0]}*{elemgroup[gr][-1][m][-1]} {mass} #{m} {gr}\n")
                else:
                    f.write(f"mass {elemgroup[gr][m][0]}*{elemgroup[gr][m][-1]} {mass} #{m} {gr}\n")

                all = [group_def[i][1] for i in range(var['ngroups'][layer]) if gr in group_def[i][0]]
                f.write(f"group {g}_all type {' '.join(all)}\n")
                
                if g == 'sub' or g == 'tip':
                    for n in ["_fix", "_thermo"]:
                        sub_group = [group_def[i][1] for i in range(var['ngroups'][layer]) if gr+n in group_def[i][0]]
                        f.write(f"group {g}{n} type {' '.join(sub_group)}\n")

            f.writelines("group mobile union tip_thermo sub_thermo\n",
                    f"pair_style hybrid {var['data']['sub']['pot_type']} {var['data']['tip']['pot_type']} {var['data']['2D']['pot_type'] * layer} lj/cut 8.0\n")
        
        potentials = {}
        for g in group:
            t+=1
            if g == '2D':
                for l in range(layer):
                    potentials[l] = [
                        group_def[i][3] if group_def[i][4]==l+1 else "NULL"
                        for i in range(1,var['ngroups'][layer])
                    ]     

                    f.write(f"pair_coeff * * {var['data']['2D']['pot_type']} {t+l+1} {var['data']['2D']['pot_path']} {'  '.join(potentials[l])} # interlayer '2D' Layer {l+1}\n")
            else:
                potentials[g]=[group_def[i][2] if any(g in group_def[i][0]) else "NULL"
                for i in range(1,var['ngroups'][layer])]
                f.write(f"pair_coeff * * sw {t} {var['data']['2D']['pot_type']}.sw {potentials[g]} # interlayer {g.capitalize()}\n")

        for t in var['data']['2D']['elements']:    
            for key in ('sub','tip'):
                for s in var['data'][key]['elements']:
                    e,sigma = LJparams(t,s)
                    if len(elemgroup['2D'][t]) == 1 and layer == 1:
                        f.write(f"pair_coeff {elemgroup[key][t][0]}*{elemgroup[key][t][-1]} {elemgroup['2D'][0][t][0]} lj/cut {e} {sigma}\n")
                    else:  
                        f.write(f"pair_coeff {elemgroup[key][t][0]}*{elemgroup[key][t][-1]} {elemgroup['2D'][0][t][0]}*{elemgroup['2D'][-1][t][-1]} lj/cut {e} {sigma}\n")
            if layer>1:
                for s in var['data']['2D']['elements']:
                    e,sigma = LJparams(s,t)
                    for l in range(layer-1):
                        t1 = f"{elemgroup['2D'][l][t][0]}*{elemgroup['2D'][l][t][-1]}"
                        t2 = f"{elemgroup['2D'][l+1][s][0]}*{elemgroup['2D'][-1][s][-1]}"
                        if elemgroup['2D'][l][t][0] == elemgroup['2D'][l][t][-1]:
                            t1 = f"{elemgroup['2D'][l][t][0]}"
                        if elemgroup['2D'][l+1][s][0] == elemgroup['2D'][-1][s][-1]:
                            t2 = f"{elemgroup['2D'][l+1][s][0]}"
                        if elemgroup['2D'][l][t][0]>elemgroup['2D'][l+1][s][0]:
                            t1, t2 = t2, t1
                        f.write(f"pair_coeff {t1} {t2} lj/cut {e} {sigma} \n")
        for s in var['data']['sub']['elements']:
            for t in var['data']['tip']['elements']:
                e,sigma = LJparams(s,t)
                f.write(f"pair_coeff {elemgroup['sub'][t][0]}*{elemgroup['sub'][t][-1]} {elemgroup['tip'][t][0]}*{elemgroup['tip'][t][-1]}  lj/cut {e} {sigma} \n")