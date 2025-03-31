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
                    m = var['2D']['elem2D'][n][0]
                    g+=1
                    group_def.update({g: [f"2D_l{l+1}_t{n}", str(g),str(m),arr[n-1],l+1]})
                    elemgroup.setdefault(l, {}).setdefault(m, []).append(g)
                    c = count

        print(group_def)
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