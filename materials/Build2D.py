# structure for 1H-MoS2, with 12 atom types
import numpy as np
from mp_api.client import MPRester
from ase import io, data
# orientation, a=armchair, z=zigzag


def build2D(var,layers):
    
    
    nlayer = layers
    x = var['2D']['x']
    y = var['2D']['y']
    natype = var['data']['2D'][4]*2


    latconst = var['data']['2D'][0].lattice.a
    thickness = var['data']['2D'][0].lattice.c/2 #distance between Mo atoms in seperate sheets
    # bond = var['data']['2D'][0].distance_matrix[1,4]
    bond = latconst / (2*np.sqrt(3/7))
    orientation = 'z'
    # # structure parameters
    b = latconst/np.sqrt(3)
    h = np.sqrt(bond*bond-b*b)
    ci = 2.0*h

    elem = var['data']['2D'][3]

    amass = {}
    aelem = {}
    atype = {}
    t=0
    
    for i in range(natype):
        if (i+2) % 3 == 0:
            amass[t] = data.atomic_masses[data.atomic_numbers[elem[0]]]
            aelem[t] = elem[0]
            atype[i] = t+1
            t+=1

    for i in range(natype):  
        if (i+2) % 3 != 0:       
            amass[t] = data.atomic_masses[data.atomic_numbers[elem[1]]]
            aelem[t] = elem[1]
            atype[i] = t+1
            t+=1


    lenxyz = [None]*3
    lenxyz[0] = x
    lenxyz[1] = y
    lenxyz[2] = nlayer*thickness

    # armchair and zigzag are different by the following coordinate transformation.
            #    0  1  0
            #    1  0  0
            #    0  0 -1
    # i.e., x --> y, y--> x, and z --> -z

    coord = [None]*nlayer
    for i in range(1,nlayer+1):
        #even layer 
        i-=1
        z0 = i * thickness

        if i % 2 == 0:
            X1 = np.array([
                [0.0*np.sqrt(3.0)*b],
                [0.5*np.sqrt(3.0)*b],
                ])

            Y1 = np.array([
                [0.0*b],
                [0.5*b],
                [0.0*b],
                [1.5*b],
                [2.0*b],
                [1.5*b],
                ])
            
        if i % 2 != 0:

            X1 = np.array([
                [0.5*np.sqrt(3.0)*b],
                [0.0*np.sqrt(3.0)*b],
                ])

            Y1 = np.array([
                [0.5*b],
                [0.0*b],
                [0.5*b],
                [2.0*b],
                [1.5*b],
                [2.0*b],
                ])

        X = np.vstack((np.tile(X1, (3,1)), np.tile(X1+latconst, (3,1)))) 
        Y = np.vstack([Y1,Y1])

        Z1 = ([
            [0.0*ci- z0],
            [-0.5*ci - z0],
            [-1*ci - z0],
        ])

        Z = np.tile(Z1,(4,1))  

        if orientation == 'z':
            coord[i] = np.hstack([X,Y,Z])
            box_y = 3.0*b
            box_x = np.sqrt(3.0)*b*2.0 
        if orientation == 'a':
            coord[i] = np.hstack([Y,X,-Z])
            box_x = 3.0*b
            box_y = np.sqrt(3.0)*b * 2.0



    # reset size
    nx = int(lenxyz[0]//box_x)
    lenxyz[0] = nx*box_x
    ny = int(lenxyz[1]//box_y)
    lenxyz[1] = ny*box_y


    # total atom number
    ntot = nlayer * natype * nx * ny

    # allocate ntot-related arrays
    xalat = [None]*ntot
    yalat = [None]*ntot
    zalat = [None]*ntot
    ntype = [None]*ntot

    k=0
    # # generate structure
    for i in range(nx):
        for j in range(ny):
            x0 = i * box_x
            y0 = j * box_y
            z0 = 0.0
            # k = i*ny*natype*nlayer + j*natype*nlayer
            for l in range(nlayer):
                for ii in range(natype):

                    xalat[k] = x0 + coord[l][ii,0]
                    yalat[k] = y0 + coord[l][ii,1]
                    zalat[k] = z0 + coord[l][ii,2]
                    ntype[k] =  atype[ii]
                    k = k + 1

    # # shift to center of simulation box
    # for i in range(ntot):
    #     zalat[i] = zalat[i] + 0.5 * lenxyz[2]

    # small shift, avoid possible boundary effects
    xmin = 1000.0
    ymin = 1000.0
    for i in range(ntot):
        if xalat[i] < xmin:
            xmin = xalat[i]

        if yalat[i] < ymin:
            ymin = yalat[i]                                                          

    for i in range(ntot):
        xalat[i] = xalat[i] - xmin + 0.1
        yalat[i] = yalat[i] - ymin + 0.1

    filename = f"{var['dir']}/system_build/{var['data']['2D'][1]}_{nlayer}.lmp"
    with open(filename, 'w') as f:
        f.writelines([
            f"LAMMPS structure data file for MoS2 with {nlayer} layers\n\n",
            f"{ntot} atoms\n",
            f"{natype} atom types\n\n",
            f"0.0 {lenxyz[0]} xlo xhi\n",
            f"0.0 {lenxyz[1]} ylo yhi\n",
            f"{-20} {lenxyz[2]+10} zlo zhi\n\n",

            "Masses \n\n",
        ])

        for i in range(natype):
            f.write(f"{i+1} {amass[i]} #{aelem[i]}\n")

        f.write("\nAtoms\n\n")
        for i in range(ntot):
            f.write(f"{i+1} {ntype[i]} {xalat[i]} {yalat[i]} {zalat[i]}\n")                               

    return aelem, natype
