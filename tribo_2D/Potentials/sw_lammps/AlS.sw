# by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

# SW parameters for AlS       , used by LAMMPS
# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3,
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q,tol

# sw2 and sw3, for Al-S -S  and S -Al-Al
Al1 S1  S1     1.000   1.618   2.032  46.910   1.000  -0.184  11.476   2.112  4  0   0.000
S1  Al1 Al1    1.000   1.618   2.032  46.910   1.000  -0.184  11.476   2.112  4  0   0.000
Al2 S2  S2     1.000   1.618   2.032  46.910   1.000  -0.184  11.476   2.112  4  0   0.000
S2  Al2 Al2    1.000   1.618   2.032  46.910   1.000  -0.184  11.476   2.112  4  0   0.000

# sw2, for Al-Al
Al1 Al2 Al2    1.000   1.280   2.735   0.000   1.000   0.000   4.575   8.388  4  0   0.000
Al2 Al1 Al1    1.000   1.280   2.735   0.000   1.000   0.000   4.575   8.388  4  0   0.000

# sw3 for Al-Al-S 
Al1 Al2 S1     1.000   0.000   0.000  26.090   1.000  -0.459   0.000   0.000  4  0   0.000
Al1 S1  Al2    1.000   0.000   0.000  26.090   1.000  -0.459   0.000   0.000  4  0   0.000
Al2 Al1 S2     1.000   0.000   0.000  26.090   1.000  -0.459   0.000   0.000  4  0   0.000
Al2 S2  Al1    1.000   0.000   0.000  26.090   1.000  -0.459   0.000   0.000  4  0   0.000

# zero terms
S1  S1  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S1  Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S1  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S1  Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  Al1 S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  Al1 S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  Al1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S2  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S2  Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S2  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S2  Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  Al2 S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  Al2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  Al2 S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  Al2 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 S1  Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 S1  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al1 S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al1 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al1 S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 S2  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 S2  Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 S2  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 S2  Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al2 S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S1  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S1  Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S1  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S1  Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  Al1 S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  Al1 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  Al1 S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  Al1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S2  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S2  Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S2  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S2  Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  Al2 S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  Al2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  Al2 S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 S1  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 S1  Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 S1  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 S1  Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al1 S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 S2  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 S2  Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al2 S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al2 S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al2 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
