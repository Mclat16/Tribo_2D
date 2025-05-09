# by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

# SW parameters for AlSe      , used by LAMMPS
# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3,
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q,tol

# sw2 and sw3, for Al-Se-Se and Se-Al-Al
Al1 Se1 Se1    1.000   1.694   2.062  41.235   1.000  -0.171  11.362   2.260  4  0   0.000
Se1 Al1 Al1    1.000   1.694   2.062  41.235   1.000  -0.171  11.362   2.260  4  0   0.000
Al2 Se2 Se2    1.000   1.694   2.062  41.235   1.000  -0.171  11.362   2.260  4  0   0.000
Se2 Al2 Al2    1.000   1.694   2.062  41.235   1.000  -0.171  11.362   2.260  4  0   0.000

# sw2, for Al-Al
Al1 Al2 Al2    1.000   1.558   2.292   0.000   1.000   0.000   4.974   3.704  4  0   0.000
Al2 Al1 Al1    1.000   1.558   2.292   0.000   1.000   0.000   4.974   3.704  4  0   0.000

# sw3 for Al-Al-Se
Al1 Al2 Se1    1.000   0.000   0.000  26.418   1.000  -0.468   0.000   0.000  4  0   0.000
Al1 Se1 Al2    1.000   0.000   0.000  26.418   1.000  -0.468   0.000   0.000  4  0   0.000
Al2 Al1 Se2    1.000   0.000   0.000  26.418   1.000  -0.468   0.000   0.000  4  0   0.000
Al2 Se2 Al1    1.000   0.000   0.000  26.418   1.000  -0.468   0.000   0.000  4  0   0.000

# zero terms
Se1 Se1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se1 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Al1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Al1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Al1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Al2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Al2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Al2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Al2 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Se1 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al1 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Se2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Se2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Se2 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al1 Al2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Al1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Al1 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Al1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Al1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Al2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Al2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Al2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Se1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Se1 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Se1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al1 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Se2 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al2 Al1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Al2 Al2 Al2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
