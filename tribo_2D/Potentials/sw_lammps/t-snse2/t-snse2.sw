# by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

# SW parameters for 1T-SnSe2  , used by LAMMPS
# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3,
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q,tol

# sw2 for bond Sn-Se and sw3 for angle Sn-Se-Se and Se-Sn-Sn
Sn  Se1 Se1    1.000   1.411   2.609  52.322   1.000   0.017   8.395   6.743  4  0   0.000
Se1 Sn  Sn     1.000   1.411   2.609  52.322   1.000   0.017   8.395   6.743  4  0   0.000
Sn  Se2 Se2    1.000   1.411   2.609  52.322   1.000   0.017   8.395   6.743  4  0   0.000
Se2 Sn  Sn     1.000   1.411   2.609  52.322   1.000   0.017   8.395   6.743  4  0   0.000

# zero terms
Sn  Sn  Sn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Sn  Sn  Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Sn  Sn  Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Sn  Se1 Sn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Sn  Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Sn  Se2 Sn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Sn  Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Sn  Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Sn  Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se1 Sn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Sn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Sn  Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Sn  Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Sn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Sn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
