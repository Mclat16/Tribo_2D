# by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

# (1). SW parameters for t-ZrSe2, used by LAMMPS, derived from VFF model analytically.

# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3, 
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q, tol

# sw2, Zr-Se bond stretching; and sw3, Zr-Se-Se angle bending
Zr  Se1 Se1    1.000   1.354   2.671  37.051   1.000   0.034   8.022   7.527  4  0   0.000
Zr  Se2 Se2    1.000   1.354   2.671  37.051   1.000   0.034   8.022   7.527  4  0   0.000

# sw2, Zr-Se bond stretching; and sw3, Se-Zr-Zr angle bending
Se1 Zr  Zr     1.000   1.354   2.671  37.051   1.000   0.034   8.022   7.527  4  0   0.000
Se2 Zr  Zr     1.000   1.354   2.671  37.051   1.000   0.034   8.022   7.527  4  0   0.000

# zero terms
Zr  Zr  Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  Zr  Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  Zr  Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  Se1 Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  Se2 Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Zr  Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Zr  Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se1 Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se1 Se2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Zr  Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Zr  Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se1 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Se1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Se2 Se2 Se2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
