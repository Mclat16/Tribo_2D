# by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

# (1). SW parameters for t-ZrS2, used by LAMMPS, derived from VFF model analytically.

# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3, 
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q, tol

# sw2, Zr-S bond stretching; and sw3, Zr-S-S angle bending
Zr  S1  S1     1.000   1.432   2.473  42.177   1.000  -0.023   8.149   5.268  4  0   0.000
Zr  S2  S2     1.000   1.432   2.473  42.177   1.000  -0.023   8.149   5.268  4  0   0.000

# sw2, Zr-S bond stretching; and sw3, S-Zr-Zr angle bending
S1  Zr  Zr     1.000   1.432   2.473  42.177   1.000  -0.023   8.149   5.268  4  0   0.000
S2  Zr  Zr     1.000   1.432   2.473  42.177   1.000  -0.023   8.149   5.268  4  0   0.000

# zero terms
Zr  Zr  Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  Zr  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  Zr  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  S1  Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  S1  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  S2  Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Zr  S2  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  Zr  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  Zr  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S1  Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S1  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S1  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S2  Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S2  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S1  S2  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  Zr  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  Zr  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S1  Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S1  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S1  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S2  Zr     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S2  S1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
S2  S2  S2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
