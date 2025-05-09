# by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

# SW parameters for 1T-PtTe2  , used by LAMMPS
# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3,
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q,tol

# sw2 for bond Pt-Te and sw3 for angle Pt-Te-Te and Te-Pt-Pt
Pt  Te1 Te1    1.000   1.667   2.229  59.607   1.000  -0.104  14.877   3.250  4  0   0.000
Te1 Pt  Pt     1.000   1.667   2.229  59.607   1.000  -0.104  14.877   3.250  4  0   0.000
Pt  Te2 Te2    1.000   1.667   2.229  59.607   1.000  -0.104  14.877   3.250  4  0   0.000
Te2 Pt  Pt     1.000   1.667   2.229  59.607   1.000  -0.104  14.877   3.250  4  0   0.000

# zero terms
Pt  Pt  Pt     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Pt  Pt  Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Pt  Pt  Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Pt  Te1 Pt     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Pt  Te1 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Pt  Te2 Pt     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Pt  Te2 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Pt  Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Pt  Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te1 Pt     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te1 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te1 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te2 Pt     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te2 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te2 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Pt  Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Pt  Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te1 Pt     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te1 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te1 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te2 Pt     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te2 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te2 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
