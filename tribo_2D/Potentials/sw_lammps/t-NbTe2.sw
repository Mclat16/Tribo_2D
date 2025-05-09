# by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

# SW parameters for 1T-NbTe2  , used by LAMMPS
# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3,
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q,tol

# sw2 for bond Nb-Te and sw3 for angle Nb-Te-Te and Te-Nb-Nb
Nb  Te1 Te1    1.000   1.094   3.328  30.968   1.000   0.174   3.123  20.560  4  0   0.000
Te1 Nb  Nb     1.000   1.094   3.328  30.968   1.000   0.174   3.123  20.560  4  0   0.000
Nb  Te2 Te2    1.000   1.094   3.328  30.968   1.000   0.174   3.123  20.560  4  0   0.000
Te2 Nb  Nb     1.000   1.094   3.328  30.968   1.000   0.174   3.123  20.560  4  0   0.000

# zero terms
Nb  Nb  Nb     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Nb  Nb  Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Nb  Nb  Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Nb  Te1 Nb     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Nb  Te1 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Nb  Te2 Nb     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Nb  Te2 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Nb  Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Nb  Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te1 Nb     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te1 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te1 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te2 Nb     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te2 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te2 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Nb  Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Nb  Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te1 Nb     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te1 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te1 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te2 Nb     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te2 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te2 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
