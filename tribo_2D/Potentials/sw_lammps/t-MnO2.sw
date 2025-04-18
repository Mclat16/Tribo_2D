# by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

# SW parameters for 1T-MnO2   , used by LAMMPS
# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3,
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q,tol

# sw2 for bond Mn-O  and sw3 for angle Mn-O -O  and O -Mn-Mn
Mn  O1  O1     1.000   1.212   2.175  60.755   1.000  -0.125   9.675   2.899  4  0   0.000
O1  Mn  Mn     1.000   1.212   2.175  60.755   1.000  -0.125   9.675   2.899  4  0   0.000
Mn  O2  O2     1.000   1.212   2.175  60.755   1.000  -0.125   9.675   2.899  4  0   0.000
O2  Mn  Mn     1.000   1.212   2.175  60.755   1.000  -0.125   9.675   2.899  4  0   0.000

# zero terms
Mn  Mn  Mn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Mn  Mn  O1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Mn  Mn  O2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Mn  O1  Mn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Mn  O1  O2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Mn  O2  Mn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Mn  O2  O1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O1  Mn  O1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O1  Mn  O2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O1  O1  Mn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O1  O1  O1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O1  O1  O2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O1  O2  Mn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O1  O2  O1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O1  O2  O2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O2  Mn  O1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O2  Mn  O2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O2  O1  Mn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O2  O1  O1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O2  O1  O2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O2  O2  Mn     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O2  O2  O1     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
O2  O2  O2     0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
