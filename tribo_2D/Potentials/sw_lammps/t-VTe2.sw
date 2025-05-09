# by Jin-Wu Jiang, jwjiang5918@hotmail.com; jiangjinwu@shu.edu.cn

# SW parameters for 1T-VTe2   , used by LAMMPS
# these entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless

# format of a single entry (one or more lines):
#   element 1, element 2, element 3,
#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q,tol

# sw2 for bond V -Te and sw3 for angle V -Te-Te and Te-V -V 
V   Te1 Te1    1.000   1.110   3.149  26.043   1.000   0.141   8.805  15.980  4  0   0.000
Te1 V   V      1.000   1.110   3.149  26.043   1.000   0.141   8.805  15.980  4  0   0.000
V   Te2 Te2    1.000   1.110   3.149  26.043   1.000   0.141   8.805  15.980  4  0   0.000
Te2 V   V      1.000   1.110   3.149  26.043   1.000   0.141   8.805  15.980  4  0   0.000

# zero terms
V   V   V      0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
V   V   Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
V   V   Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
V   Te1 V      0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
V   Te1 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
V   Te2 V      0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
V   Te2 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 V   Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 V   Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te1 V      0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te1 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te1 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te2 V      0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te2 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te1 Te2 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 V   Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 V   Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te1 V      0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te1 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te1 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te2 V      0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te2 Te1    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
Te2 Te2 Te2    0.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000  4  0   0.000
