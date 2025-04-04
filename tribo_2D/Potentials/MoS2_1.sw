#
# Format of parameter file:
#
# First line: number of species
#
# Lines 2~4: parameters for 2-body interaction
#   species1 species2 A B p q sigma gamma cutoff
#
# Lines 5~6: parameters for 3-body interaction
#   species1 species2 species3 lambda cos_beta0 cutoff_jk
#
# species is valid KIM API particle species string
# A and lambda in [eV]
# sigma, gamma, cutoff, and cutoff_jk in [Angstrom]
# others are unitless
#

#
# Parameters are taken from:
#
# M. Wen, S. N. Shirodkar, P. Plechac, E. Kaxiras, and E. B. Tadmor, "Stillinger-Weber
# potential for MoS2: Parameterization and sensitivity analysis," J. Appl. Phys.,
# 122, 244301 (2017).
#

2

Mo  Mo   3.9781804791   0.4446021306  5 0  2.85295  1.3566322033  5.54660
Mo  S   11.3797414404   0.5266688197  5 0  2.17517  1.3566322033  4.02692
S   S    1.1907355764   0.9015152673  5 0  2.84133  1.3566322033  4.51956

Mo  S   S   7.4767529158  0.1428569579923222  3.86095
S   Mo  Mo  8.1595181220  0.1428569579923222  5.54660


