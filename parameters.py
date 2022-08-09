# import NumPy package
import numpy as np

def parameters:
  # output the 32 parameters for the aspartate/beta alanine pathway (see 'parameters' sect of paper)

  # small dimensionless parameter epsilon
  ep = 1e-2

  # use value s_0 = 1mM
  s_0 = 1e-3

  # initial conditions. first six are substrates, then the two reactants, then four enzymes and complexes
  # ordering is (S1,...S6,R1,R2,E1,C1,...,E4,C4,P)
  y0 = np.zeros((8, 1), dtype=float)
  y0[0] = s_0

  # reaction rates (see table 2 of paper)
  k = np.zeros((20, 1), dtype=float)

  # k1 saccharomyces cerevisiae
  k[0] = 10
  k1 = k[1]

  # k2 enterococcus faecalis
  k[1] = k1 * 0.2 / ep

  # k_2 enterococcus faecalis
  k[2] = k1 * 3 / ep

  # k3 saccharomyces cerevisiae
  k[3] = k1 * 0.6

  # k4 enterococcus faecalis
  k[4] = k1 * 1

  # k5 methanosarcina mazei
  k[5] = k1 * 2

  # k6 saccharomyces cerevisiae
  k[6] = k1 * 0.4

  # k7 saccharomyces cerevisiae
  k[7] = 2 * k1 * 0.1

  # k_1^M saccharomyces cerevisiae
  k[8] = 0.65

  # k_1^i saccharomyces cerevisiae
  k[9] = 1.4 * ep

  # k_2^M enterococcus faecalis
  k[10] = 1

  # k_{-2}^M enterococcus faecalis
  k[11] = ep

  # k_{3,a}^M saccharomyces cerevisiae
  k[12] = 1.5 * ep

  # k_{3,b}^M saccharomyces cerevisiae
  k[13] = 0.3 * ep

  # k_3^i saccharomyces cerevisiae
  k[14] = ep

  # k_4^M enterococcus faecalis
  k[15] = 2 * ep

  # k_5^M methanosarcina mazei
  k[16] = 0.1

  # k_6^M saccharomyces cerevisiae
  k[17] = 0.9

  # k_7^M saccharomyces cerevisiae
  k[18] = 0.2

  # A
  k[19] = 1

  k[8:18] = s_0 * k[8:18]
  
  return k;
  return y0;
  return ep;
