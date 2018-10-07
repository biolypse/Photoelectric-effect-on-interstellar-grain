import numpy as np
C = 299792458
E0 = 8
H = 6.62607e-34
GRAIN_SIZE = 1e-7
LA = 1e-8
LE = 1e-9
Z_ELECTRON = -1
RHO_C = 1.9
M_C = 12.0107*1.66e-27
N_C = (np.pi * GRAIN_SIZE**2 * RHO_C) / (M_C)
IONIZATION = 4.4 + (6 + 0.5) * (25.1) / (N_C**(0.5))

COLOR = ["blue", "red", "yellow"]
