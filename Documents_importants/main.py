import fonctions as f
from matplotlib import pyplot as plt
import numpy as np
from random import randrange

matrix = f.load("Grain_N10_S2p0_B5p0.txt")
nb_rows, nb_column = matrix.shape

x_init = randrange(0, nb_rows)

row = matrix[x_init]

energy = randrange(3, 16)
