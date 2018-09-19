################################################################################
################################### IMPORT #####################################
################################################################################
import fonctions as f
import constants as c
from matplotlib import pyplot as plt
import numpy as np
from random import randrange
from random import uniform

################################################################################
################################ INITIALIZATION ################################
################################################################################
energy = uniform(3, 16)
da = f.non_uniform_generator_exp()
grain_size = uniform(c.GRAIN_MIN_RADIUS, c.GRAIN_MAX_RADIUS)




################################################################################
matrix = f.load("Grain_N50_S0p1_B5p0.txt")
nb_rows, nb_column = matrix.shape


photon_init_position = randrange(0, nb_rows)
row = matrix[photon_init_position]

################################################################################
################################ PRINT TEST ####################################
################################################################################
print("La taille du grain est de : ", grain_size, "m")
print("le photon arrive sur la ligne n°", photon_init_position)


contact_position, exit = f.x_contact(row)
if not exit :
    print("le photon rencontre le grain en position", contact_position)
    print("le photon est absorbé au pixel ",f.x_absorption(da, matrix,
          grain_size, contact_position))

    # ejection_test = f.non_uniform_generator_th()
    # ejection = f.eject_electron(ejection_test, grain_size)
    # print(ejection)
