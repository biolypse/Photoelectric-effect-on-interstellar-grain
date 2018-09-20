################################################################################
################################### IMPORT #####################################
################################################################################
import fonctions as f
import constants as c
from matplotlib import pyplot as plt
import numpy as np
import random as rand
################################################################################
################################ INITIALIZATION ################################
################################################################################
angles = [0, 45, 90, 135, 180, 225, 270, 315]
ejection_direction = rand.choice(angles)
ejection_test = f.non_uniform_generator_th()
energy = rand.uniform(3, 16)
da = f.non_uniform_generator_exp(c.LA)
de = f.non_uniform_generator_exp(c.LE)
grain_size = 2 * (rand.uniform(c.GRAIN_MIN_RADIUS, c.GRAIN_MAX_RADIUS))
matrix = f.load("Grain_N500_S0p1_B5p0.txt")
nb_rows, nb_column = matrix.shape
pixel_size = grain_size / nb_rows





################################################################################
photon_init_position = rand.randrange(0, nb_rows)
row_photon = matrix[photon_init_position]


################################################################################
############################# FIRST PRINT TESTS ################################
################################################################################
#print("La taille du grain est de : ", grain_size, "m")
#print("La taille des pixels est de : ", pixel_size)
#print("le photon arrive sur la ligne n°", photon_init_position)

################################################################################
contact_pixel, exit = f.contact(row_photon)

if not exit :
    absorption_pixel = f.absorption(da, pixel_size, grain_size, contact_pixel)
    is_ejected = f.eject_electron(ejection_test, grain_size)


################################################################################
############################# SECOND PRINT TESTS ###############################
################################################################################
    #print("Le photon est absorbé à la colonne n°", absorption_pixel)
    #print("distance parcouru par l'électron", de)
    #print("L'électron est éjecté? ",is_ejected)
