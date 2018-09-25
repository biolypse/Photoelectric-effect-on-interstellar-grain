# -*- coding: utf-8 -*-
"""
@author : Nicolas
"""
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
matrix = f.load("Grain_N50_S0p1_B5p0.txt")
nb_rows, nb_column = matrix.shape
diameter = f.pixel_size_calculator(matrix)
grain_size = 2 * (c.GRAIN_MAX_RADIUS)
pixel_size = grain_size / diameter

################################################################################
################################## TREATMENT ###################################
################################################################################
with open("results.txt", "a") as file :
    for i in range(1000000):
        da = f.non_uniform_generator_exp(c.LA)
        de = f.non_uniform_generator_exp(c.LE)
        da_in_pixel = int(da / pixel_size)
        de_in_pixel = int(de / pixel_size)
        de_in_pixel_diag = int((de / (np.sqrt(2) * pixel_size)))

        energy = rand.uniform(3, 16)

        angles = [0, 45, 90, 135, 180, 225, 270, 315]
        ejection_direction = rand.choice(angles)
        ejection_energy = f.non_uniform_generator_th()

        photon_init_position = rand.randrange(0, nb_rows)
        row_photon = matrix[photon_init_position]
        contact_pixel, touch = f.contact(row_photon)

################################# PRINT TEST N°1 ###############################
        #print("La taille du grain est de : ", grain_size, "m")
        #print("Le diamètre du grain en pixel est de : ", diameter)
        #print("La taille des pixels est de : ", pixel_size, "m")
        #print("le photon arrive sur la ligne n°", photon_init_position)

        if touch :
            absorption_column, is_absorbed = f.absorption(da, pixel_size, grain_size, contact_pixel)
            if is_absorbed:
                is_ejected, ionization_energy = f.eject_electron(ejection_energy, grain_size)
                if is_ejected:
                    try :
                        is_free = f.freedom(ejection_direction, de, de_in_pixel, pixel_size, absorption_column, photon_init_position, matrix)
                    except IndexError:
                        pass
                    if is_free:
                        print("freeeeeeeedom")
                        kinetic_energy = energy - ionization_energy
                        file.write("{}\n".format(kinetic_energy))

################################################## PRINT TEST N°2 ####################################################
                        #print("Le photon rencontre le grain à la colonne n°", contact_pixel)
                        #print("Le photon est absorbé à la colonne n°", absorption_column)
                        #print("La distance parcouru par l'électron est de : ", de, "m sois ", de_in_pixel, "pixels")
                        #print("La distance hypothétiquement parcouru en diagonale en pixel est de", de_in_pixel_diag)
                        #print("L'électron est éjecté? True/False : ",is_ejected)
                        #print("L'électron est éjecté selon un angle de : ", ejection_direction, "°")

                    else :
                        stock = -3
                        file.write("{}\n".format(stock))
                else :
                    stock = -2
                    file.write("{}\n".format(stock))
            else :
                stock = -1
                file.write("{}\n".format(stock))
        else :
            stock = -1
            file.write("{}\n".format(stock))

f.histogramme("results.txt")
