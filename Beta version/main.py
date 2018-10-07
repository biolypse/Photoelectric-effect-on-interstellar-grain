# -*- coding: utf-8 -*-
"""
@author : Nico
"""
################################################################################
################################### IMPORT #####################################
################################################################################
import fonctions as f
import random as rand
import constants as c
from matplotlib import pyplot as plt
import numpy as np

################################################################################
################################ INITIALIZATION ################################
################################################################################
matrix = f.load("Grain_N100_S1p0_B3p0.txt")
dim1, dim2 = matrix.shape
GRAIN_SIZE = c.GRAIN_SIZE * 2
pixel_size = GRAIN_SIZE / dim1
touche = 0
absorption = 0
ejection = 0
libere = 0

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
        photon_init_position = rand.randrange(0, dim1)
        row_photon = matrix[photon_init_position]
        angle = rand.randrange(0,360)

        contact_pixel, touch = f.contact(row_photon)

################################################################ PRINT TEST N°1 ####################################################
        #print("Le photon arrive sur la ligne", photon_init_position, "Le photon rencontre le grain à la colonne n°", contact_pixel)
        #print("La taille des pixels est de : ", pixel_size, "m")
        #print("Le photon parcours", da,"m")
        #print("En terme de pixel", da_in_pixel)

        if touch :
            touche += 1
            try :
                absorption_column, is_absorbed = f.absorption(da, pixel_size, GRAIN_SIZE, contact_pixel, matrix, photon_init_position)
            except IndexError:
                pass
            #print("Le photon a touché")
            if is_absorbed :
                absorption += 1
                #print("Le photon est absorbé")
                energy = rand.uniform(3,15)
                Y = 0.5 * (1 + np.tanh((energy - c.E0) / 2))
                pile_face = rand.uniform(0,1)
                if pile_face <= Y:
                    is_ejected = True
                else :
                    is_ejected = False

                if is_ejected :
                    ejection += 1
                    #print("L'électron est éjecté")
                    try :
                        is_free = f.freedom(de, angle, pixel_size, matrix, photon_init_position, absorption_column)
                    except IndexError :
                        pass
                    if is_free :
                        libere += 1
                        #print("FREEEEEEEEEEEEEEEDOM")
                        kinetic_energy = energy - c.IONIZATION
                        file.write("{}\n".format(kinetic_energy))

                    else :
                        stock = -3
                        file.write("{}\n".format(stock))
                else :
                    # print("L'électron n'est pas éjecté")
                    stock = -2
                    file.write("{}\n".format(stock))
            else :
                # print("Le photon n'a pas été absorbé")
                stock = -1
                file.write("{}\n".format(stock))
        else :
            # print("Le photon n'a pas touché")
            stock = -1
            file.write("{}\n".format(stock))

################################## PRINT TEST N°2###############################
print("Nombres de photon qui ont touché le grain : ", touche)
print("Nombres de photon qui ont été absorbé : ", absorption)
print("Nombres d'électron qui ont été créé : ", ejection)
print("Nombres d'électron qui sont sortis du grain", libere)
f.histogramme("results.txt")
