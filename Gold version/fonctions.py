################################################################################
################################### IMPORT #####################################
################################################################################
import constants as c
from matplotlib import pyplot as plt
import numpy as np
import random as rand
import math as m
import time
import sys

################################################################################
################################# LOADING ######################################
################################################################################
# def load(name):
#     return np.loadtxt("Grain_File/name", dtype = 'int')

################################################################################
################################# CALCUL #######################################
################################################################################
def non_uniform_generator_exp(para):
    #x = []
    #for i in range(100000):
    a = np.random.exponential(para)
        #x.append(a)
        #i += 1
    # plt.hist(x, bins = 1000)
    # plt.show()
    return a

################################################################################
################################# TREATMENT ####################################
################################################################################
def contact(row):
    if sum(row) == 0:
        touch = False
        return None, touch
    else:
        touch = True
        i = 0
        while row[i] != 1:
            i += 1
        # for i in range(len(row)):
        #     if row[i] == 1:
        return i, touch


def absorption(da, pixel_size, GRAIN_SIZE, contact_position, matrix, photon_init_position):
    if c.LA < GRAIN_SIZE:
        da_in_pixel = int(da / pixel_size)
        absorption_pixel = int(contact_position + da_in_pixel)
        if matrix[photon_init_position, absorption_pixel] == 1:
            is_absorbed = True
            return absorption_pixel, is_absorbed
        else :
            is_absorbed = False
            return absorption_pixel, is_absorbed
    else :
        is_absorbed = False
        return None, is_absorbed



# def freedom(de, angle, matrix, photon_init_position, absorption_column, pixel_size):
#     rad = m.radians(angle)
#     step = de / 10000
#     lenght = 0
#     hole = 0
#     dist = [0]
#     total_dist = 0
#     while lenght < de:
#         lenght += step
#         if 0 <= angle <90:
#             ligne = photon_init_position + int((np.sin(rad) * lenght) / pixel_size)
#             colonne = absorption_column - int((np.cos(rad) * lenght) / pixel_size)
#         elif 90 <= angle < 180:
#             ligne = photon_init_position + int((np.sin(rad) * lenght) / pixel_size)
#             colonne = absorption_column + int((np.cos(rad) * lenght) / pixel_size)
#
#         elif 180 <= angle < 270:
#             ligne = photon_init_position - (int((np.sin(rad) * lenght) / pixel_size))
#             colonne = int((np.cos(rad) * lenght) / pixel_size) + absorption_column
#         else :
#             ligne = photon_init_position - int((np.sin(rad) * lenght) / pixel_size)
#             colonne = absorption_column - int((np.cos(rad) * lenght) / pixel_size)
#         if matrix[ligne, colonne] == hole:
#             dist.append(lenght)
#             hole = (hole + 1) % 2
#
#     try :
#         for i in range(0, len(dist), 2):
#             total_dist += dist[i+1] - dist[i]
#     except :
#         pass
#     return total_dist




def freedom(de, angle, pixel_size, matrix, photon_init_position, absorption_column):
    rad = m.radians(angle)
    colonne = int((np.cos(rad) * de) / pixel_size)
    ligne = int((np.sin(rad) * de) / pixel_size)
    if 0 < angle <= 90:
        if matrix[photon_init_position - ligne, absorption_column + colonne] == 0:
            return True
        else :
            return False
    elif 90 < angle <=180:
        if matrix[photon_init_position - ligne, absorption_column - colonne] == 0:
            return True
        else :
            return False
    elif 180 < angle <= 270:
        if matrix[photon_init_position + ligne, absorption_column - colonne] ==0:
            return True
        else :
            return False
    elif 270 < angle <= 360:
        if matrix[photon_init_position + ligne, absorption_column + colonne] == 0:
            return True
        else :
            return False


################################################################################
################################# EXPLOITATION #################################
################################################################################
def histogramme(name):
    results = np.loadtxt(name)
    plt.hist(results, bins = 30, range = (0,13), color = "red", edgecolor = "black")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Number of photo_electron")
    plt.show()


################################################################################
################################################################################
################################# MAIN FUNCTION ################################
################################################################################

def main_function(number_of_grain, GRAIN_SIZE, sigma, S, p):
################################################################################
################################ INITIALIZATION ################################
################################################################################
    for i in range (number_of_grain):
        # if choice == "1" or choice == "4":
        matrix = np.loadtxt("Grain_Files/Grain_N100_S{}p{}_B3p0_{}.txt".format(S, p, i))
        dim1, dim2 = matrix.shape
        GRAIN_SIZE = c.GRAIN_SIZE * 2
        pixel_size = GRAIN_SIZE / dim1
    ################################################################################
    ################################## TREATMENT ###################################
    ################################################################################
        with open("results_Grain_N100_S{}p{}_B3p0_{}.txt".format(S, p, i), "a") as file :
            for i in range(1000000):
                j = i / 10000
                sys.stdout.write("\r%s%%" % j )
                sys.stdout.flush()
                da = non_uniform_generator_exp(c.LA)
                de = non_uniform_generator_exp(c.LE)
                da_in_pixel = int(da / pixel_size)
                de_in_pixel = int(de / pixel_size)
                de = de_in_pixel * pixel_size
                de_in_pixel_diag = int((de / (np.sqrt(2) * pixel_size)))
                photon_init_position = rand.randrange(0, dim1)

                row_photon = matrix[photon_init_position]

                angle = rand.randrange(0,360)

                contact_pixel, touch = contact(row_photon)

                if touch :
                    try :
                        absorption_column, is_absorbed = absorption(da, pixel_size, GRAIN_SIZE, contact_pixel, matrix, photon_init_position)
                    except IndexError:
                        pass

                    if is_absorbed :
                        energy = rand.uniform(3,15)
                        Y = 0.5 * (1 + np.tanh((energy - c.E0) / 2))
                        pile_face = rand.uniform(0,1)
                        if pile_face <= Y:
                            is_ejected = True
                        else :
                            is_ejected = False

                        if is_ejected :
                            try :
                                is_free = freedom(de, angle, pixel_size, matrix, photon_init_position, absorption_column)
                            except IndexError :
                                pass

                            if is_free :
                                kinetic_energy = energy - c.IONIZATION
                                file.write("{}\n".format(kinetic_energy))

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
