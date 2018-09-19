################################################################################
################################### IMPORT #####################################
################################################################################
import constants as c
from matplotlib import pyplot as plt
import numpy as np
from random import uniform
################################################################################
################################# FUNCTIONS ####################################
################################################################################
def load(name):
    return np.loadtxt(name)


def non_uniform_generator_exp():
    #x = []
    #for i in range(100000):
    a = (np.random.exponential(c.LA))
        #x.append(a)
        #i += 1
    # plt.hist(x, bins = 1000)
    # plt.show()
    return a


def non_uniform_generator_th():
    # x = []
    # fx = []
    for i in range(1000):
        a = uniform(3, 16)
        u = uniform(0, 1)
        while u > 0.5 * (1 + np.tanh((a - c.E0) / 2)):
            a = uniform(3, 16)
            u = uniform(0, 1)
        return a
        # b = 0.5 * (1 + np.tanh((a - 8) / 2))
        # x.append(a)
        # fx.append(b)

    # plt.plot(x, fx, '.')
    # plt.show()


def x_contact(row):
    if sum(row) == 0:
        print("Le photon n'a pas rencontré le grain")
        exit = True
        return None, exit
    else:
        exit = False
        i = 0
        while row[i] != 1:
            i += 1
        # for i in range(len(row)):
        #     if row[i] == 1:
        return i, exit


def x_absorption(da, nb_rows, grain_size, contact_position):
    pixel_size = grain_size / len(nb_rows)
    print("taille d'un pixel", pixel_size)
    if c.LA < grain_size:
        absorption = np.round_(da / pixel_size)
        absorption_pixel = int(contact_position + absorption)
        return(absorption_pixel)
    else :
        print("Il n'y a pas d'absorption")


# def eject_electron(ejection_test, grain_size):
#     ejection = True
#     nc = ((4 / 3) * np.pi * (grain_size ** 3) * c.RHO_C) / (c.M_C)
#     gap = 4.4 + (c.Z_ELECTRON + 0.5) * (11.1) / (nc**(1 / 3))
#     if ejection_test > nc:
#         ejection = True
#     else:
#         ejection = False
#     return ejection
