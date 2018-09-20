################################################################################
################################### IMPORT #####################################
################################################################################
import constants as c
from matplotlib import pyplot as plt
import numpy as np
from random import uniform
################################################################################
################################# LOADING ######################################
################################################################################
def load(name):
    return np.loadtxt(name, dtype = 'int')

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

################################################################################
################################# TREATMENT ####################################
################################################################################
def contact(row):
    if sum(row) == 0:
        #print("Le photon n'a pas rencontré le grain")
        exit = True
        print("Le photon n'a pas rencontré le grain")
        return None, exit
    else:
        exit = False
        i = 0
        while row[i] != 1:
            i += 1
        # for i in range(len(row)):
        #     if row[i] == 1:
        return i, exit


def absorption(da, pixel_size, grain_size, contact_position):
    if c.LA < grain_size:
        da_in_pixel = np.round_(da / pixel_size)
        absorption_pixel = int(contact_position + da_in_pixel)
        return absorption_pixel
    else :
        print("Il n'y a pas eu d'absorption")


def eject_electron(ejection_test, grain_size):
    nc = (np.pi * (grain_size ** 2) * c.RHO_C) / (c.M_C)
    gap = 4.4 + (c.Z_ELECTRON + 0.5) * (25.1) / (nc**(1 / 2))
    if ejection_test > gap:
        ejection = True
    else:
        ejection = False
    return ejection
