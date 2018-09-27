################################################################################
################################### IMPORT #####################################
################################################################################
import constants as c
from matplotlib import pyplot as plt
import numpy as np
import random as rand
import math as m

################################################################################
################################# LOADING ######################################
################################################################################
def load(name):
    return np.loadtxt(name, dtype = 'int')

################################################################################
################################# CALCUL #######################################
################################################################################
def photon_energy():
    # x = []
    # fx = []
    # for i in range(1000):
    l = rand.uniform(100e-9, 400e-9)
    energy = ((c.H * c.C) / l) / 1.60e-19
    #     x.append(l)
    #     fx.append(energy)
    # plt.plot(x, fx, '.')
    # plt.show()
    return energy


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
