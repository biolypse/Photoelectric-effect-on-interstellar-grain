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


def pixel_size_calculator(matrix):
    diameter = 0
    for line in matrix:
        if sum(line) > diameter:
            diameter = sum(line)
    return diameter


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
        touch = False
        print("Le photon n'a pas rencontré le grain")
        return None, touch
    else:
        touch = True
        i = 0
        while row[i] != 1:
            i += 1
        # for i in range(len(row)):
        #     if row[i] == 1:
        return i, touch


def absorption(da, pixel_size, grain_size, contact_position):
    if c.LA < grain_size:
        is_absorbed = True
        da_in_pixel = int(da / pixel_size)
        absorption_pixel = int(contact_position + da_in_pixel)
        return absorption_pixel, is_absorbed
    else :
        is_absorbed = False
        print("Il n'y a pas eu d'absorption")
        return None, is_absorbed


def eject_electron(ejection_energy, grain_size):
    nc = (np.pi * (grain_size ** 2) * c.RHO_C) / (c.M_C)
    gap = 4.4 + (c.Z_ELECTRON + 0.5) * (25.1) / (nc**(1 / 2))
    if ejection_energy > gap:
        ejection = True
    else:
        ejection = False
    return ejection, gap


def freedom(ejection_direction, de, de_in_pixel, pixel_size, absorption_column, photon_init_position, matrix):
    de_in_pixel_diag = int((de / (np.sqrt(2) * pixel_size)))
    if ejection_direction == 0:
        if matrix[photon_init_position, absorption_column + de_in_pixel] == 1:
            return False
        else :
            return True
    elif ejection_direction == 90:
        if matrix[photon_init_position - de_in_pixel, absorption_column] == 1:
            return False
        else :
            return True
    elif ejection_direction == 180:
        if matrix[photon_init_position, absorption_column - de_in_pixel] == 1:
            return False
        else :
            return True
    elif ejection_direction == 270:
        if matrix[photon_init_position + de_in_pixel, absorption_column] == 1:
            return False
        else :
            return True
    elif ejection_direction == 45:
        if matrix[photon_init_position - de_in_pixel_diag, absorption_column + de_in_pixel_diag] == 1:
            return False
        else :
            return True
    elif ejection_direction == 135:
        if matrix[photon_init_position - de_in_pixel_diag, absorption_column - de_in_pixel_diag] == 1:
            return False
        else :
            return True
    elif ejection_direction == 225:
        if matrix[photon_init_position + de_in_pixel_diag, absorption_column - de_in_pixel_diag] == 1:
            return False
        else :
            return True
    elif ejection_direction == 315:
        if matrix[photon_init_position + de_in_pixel_diag, absorption_column + de_in_pixel_diag] == 1:
            return False
        else :
            return True
