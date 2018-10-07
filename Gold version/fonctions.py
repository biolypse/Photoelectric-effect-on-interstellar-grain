# -*- coding: utf-8 -*-
"""
@author : Nicolas
"""

import constants as c
from matplotlib import pyplot as plt
import numpy as np
import random as rand
import math as m
from GRF_routines import addGRF
import time
import sys

################################################################################
################################# GRAIN GENERATOR ##############################
################################################################################
def GenGrain(sigma_dens, number_of_grain=0, test = False):
    # Tunable parameters
    # - Main parameters
    N = 100           #size of the grain image
    # sigma_dens = 0.9     Width of the density distribution
    beta = 3.0           # Slope of the power spectrum (=probability density function=PDF) for k>Kmin

    # - Secondary parameters
    Kmin = 3             # Minimum wavenumber for PDF slope = beta
    beta_in = 4          # Slope of PDF between k=0 and Kmin
    Npad = 30            # Padding around the map to get rid of the peak near k=0
    Adens = 1.           # Mean density - should not impact dust grain

    doplot = 0           # Visual check of grain calculation: 0= desactivated, 1= only final plot, 2= animation of search for the main group + final plot


    # while loop to make sure the peak is sufficiently close to the center of the cube
    tol = 0.4            # Tolerance for location of the density maximum - between 0 (very tolerant) and 0.49999 (very tough)
    mmax = [0,0]
    ii = 0
    while ((mmax[0]<N*tol) | (mmax[0]>N*(1-tol)) | (mmax[1]<N*tol) | (mmax[1]>N*(1-tol))):

       # Build a fractal cube using Gaussian Random Field
       cube = np.zeros((N,N))
       cube = addGRF(cube, sigma_dens, Adens, beta, Npad, Kmin, beta_in)

       # Search the location of the max value
       mmax = np.nonzero(cube==np.amax(cube))
       # print(mmax)
       ii += 1
       if (ii>1000):
          print("Failed generating a GRF with this tolerance value:", tol)
          print("Try a lower value. It should be between 0 (very tolerant) and 0.49999 (very tough).")
          exit()


    # Roll the cube to put the max at the center
    cube = np.roll(cube, N//2-mmax[0][0], axis=0)
    cube = np.roll(cube, N//2-mmax[1][0], axis=1)

    # Apodize the cube
    I, J = np.indices(np.shape(cube))
    Rad = (I-N/2)**2+(J-N/2)**2
    sig = N*0.1
    cube *= np.exp(-Rad/(2*sig**2))

    # Compute the location of the mass center
    cdm = [np.sum(J*cube)/np.sum(cube),np.sum(I*cube)/np.sum(cube)]
    # print(cdm)

    # Roll the cube to put the mass center at the center
    cube = np.roll(cube, int(N/2-cdm[0]), axis=0)
    cube = np.roll(cube, int(N/2-cdm[1]), axis=1)

    # Threshold
    grain = np.zeros(np.shape(cube))
    m = np.nonzero(cube>np.percentile(cube,70))
    grain[m] = 1.

    # Identify the pixels attached to the main region
    def group(grain, x0,y0):
       # Warning: no precaution for edges
       goon = 1
       progress = np.zeros(np.shape(grain))
       progress[x0,y0] = 2
       if (doplot==2): plt.ion()
       while (goon==1):
          goon = 0
          m = np.nonzero(progress==2)   # search last included friends
          for i in range(len(m[0])):
             i1, j1 = m[0][i], m[1][i]
             # Search only touching pixels
             if (0):
                searchind = [[i1+1,j1], [i1+1,j1-1], [i1,j1-1], [i1-1,j1-1], \
                             [i1-1,j1], [i1-1,j1+1], [i1,j1+1], [i1+1,j1+1]]
             # Search touching pixels + 1
             if (1):
                searchind = [[i1+1,j1], [i1+1,j1-1], [i1,j1-1], [i1-1,j1-1], \
                             [i1-1,j1], [i1-1,j1+1], [i1,j1+1], [i1+1,j1+1], \
                             [i1+2,j1], [i1+2,j1-1], [i1+2,j1-2], [i1+1,j1-2], \
                             [i1,j1-2], [i1-1,j1-2], [i1-2,j1-2], [i1-2,j1-1], \
                             [i1-2,j1], [i1-2,j1+1], [i1-2,j1+2], [i1-1,j1+2], \
                             [i1,j1+2], [i1+1,j1+2], [i1+2,j1+2], [i1+2,j1+1]]
             for inds in searchind:
                inds = np.asarray(inds)
                if ((np.all(inds>=0)) & (np.all(inds<N))):
                   if (grain[inds[0],inds[1]]==1):
                      if (progress[inds[0],inds[1]]==0): # not found before
                         progress[inds[0],inds[1]] = 3   # 3 for freshly identified
                         goon = 1
          if (doplot==2): plt.imshow(progress, origin='low', cmap='bone', \
                                    interpolation='nearest')
          if (doplot==2): plt.pause(0.0001)
          progress[m] = 1
          progress[np.nonzero(progress==3)] = 2

       progress[np.nonzero(progress>1)] = 1
       plt.ioff()
       return progress

    mmax = np.nonzero(cube==np.amax(cube))
    grain2 = group(grain, mmax[0][0], mmax[1][0])

    # Write down the grain image in ASCII format (very inefficient format, but easy to control)
    if test:
        np.savetxt("Grain_Files/Grain_N%i_S%ip%i_B%ip%i_%i.txt" % (N, int(sigma_dens), \
                                                   int((sigma_dens-int(sigma_dens))*10), \
                                                   int(beta), \
                                                   int((beta-int(beta))*10),number_of_grain), \
                                                   grain2, fmt='%i')
    else:
        np.savetxt("Grain_Files/Grain_N%i_S%ip%i_B%ip%i.txt" % (N, int(sigma_dens), \
                                                       int((sigma_dens-int(sigma_dens))*10), \
                                                       int(beta), \
                                                       int((beta-int(beta))*10)), \
                                                       grain2, fmt='%i')


################################################################################
############################## CALCUL FUNCTION #################################
################################################################################
def non_uniform_generator_exp(parameter):
    """
    Generate a random number following an exponential law of probability

    Args :
        parameter(float) : parameter for the exponential law
    Return :
        a(float) : random number
    """
    a = np.random.exponential(parameter)
    return a

################################################################################
############################ TREATMENT FUNCTION ################################
################################################################################
def contact(row):
    """
    Verifie if the photon arrive on the grain

    Args :
        row(int) : photon input line
    Returns :
        touche(bool) : True if the photon has touch
        i(int) : column where the photon touch the grain
    """
    if sum(row) == 0:
        touch = False
        return None, touch
    else:
        touch = True
        i = 0
        while row[i] != 1:
            i += 1
        return i, touch


def absorption(da, pixel_size, GRAIN_SIZE, contact_position, matrix, photon_init_position):
    """
    Find the position where the photon is absorbe in the Grain_N

    Args :
        da(float) : distance traveled by the photon before being absorbe
        pixel_size(float) : size of the pixels
        GRAIN_SIZE(float) : size of the interstellar grain
        contact_position(int) : column where the photon touch the grain
        photon_init_position(int) : photon input line

    Returns :
        is_absorbed(bool) : True if the phton is absorbed by the grain
        absorption_pixel(int) : column where the photon is absorbed
    """
    if c.LA < GRAIN_SIZE:
        da_in_pixel = int(da / pixel_size)
        absorption_pixel = int(contact_position + da_in_pixel)
        if matrix[photon_init_position, absorption_pixel] == 1:
            is_absorbed = True
            return absorption_pixel, is_absorbed
        else:
            is_absorbed = False
            return absorption_pixel, is_absorbed
    else:
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
    """

    """
    rad = m.radians(angle)
    colonne = int((np.cos(rad) * de) / pixel_size)
    ligne = int((np.sin(rad) * de) / pixel_size)
    if 0 < angle <= 90:
        if matrix[photon_init_position - ligne, absorption_column + colonne] == 0:
            return True
        else:
            return False
    elif 90 < angle <= 180:
        if matrix[photon_init_position - ligne, absorption_column - colonne] == 0:
            return True
        else:
            return False
    elif 180 < angle <= 270:
        if matrix[photon_init_position + ligne, absorption_column - colonne] == 0:
            return True
        else:
            return False
    elif 270 < angle <= 360:
        if matrix[photon_init_position + ligne, absorption_column + colonne] == 0:
            return True
        else:
            return False

################################################################################
################################################################################
################################# MAIN FUNCTION ################################
################################################################################
################################################################################
def main_function(number_of_grain, GRAIN_RADIUS, sigma, choice):
    choice = int(choice)
    S = []
    p = []
    for j in range(number_of_grain):
        if choice == 1 or choice == 4:
            S = int(sigma)
            wait = str(sigma)
            wait = wait.split(".")
            p = wait[1]
            GenGrain(sigma, j, True)
            matrix = np.loadtxt("Grain_Files/Grain_N100_S{}p{}_B3p0_{}.txt".format(S, p, j))
        else:
            S.append(int(sigma[j]))
            wait = str(sigma[j])
            wait = wait.split(".")
            p.append(wait[1])
            GenGrain(sigma[j])
            matrix = np.loadtxt("Grain_Files/Grain_N100_S{}p{}_B3p0.txt".format(S[j], p[j]))

        dim1, dim2 = matrix.shape

        if choice == 1 or choice == 4:
            if choice == 4:
                GRAIN_SIZE = GRAIN_RADIUS[j] * 2
            else:
                GRAIN_SIZE = c.GRAIN_SIZE * 2

            pixel_size = GRAIN_SIZE / dim1
            with open("results_Grain_N100_S{}p{}_B3p0_{}.txt".format(S, p, j), "w") as file:
                for i in range(1000000):
                    da = non_uniform_generator_exp(c.LA)
                    de = non_uniform_generator_exp(c.LE)
                    photon_init_position = rand.randrange(0, dim1)
                    row_photon = matrix[photon_init_position]
                    angle = rand.randrange(0, 360)
                    contact_pixel, touch = contact(row_photon)

                    if touch:
                        try:
                            absorption_column, is_absorbed = absorption(da, pixel_size, GRAIN_SIZE, contact_pixel, matrix, photon_init_position)
                        except IndexError:
                            pass

                        if is_absorbed:
                            energy = rand.uniform(3, 15)
                            Y = 0.5 * (1 + np.tanh((energy - c.E0) / 2))
                            pile_face = rand.uniform(0, 1)

                            if pile_face <= Y:
                                is_ejected = True
                            else:
                                is_ejected = False

                            if is_ejected:
                                try:
                                    is_free = freedom(de, angle, pixel_size, matrix, photon_init_position, absorption_column)
                                except IndexError:
                                    pass

                                if is_free:
                                    kinetic_energy = energy - c.IONIZATION
                                    file.write("{}\n".format(kinetic_energy))

                                else:
                                    stock = -3
                                    file.write("{}\n".format(stock))
                            else:
                                stock = -2
                                file.write("{}\n".format(stock))
                        else:
                            stock = -1
                            file.write("{}\n".format(stock))
                    else:
                        stock = -1
                        file.write("{}\n".format(stock))
        else:
            GRAIN_SIZE = c.GRAIN_SIZE * 2
            pixel_size = GRAIN_SIZE / dim1
            with open("results_Grain_N100_S{}p{}_B3p0.txt".format(S[j], p[j]), "w") as file:
                for i in range(1000000):
                    da = non_uniform_generator_exp(c.LA)
                    de = non_uniform_generator_exp(c.LE)
                    photon_init_position = rand.randrange(0, dim1)
                    row_photon = matrix[photon_init_position]
                    angle = rand.randrange(0, 360)
                    contact_pixel, touch = contact(row_photon)

                    if touch:
                        try:
                            absorption_column, is_absorbed = absorption(da, pixel_size, GRAIN_SIZE, contact_pixel, matrix, photon_init_position)
                        except IndexError:
                            pass

                        if is_absorbed:
                            energy = rand.uniform(3, 15)
                            Y = 0.5 * (1 + np.tanh((energy - c.E0) / 2))
                            pile_face = rand.uniform(0, 1)

                            if pile_face <= Y:
                                is_ejected = True
                            else:
                                is_ejected = False

                            if is_ejected:
                                try:
                                    is_free = freedom(de, angle, pixel_size, matrix, photon_init_position, absorption_column)
                                except IndexError:
                                    pass

                                if is_free:
                                    kinetic_energy = energy - c.IONIZATION
                                    file.write("{}\n".format(kinetic_energy))

                                else:
                                    stock = -3
                                    file.write("{}\n".format(stock))
                            else:
                                stock = -2
                                file.write("{}\n".format(stock))
                        else:
                            stock = -1
                            file.write("{}\n".format(stock))
                    else:
                        stock = -1
                        file.write("{}\n".format(stock))
    return S, p
