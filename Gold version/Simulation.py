#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# license: CC BY-SA
# author: Nico
# version: 1.0
# email: nico.bellemont@gmail.com

import fonctions as f
from matplotlib import pyplot as plt
import numpy as np
import constants as c

################################################################################
# initializing false values to trigger the verifications
choice = "0"
fractale_transit = "initialization"
number_of_grain = "initialization"
sigma_dens = "initialization"
sphere_transit = "initialization"
################################################################################

while choice not in "1234":  # Verify if the choice is one of the four expected
    choice = str(input("Which simulation do you want to realise ? \n1 : \
Several different grains but with the same fractal properties\n2 : A fractal \
grain and a spherical grain\n3 : Several grains with different fractal \
properties\n4 : Several grain of different sizes but with the same fractal \
properties\n"))


if choice == "1":
    while str(number_of_grain) not in "0123456789":  # Verify if the number is an positive interger
        number_of_grain = input("How much grain do you want to compare?\n")
    number_of_grain = int(number_of_grain)  # Cast the string into an integer

    while f.is_int(sigma_dens) is False or float(sigma_dens) < 0:  # Verify if the sigma parameter is a positive float
        sigma_dens = input("Which value do you want to give to the fractal \
parameter (> 0.1)?\n")
    sigma_dens = float(sigma_dens)  # Cast the string into a float

    S, p = f.main_function(number_of_grain, 0, sigma_dens, choice)  # Calling the main function

    for i in range(number_of_grain):  # Loop that create and plot histogram for each grain
        plt.xlabel("Energy (eV)")
        plt.ylabel("Number of photo_electron who escaped from the grain")
        results = np.loadtxt("Results_Files/results_Grain_N100_S{}p{}_B3p0_{}. \
txt".format(S, p, i))  # Loading the file who contains the energy value
        plt.subplot(1, number_of_grain, i+1)  # Craeting subplot in order to have hist side by side
        plt.hist(results, bins=30, range=(0, 13), edgecolor="black",
                 color="red", label="grain n°{}".format(i))  # The range corresponding to the range of photons energies
        plt.legend(loc="upper left", bbox_to_anchor=[0, 1], fancybox=True)
    plt.suptitle("Different grain with the same fractale parameters (1000000\
photons were sent on each grain)")
    plt.show()


elif choice == "2":
    sigma_dens = []  # List that will contains all the values of sigma
    legend = []  # List that will contains all the legend
    number_of_grain = 2

    while f.is_int(fractale_transit) is False or float(fractale_transit) < 0:
        fractale_transit = input("Which value do you want to assign to sigma \
for the fractal grain ?(Recommended : > 0.4)\n")
    sigma_dens.append(float(fractale_transit))  # Adding sigma for a fractal grain to the sigma list

    while f.is_int(sphere_transit) is False or float(sphere_transit) < 0:
        sphere_transit = input("Which value do you want to assign to sigma for \
the spherical grain ?(Recommended : < 0.4)\n")
    sigma_dens.append(float(sphere_transit))  # Adding sigma for a spherical grain to the sigma list

    sigma_dens = sorted(sigma_dens, reverse=True)  # Sorted in descending order

    S, p = f.main_function(number_of_grain, 0, sigma_dens, choice)

    for i in range(number_of_grain):  # Hist are ploted one over the other, thanks to the sorting
        results = np.loadtxt("Results_Files/results_Grain_N100_S{}p{}_B3p0.txt"
                             .format(S[i], p[i]))
        plt.hist(results, bins=30, range=(0, 13), edgecolor="black",
                 label="sigma = {}".format(sigma_dens[i]))
        plt.xlabel("Energy (eV)")
        plt.ylabel("Number of photo_electron who escaped from the grain")
    plt.title("A fractal and a sperical grain (1000000 photons were sent on \
each grain)")
    plt.legend()
    plt.show()


elif choice == "3":
    sigma_dens = []

    while str(number_of_grain) not in "0123456789":
        number_of_grain = input("How much grain do you want to compare?\n")
    number_of_grain = int(number_of_grain)

    for i in range(number_of_grain):
        transit = input("Which value do you want to assign to the fractal \
parameter n°{}?\n".format(i+1))
        while f.is_int(transit) is False or float(transit) < 0:
            transit = input("Which value do you want to assign to the fractal \
parameter n°{}?\n".format(i+1))
        sigma_dens.append(float(transit))

    sigma_dens = sorted(sigma_dens, reverse=True)

    S, p = f.main_function(number_of_grain, 0, sigma_dens, choice)

    for i in range(number_of_grain):
        results = np.loadtxt("Results_Files/results_Grain_N100_S{}p{}_B3p0.txt"
                             .format(S[i], p[i]))
        plt.hist(results, bins=30, range=(0, 13), color=c.COLOR[i],
                 edgecolor="black", label="sigma = {}"
                 .format(sigma_dens[i]))
        plt.xlabel("Energy (eV)")
        plt.ylabel("Number of photo_electron who escaped from the grain")
    plt.title("Different grain with different fractal properties (1000000 \
photons were sent on each grain)")
    plt.legend()
    plt.show()


elif choice == "4":
    GRAIN_RADIUS = []  # List that will contains all the values of the grain's radius
    while str(number_of_grain) not in "0123456789":
        number_of_grain = input("How much grain do you want to compare?\n")
    number_of_grain = int(number_of_grain)

    while f.is_int(sigma_dens) is False or float(sigma_dens) < 0:
        sigma_dens = input("Which value do you want to assign to sigma ?\n")
    sigma_dens = float(sigma_dens)

    for i in range(number_of_grain):
        transit = input("Which value do you want to assign to the grain radius \
n°{}? (between 1e-7 and 1e-8 for acceptable results)\n".format(i+1))
        while f.is_int(transit) is False or float(transit) < 0:
            transit = input("Which value do you want to assign to the grain \
radius n°{}? (between 1e-7 and 1e-8 for acceptable results)\n".format(i+1))

        GRAIN_RADIUS.append(float(transit))

    GRAIN_RADIUS = sorted(GRAIN_RADIUS)

    S, p = f.main_function(number_of_grain, GRAIN_RADIUS, sigma_dens, choice)

    for i in range(number_of_grain):
        results = np.loadtxt("Results_Files/results_Grain_N100_S{}p{}_B3p0_{}.\
txt".format(S, p, i))
        plt.hist(results, bins=30, range=(0, 13), edgecolor="black",
                 label="Radius = {}".format(GRAIN_RADIUS[i]))
        plt.xlabel("Energy (eV)")
        plt.ylabel("Number of photo_electron who escaped from the grain")
    plt.suptitle("Different grain size with the same fractale parameters \
(1000000 photons were sent on each grain)")
    plt.legend()
    plt.show()
