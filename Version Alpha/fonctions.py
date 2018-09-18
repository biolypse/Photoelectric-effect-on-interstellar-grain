################################################################################
################################### IMPORT #####################################
################################################################################
import numpy as np
from random import uniform
from matplotlib import pyplot as plt


################################################################################
################################# VARIABLES ####################################
################################################################################
La = 100e-10



################################################################################
################################# FUNCTIONS ####################################
################################################################################
def load(name):
    return np.loadtxt(name)

def non_uniform_generator_P():
    x = []
    for i in range(100000):
        a = (np.random.exponential(La))
        x.append(a)
        i += 1
    plt.hist(x, bins = 1000)
    plt.show()


def non_uniform_generator_Y():
    x = []
    fx = []
    for i in range(1000):
        a = uniform(3, 16)
        u = uniform(0, 1)
        while u > 0.5 * (1 + np.tanh((a - 8) / 2)):
            a = uniform(3, 16)
            u = uniform(0, 1)
        b = 0.5 * (1 + np.tanh((a - 8) / 2))
        x.append(a)
        fx.append(b)

    plt.plot(x, fx, '.')
    plt.show()

def x_enter(row):
    for i in range(len(row)):
        if row[i] == 1:
            return i
