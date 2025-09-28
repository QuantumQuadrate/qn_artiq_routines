import numpy as np
from . import FitBase


def parameter_initialiser(x, y, p):
    p["C"] = 1.0
    p["f0"] = np.mean(x)
    p["sigma"] = 50e3

def fitting_function(x, p):
    y = (p["C"] * np.exp(-0.5 * ((x - p["f0"]) / p["sigma"]) ** 2))

    # y = np.where(x <= p['x0'], p['k1'] * (x - p['x0']) + p['d'],
    #              p['k2'] * (x - p['x0']) + p['d'])
    return y


resonance_positive_dip = FitBase.FitBase(['C', 'f0', 'sigma'],
                             fitting_function,
                             parameter_initialiser=parameter_initialiser)
