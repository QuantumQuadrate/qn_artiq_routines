#### libraries
import numpy as np
# np.seterr(all='raise') # raise errors, don't just print warnings
from numpy import *
from numpy.linalg import eig
from numpy.random import normal
from scipy.optimize import curve_fit
from random import random as rand



#### local files
from utilities.physconsts import *
from utilities.rbconsts import *
from utilities.rbensemble import RbEnsemble as ensemble
from utilities.dipole_trap import dipole_trap

def start_sim(Tatom, base_retention, tlist, units=1):
        """ A dipole trap object with the beams potential and distribution of
            atoms specified by Tatom.
            'wx': x beam waist in focal plane (z=0)
            'wy': Assumed equal to wx by default
            'Tdepth'
            'Tatom'
        """
        w0 = 2.5e-6  # [m]
        TFORT = 1.5e-3  # [K]
        steps = 100
        lmda =1.064e-6
        tempexp = dipole_trap(w0, lmda, TFORT, Tatom)

        tlist, reten = tempexp.drop_recap(tlist = tlist,T = Tatom, base_retention=base_retention, events = 100)



        return reten
