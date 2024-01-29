#### libraries
import numpy as np
# np.seterr(all='raise') # raise errors, don't just print warnings
from numpy import *
from numpy.linalg import eig
from numpy.random import normal
from scipy.optimize import curve_fit
from random import random as rand
from scipy import stats
from scipy.stats import poisson
from matplotlib import pyplot as plt
import time
#### local files
from utilities.physconsts import *
from utilities.rbconsts import *
from utilities.rbensemble import RbEnsemble as ensemble
from utilities.dipole_trap import dipole_trap


def retention_at_t(t, T=None, base_retention=0.9):
        """ Procedure for simulating a release ("drop") and recapture experiment
            to deduce the temperature of actual atoms in such an experiment.

            Based on code by Mark, with some corrections
            'wx': waist
            'Tdepth': FORT temperature depth
            'T': atom temp
            'tmax': max time in units us
            'steps': number of FORT drop outs
            'events': number of release-recapture events per data pt
            'wy': optional waist for eliptical FORT
        """
        T = abs(T)
        events = 20000
        factor = 1
        w0 = 2.5e-6 * factor  # [m]
        Tdepth = 1.5e-3 * factor  # [K]
        Tatom = 5e-5  # * factor# [K]
        steps = 10
        lmda = 1.064e-6 * factor
        wx = w0
        wy = wx
        TFORT = Tdepth
        umax = kB * Tdepth
        zR = pi * wx ** 2 / lmda

        omega_r = (1 / sqrt((wx ** 2 + wy ** 2) / 2)) * sqrt(2 * kB * TFORT / mRb)  # radial trap frequency

        omega_z = (1 / zR) * sqrt(2 * kB * TFORT / mRb)  # axial trap frequency

        # Tdepth = Tdepth
        def U(x, y, z):
                ww = (1 + z ** 2 / zR ** 2)
                return -umax * exp(-2 * x ** 2 / (wx ** 2 * ww) - 2 * y ** 2 / (wy ** 2 * ww)) / ww

        dx = dy = sqrt(kB * T / (mRb * omega_r ** 2))
        dx
        dz = sqrt(kB * T / (mRb * omega_z ** 2))
        zlist = normal(0, dz, size=events)
        xlist = normal(0, dx, size=events)
        ylist = normal(0, dy, size=events)

        atoms = ensemble(T)

        vlist = atoms.sampling_maxboltzv(events, [0, 2])  # speeds

        vxlist = empty(events)
        vylist = empty(events)
        vzlist = empty(events)

        for i in range(events):
                ex = 2 * rand() - 1
                ey = 2 * rand() - 1
                ez = 2 * rand() - 1
                v = vlist[i]
                A = sqrt(ex ** 2 + ey ** 2 + ez ** 2)
                vxlist[i] = ex * v / A
                vylist[i] = ey * v / A
                vzlist[i] = ez * v / A
                # print(f"A:{i}::",A)

        t = t * 1e-6

        if events is None:
                events = 20000
        if base_retention is None:
                base_retention = 1  # the retention baseline with no fort drop

        # xlist, ylist, zlist = xdist(T,events)
        # vzlist, vxlist, vylist = vdist(T,events)
        # print(xlist)
        escape = np.float64(np.empty(len(t)))
        nhot = 0  # this is an untrapped atom

        for i in range(events):
                hot = 0
                KE = .5 * mRb * ((vxlist[i] - g * t) ** 2 + vylist[i] ** 2
                                 + vzlist[i] ** 2)
                PE0 = U(xlist[i], ylist[i], zlist[i])
                PE = U(xlist[i] + t * vxlist[i] + .5 * g * (t) ** 2,
                       ylist[i] + t * vylist[i],
                       zlist[i] + t * vzlist[i])

                # print(f"t = {t}")

                init_E = KE + PE0
                hot_E = KE + PE

                escape += (np.float64((np.greater(hot_E, 0))))

        retention = base_retention * (1 - escape / np.float64(events))

        print(f"finished. T={T * 1e6} [uK], r = {base_retention}, b = , c = ")

        return retention

def temp(tlist, retention, p0 = None):
        if p0 == None:
                p0 = [4e-5,retention[0]]

        print(tlist, retention)
        popt, pcov = curve_fit(retention_at_t, tlist, retention, p0=p0, absolute_sigma=False,
                                                   maxfev=1000000000, epsfcn=1)
        Topt, ropt = popt
        modeled_y = retention_at_t(tlist,Topt,ropt)
        return Topt, ropt, modeled_y


def atom_loading_fit(xdata=None, p0 = None):
        ret_args = {}
        print("Proceding with atom_loading_fit at time:", time.time())

        if False:
                print("Incorrect paratemeters passed: exiting with return -1")
                return -1
        else:
                iteration = 0

                def double_poissonian_model(counts, a=100, b=100, ma=1, mb=1, wa=1, wb=1):
                        """
                                histogram fit. neglect two atom loading
                                a,b: amplitudes of the distribution
                                ma,b: the x offsets
                                wa,b: the dist widths
                        """
                        #print("Proceding with poisson: x")

                        return a * exp(-((counts - ma) / wa) ** 2) + b * exp(-((counts - mb) / wb) ** 2)
                print("Entered fit:")
                def otsu_intraclass_variance(image, threshold):
                        """
                                Otsu’s intra-class variance.
                                If all pixels are above or below the threshold, this will throw a warning that can safely be ignored.
                        """
                        return np.nansum([
                                np.mean(cls) * np.var(image, where=cls)
                                #   weight   ·  intra-class variance
                                for cls in [image >= threshold, image < threshold]
                        ])
                otsu_threshold = min(range(int(np.min(xdata)) + int(10), int(np.max(xdata))),
                                     key=lambda th: otsu_intraclass_variance(xdata, th)
                                     )
                print(f"Otsu Threshold is {otsu_threshold}")
                domain = [0, 300]
                counts_pruned = np.array([x for x in xdata if x < domain[1]])
                ypts, bins, _ = plt.hist(counts_pruned * 0.01, bins=50)
                plt.close()
                xpts = np.linspace(min(counts_pruned) * 0.01, max(counts_pruned) * 0.01, len(bins) - 1)
                atoms_loaded = np.sum(counts_pruned[:] >= otsu_threshold)
                print("fitting")
                if p0!=None:
                        popt, pcov = curve_fit(double_poissonian_model, xpts, ypts, p0=p0, ftol=1e-15, maxfev=1000000,
                                               xtol=1e-15,)
                else:
                        popt, pcov = curve_fit(double_poissonian_model, xpts, ypts, ftol=1e-15, maxfev=1000000,
                                               xtol=1e-15)


                # rint(popt)
                popt[2] = popt[2]*100
                popt[3] = popt[3] * 100
                popt[4] = popt[4] * 100
                popt[5] = popt[5] * 100
                y_dat = double_poissonian_model(xpts * 100, *popt)
                x_dat = xpts * 100
                ret_args['atoms_loaded'] = atoms_loaded
                ret_args['otsu_threshold'] = otsu_threshold
                ret_args['opt_params'] = popt
                ret_args['modeled_func'] = x_dat, y_dat

                return ret_args


def start_modeling(model = "temperature", args = None):
        if model == "temperature":
                return temp(*args)

        elif model == "count_dist":
                """ 
                *args = (xdata, p0)
                p0 should be provided with good guess to ensure consistency
                """
                return atom_loading_fit(*args)
        else:
                print(f"Model :{model} was not found")
                return -1