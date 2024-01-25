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


#### local files
from utilities.physconsts import *
from utilities.rbconsts import *
from utilities.rbensemble import RbEnsemble as ensemble
from utilities.dipole_trap import dipole_trap


def retention_at_t_3(t, T=None, base_retention=0.9):
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


def sample_cts_fit(measurements, mu_bg = 100, mu_sg = 100, ):
        """generator function to sample counts from background + single atom distribution"""

        n = measurements

        mu_bg = mu_bg  # the mean background
        mu_sig = 100  # the signal mean

        def count_dist(x_data, bg, sig):
                return 1.5 * poisson.pmf(x_data, bg) + poisson.pmf(x_data - sig, bg)

        domain = [0, 500]  # assume we don't measure fewer than 0 or more than 500 counts
        x1, x2 = domain

        fmax = 1.5 * poisson.pmf(mu_bg, mu_bg)  # the maximum
        counts_per_bin = 10
        # these arrays exist so we could plot the points (x,f) and (x,y) to verify the distribution sampling
        y_dist = np.empty(n)
        f_dist = np.empty(n)
        x_dist = np.empty(n)  # this is the distribution we want
        j = 0  # dist index

        while j < n:

                x = int((x2 - x1) * np.random.rand() + 0.5)  # rand val on domain of f(x)
                f = count_dist(x, mu_bg, mu_sig)
                y = np.random.rand() * fmax  # rand val on range of f(x)
                if y <= f:
                        y_dist[j] = y
                        f_dist[j] = f
                        x_dist[j] = x  # x vals with approximate gaussian pdf
                        j += 1

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

        otsu_threshold = min(
                range(int(np.min(x_dist)) + int(counts_per_bin), int(np.max(x_dist))),
                key=lambda th: otsu_intraclass_variance(x_dist, th)
        )
        error = []
        load_error = np.array([1 / np.sqrt(n) if n > 0 else 0 for n in y_dist])
        print(len(load_error))
        atoms_loaded = np.sum(x_dist[:] >= otsu_threshold)

def temperature(tlist, retention, p0 = None):
        if p0 == None:
                p0 = [4e-5,retention[0]]

        print(tlist, retention)
        popt, pcov = curve_fit(retention_at_t_3, tlist, retention, p0=p0, absolute_sigma=False,
                                                   maxfev=1000000000, epsfcn=1)
        Topt, ropt = popt
        modeled_y = retention_at_t_3(tlist,Topt,ropt)
        return Topt, ropt, modeled_y

def count_distribution(xdata, ydata, p0 = None):


        return

def start_modeling(model = "temperature", args = None):
        if model == "temperature":
                return temperature(*args)
        elif model == "count_dist":
                return count_distribution(*args)
        else:
                print(f"Model :{model} was not found")
                return -1