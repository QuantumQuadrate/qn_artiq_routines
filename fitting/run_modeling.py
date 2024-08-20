"""

"""

from numpy import *
import numpy as np
from numpy.random import normal
from scipy.optimize import curve_fit
from random import random as rand
from matplotlib import pyplot as plt
import time
#### local files
from utilities.physics.rbconsts import *
from utilities.physics.rbensemble import RbEnsemble as ensemble


def release_recap_retention_at_t(t, T, base_retention, Tdepth=1e-3, wx=0.7e-6, wy=None, lmda=8.52e-7, events=1000):
        """ Procedure for simulating a release ("drop") and recapture experiment
        to deduce the temperature of actual atoms in such an experiment.

        Based on code by Mark, with some corrections
        't': time at which to evaluate the simulation (us). can be scalar or 1D numpy array
        'wx': waist
        'Tdepth': FORT temperature depth (K)
        'T': atom temp (K)
        'tmax': max time in units us
        'steps': number of FORT drop outs
        'events': number of release-recapture events per photocounts pt
        'wy': optional waist for elliptical FORT
        'lmda': optional wavelength of the FORT

        Note: the time entered in microseconds empirically gives a much more reliable fit value
        compared to entering the values in seconds, likely due to precision loss somewhere.
        """

        if wy is None:
                wy = wx

        # print(f"scipy is trying T={T*1e6:.2f}K, base_retention={base_retention:.2f}...")
        # print(f"scipy is trying wx = {wx}, wy = {wy}, lmda = {lmda}, Tdepth = {Tdepth}, ")
        umax = kB * Tdepth
        zR = pi * wx ** 2 / lmda

        omega_r = (1 / sqrt((wx ** 2 + wy ** 2) / 2)) * sqrt(2 * kB * Tdepth / mRb)  # radial trap frequency

        omega_z = (1 / zR) * sqrt(2 * kB * Tdepth / mRb)  # axial trap frequency

        # define potential energy function
        def U(x, y, z):
                ww = (1 + z ** 2 / zR ** 2)
                return -umax * exp(-2 * x ** 2 / (wx ** 2 * ww) - 2 * y ** 2 / (wy ** 2 * ww)) / ww

        dx = dy = sqrt(kB * T / (mRb * omega_r ** 2))
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

        escape = float64(empty(len(t)))

        for i in range(events):
                KE = .5 * mRb * ((vxlist[i] - g * t) ** 2 + vylist[i] ** 2
                                 + vzlist[i] ** 2)
                PE0 = U(xlist[i], ylist[i], zlist[i])
                PE = U(xlist[i] + t * vxlist[i] + .5 * g * (t) ** 2,
                       ylist[i] + t * vylist[i],
                       zlist[i] + t * vzlist[i])

                # print(f"t = {t}")

                init_E = KE + PE0
                hot_E = KE + PE

                escape += ((float64((greater(hot_E, 0)))) - (float64((greater(init_E, 0)))))

        retention = base_retention * (1 - escape / float64(events))

        # print(f"finished. T={T * 1e6} [uK], r = {base_retention}")

        return retention

def get_release_recap_fit_result(tlist, retention, p0=None, bounds=None, retention_at_t_kwargs={}):
        """

        :param tlist:
        :param retention:
        :param p0:
        :param bounds: array of 2-tuple of boundaries for temp in K and retention with 1st (2nd) tuple minima (maxima)
        :param retention_at_t_kwargs: keyword arguments for release_recap_retention_at_t which is used as the model
                for the fit. these arguments might specify, e.g., the trap parameters
        :return: tuple of popt from curve_fit, y points from model evaluated with fit params
        """
        if p0 is None:
                p0 = [4e-5, retention[0]]
        if bounds is None:
                lower_bounds = np.array([1e-7, 0.0])
                upper_bounds = np.array([5e-4, 1.0])
                bounds = (lower_bounds, upper_bounds)

        model = lambda t, T, r: release_recap_retention_at_t(t, T, r, **retention_at_t_kwargs)

        popt, pcov = curve_fit(model, tlist, retention, p0=p0, absolute_sigma=False, bounds=bounds,
                               x_scale=(1e-6,1), diff_step=0.1)

        modeled_y = release_recap_retention_at_t(tlist, *popt, **retention_at_t_kwargs)
        return popt, modeled_y


def atom_loading_fit(xdata=None, p0=None, bin_count=40, measurements=500):
        ret_args = {}
        print("Proceding with atom_loading_fit at time:", time.time())
        factor = 1 / 100
        if False:
                print("Incorrect parameters passed: exiting with return -1")
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
                print(f"Otsu Threshold is at {otsu_threshold} counts")
                domain = [0, 300]
                counts_pruned = np.array([x for x in xdata if x < domain[1]])
                ypts, bins, _ = plt.hist(counts_pruned * factor, bins=bin_count)
                plt.show()
                xpts = np.linspace(min(counts_pruned) * factor, max(counts_pruned) * factor, len(bins) - 1)
                atoms_loaded = np.sum(counts_pruned[:] >= otsu_threshold)
                print("Fitting to integrated counts")
                if p0!=None:
                        p0[2] = p0[2] * factor
                        p0[3] = p0[3] * factor
                        p0[4] = p0[4] * factor
                        p0[5] = p0[5] * factor
                        popt, pcov = curve_fit(double_poissonian_model, xpts, ypts, p0=p0, ftol=1e-15, maxfev=10000000,
                                               xtol=1e-15,)
                else:
                        popt, pcov = curve_fit(double_poissonian_model, xpts, ypts, ftol=1e-15, maxfev=10000000,
                                               xtol=1e-15)


                # rint(popt)
                popt[2] = popt[2] * (1/factor)
                popt[3] = popt[3] * (1/factor)
                popt[4] = popt[4] * (1/factor)
                popt[5] = popt[5] * (1/factor)
                y_dat = double_poissonian_model(xpts * (1/factor), *popt)
                x_dat = xpts * (1/factor)
                ret_args['atoms_loaded'] = atoms_loaded
                ret_args['otsu_threshold'] = otsu_threshold
                ret_args['opt_params'] = popt
                ret_args['modeled_func'] = x_dat, y_dat

                return ret_args


def start_modeling(model = "temperature", args=None):
        starting_time = time.time()
        print(f"Attempting to run: {model} at {starting_time}")

        if model == "temperature":
                """
                *args = (xdata, retention, p0,)
                xdata should be time in terms of micro seconds
                if p0 is not provided, p0 = (40(uK), retention[0], Tdepth = 1e-3, wx = 0.7e-6, wy =0.7e-6 )
                """
                ret_args = get_release_recap_fit_result(*args)
                print(f"Completed: {model} after {(time.time() - starting_time)} seconds")
                return ret_args

        elif model == "count_dist":
                """ 
                *args = (xdata, p0, bin_count,)
                p0 should be provided with good guess to ensure consistency, will be scaled down but rescaled
                on return of optimal parameters
                bin_count = number of points the count photocounts will be summed to 
                """
                ret_args = atom_loading_fit(*args)
                print(f"Completed: {model} after {time.time()-starting_time} seconds")
                return ret_args
        else:
                print(f"Model :{model} was not found after {time.time()-starting_time} seconds")
                return -1