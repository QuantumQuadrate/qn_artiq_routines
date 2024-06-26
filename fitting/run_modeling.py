#### libraries
# np.seterr(all='raise') # raise errors, don't just print warnings
from numpy import *
from numpy.random import normal
from scipy.optimize import curve_fit
from random import random as rand
from matplotlib import pyplot as plt
import time
#### local files
from utilities.physics.rbconsts import *
from utilities.physics.rbensemble import RbEnsemble as ensemble


def retention_at_t(t, T=None, base_retention=0.9, Tdepth = 1e-3, wx = 0.7e-6, wy =0.7e-6, lmda = 8.52e-7):
        """ Procedure for simulating a release ("drop") and recapture experiment
            to deduce the temperature of actual atoms in such an experiment.

            Based on code by Mark, with some corrections
            'wx': waist
            'Tdepth': FORT temperature depth
            'T': atom temp
            'tmax': max time in units us
            'steps': number of FORT drop outs
            'events': number of release-recapture events per photocounts pt
            'wy': optional waist for elliptical FORT
            'lmda': optional wavelength of the FORT
        """
        T = abs(T)
        events = 20000

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
                base_retention = 1  #the retention baseline with no fort drop

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

        print(f"finished. T={T * 1e6} [uK], r = {base_retention}")

        return retention

def temp(tlist, retention, p0 = None):
        if p0 == None:
                p0 = [4e-5, retention[0], 1e-3, 0.7e-6, 0.7e-6]


        popt, pcov = curve_fit(retention_at_t, tlist, retention, p0=p0, absolute_sigma=False,
                                                   maxfev=1000000000, epsfcn=1)

        modeled_y = retention_at_t(tlist, *popt)
        return popt, modeled_y


def atom_loading_fit(xdata=None, p0 = None, bin_count = 40, measurements = 500):
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


def start_modeling(model = "temperature", args = None):
        starting_time = time.time()
        print(f"Attempting to run: {model} at {starting_time}")

        if model == "temperature":
                """
                *args = (xdata, retention, p0,)
                xdata should be time in terms of micro seconds
                if p0 is not provided, p0 = (40(uK), retention[0], Tdepth = 1e-3, wx = 0.7e-6, wy =0.7e-6 )
                """
                ret_args = temp(*args)
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