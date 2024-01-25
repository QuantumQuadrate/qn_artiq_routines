from artiq.experiment import *
import numpy as np
from numpy import linspace,empty,sin,cos,log,exp,sqrt,meshgrid,array
from numpy.random import normal
from scipy.optimize import curve_fit
import numpy as np
from collections import namedtuple
from math import pi,e
from random import random as rand
import time

from utilities.physconsts import *

from utilities.rbconsts import *
from utilities.dipole_trap import dipole_trap
from utilities.rbensemble import RbEnsemble as ensemble
#from utilities.dropsimulation import *
from repository.run_modeling import *
class AtomExperimentSim(EnvExperiment):
    """
    Atom Temp Test
    """

    def build(self):
        Variable = namedtuple("Variable", "name value value_type kwargs")

        self.setattr_argument("measurements", NumberValue(42, step=1))
        self.setattr_argument("bins", NumberValue(42, step=1))
        self.setattr_argument("bin width", NumberValue(42, step=1))


    def run(self):
        w0 = 2.5e-6  # [m]
        TFORT = 1.5e-3  # [K]
        Tatom = 5e-5  # [K]
        steps = 100
        lmda =1.064e-6
        tlist = linspace(0, steps, steps)  # time [us]

        tempexp = dipole_trap(w0, lmda, TFORT, Tatom)
        # tempexp.distplot(100) # example of atoms in xz plane of FORT at t=0
        tlist_1, reten_1 = tempexp.drop_recap(tlist, base_retention=.95)
        x, y = tempexp.distplot(100)

        self.set_dataset("x_trap", x, broadcast=True)
        self.set_dataset("y_trap", y, broadcast=True)
        self.set_dataset("tlist_1", tlist_1, broadcast=True)
        self.set_dataset("ret_1", reten_1, broadcast=True)
        # try to fit the datas
        real_t_data = array([0.00000000000000000000000000, 2.0, 4.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0,
                        20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0])
        real_ret_data = [0.94079207920792083, 0.94999999999999996, 0.95000000000000005,
                    0.89000000000000001, 0.88999999999999995, 0.90000000000000002,
                    0.93999999999999995, 0.89999999999999998, 0.85999999999999999,
                    0.81000000000000005, 0.74000000000000004, 0.50000000000000001,
                    0.28712871287128711, 0.27000000000000002, 0.16, 0.10000000000000001,
                    0.089999999999999997]

        """w0 = 2.5e-6  # [m]
        TFORT = 1.5e-3  # [K]
        Tatom = 3.e-5  # [K] # this will be a fit guess
        steps = 100
        tempexp = dipole_trap(w0,lmda, TFORT, Tatom)

        tlist, reten = tempexp.drop_recap(tlist, base_retention=.95)"""

        def drop_recap(tlist, T=None, events=None, base_retention=None, Tdepth=None,
                       progress=False, b=1, c=0):
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
            zR = pi * wx ** 2 / lmda

            tlist = 1e-6 * tlist

            if events is None:
                events = 2000
            if base_retention is None:
                base_retention = 1  # the retention baseline with no fort drop

            retention = empty(len(tlist))

            xlist, ylist, zlist = xdist(T, events)
            vzlist, vxlist, vylist = vdist(T, events)
            # print(xlist)
            g_scaled = g * 1e3
            for j, t in enumerate(tlist):

                escape = 0
                nhot = 0  # this is an untrapped atom

                for i in range(events):
                    hot = 0
                    KE = .5 * mRb * ((vxlist[i] - g * t) ** 2 + vylist[i] ** 2
                                     + vzlist[i] ** 2)
                    PE0 = U(xlist[i], ylist[i], zlist[i])
                    PE = U(xlist[i] + t * vxlist[i] + .5 * g * (t) ** 2,
                           ylist[i] + t * vylist[i],
                           zlist[i] + t * vzlist[i])

                    if KE + PE0 > 0:
                        hot = 1
                    nhot += hot
                    if KE + PE > 0:
                        escape += 1 - hot
                retention[j] = base_retention * (1 - escape / events)

                if progress is not False:
                    if j % 10 == 0:
                        print(f"timestep {j}: t = {t * 1e6:.0f} [us], ret = {retention[j]:.2f}")

            retention = (retention * b - c)
            tlist *= 1e6
            plt.plot(tlist, retention, label=f'{T / 1e-6:.0f} uK')  # show time in units [us]
            #         plt.xlabel("time [us]")
            #         plt.ylabel("retention")
            #         plt.ylim((0,1))
            #         plt.legend()
            #         plt.show()

            print(f"finished. T={T * 1e6} [uK], r = {base_retention}")
            return tlist, retention





        p0 = [4.e-5]
        print("p0::::",type(p0))
        upper_bounds = [9e-4]
        lower_bounds = [1e-5]



        TLIST = linspace(0, 160, 20)  # time [us]
        R1 =retention_at_t_3(TLIST, 5e-5, base_retention=0.9)
        args = {}
        args = TLIST, R1
        Topt, ropt, modeled_r = start_modeling("temperature", args)

        sigma = np.ones((len(modeled_r))) / sqrt(72)

        self.set_dataset("optimal_temp_k", Topt, broadcast=True)
        self.set_dataset("r_opt", ropt, broadcast=True)
        self.set_dataset("modeled_r", modeled_r, broadcast=True)
        self.set_dataset("sigma", sigma, broadcast=True)
        self.set_dataset("model_t_data", TLIST, broadcast=True)
        self.set_dataset("model_ret_data", R1, broadcast=True)



    def create_applet(self):
        return -1

    def analyze(self):

        return -1

    def ct_rate_simulation(self):
        s = 14 * 6  # they say s ~ 9 but later
        s = 0.8 * 6  # 6 beams, assumed non-interfering
        gamma = 2 * np.pi * 6.06e6
        det = 2.5 * gamma
        r = gamma * 0.5 * s / (1 + 4 * (
                    det / gamma) ** 2 + s)  # single beam scattering rate, two level atom. Mark's notes, Eq. 5.21
        # NA = 0.43
        # eta_lens = np.arcsin(NA)/(2*np.pi)
        eta_lens = 0.02
        eta_fiber = 0.6  # from their simulation of 780 mode overlap with the fiber
        eta_det = 0.58  # APD
        eta_AOM = 0.7
        eta_bp = 0.9 * 0.94 * 0.95 * 0.96  # filter, bandpass, dichroic, MM fiber
        other = 0.9  # spliced fiber after damage

        # expected
        eta_total = eta_lens * eta_fiber * eta_det * eta_AOM * eta_bp * other
        print("expected collection efficiency", eta_total)
        cnt_rate = r * eta_total
        print(f"{cnt_rate:.2e} counts/s:")

        # measured
        eta_total = 0.0031
        print("measured collection efficiency", eta_total)
        cnt_rate = r * eta_total
        print(f"{cnt_rate:.2e} counts/s:")

        cnt_rate*0.01 # counts during exposure

        # they estimate 2.3% of photons collected by the lens for isotropic
        # emission. how?
        NA = 0.43
        eta_lens = np.arcsin(NA) / (2 * np.pi)
        eta_lens

        """
        Networking Experiment
        """
        s = 11 * 6  # 6 beams, assumed non-interfering
        ac_shift = 0
        det = 2 * np.pi * (15.5e6 + ac_shift)
        gamma = 2 * np.pi * 6.06e6
        r = gamma * 0.5 * s / (1 + 4 * (
                    det / gamma) ** 2 + s)  # single beam scattering rate, two level atom. Mark's notes, Eq. 5.21
        NA = 0.61
        eta_lens = np.arcsin(NA) / (2 * np.pi)
        eta_fiber = 0.7
        eta_det = 0.7
        eta_bp = 0.9  # bandpasses
        print(f"{r * eta_lens * eta_fiber * eta_bp * eta_det:.2e} counts/s (trap off):")
        ac_shift = 20e6
        det = 2 * np.pi * (15.5e6 + ac_shift)
        r = gamma * 0.5 * s / (1 + 4 * (
                    det / gamma) ** 2 + s)  # single beam scattering rate, two level atom. Mark's notes, Eq. 5.21
        print(f"{r * eta_lens * eta_fiber * eta_bp * eta_det:.2e} counts/s (trap on):")

        # estimate s for the MOT beams
        eta_7030 = 0.7
        eta_AOM = 0.7 * 0.7  # the AOMs are operated at roughly 70% of their max eff., which is about 70%
        MOTchip_mW = 1  # before the beamsplitters and fiber launchers for MOT1-4
        wMOT_cm = 0.06  # 600 um waist
        P_MOT_mW = (eta_AOM * eta_7030 * eta_fiber * MOTchip_mW)  # single beam power (mW/cm^2)
        I_MOT = 2 * P_MOT_mW / (np.pi * wMOT_cm ** 2)
        print("saturation parameter for one MOT beam")
        print(I_MOT / 3.576)
        print(P_MOT_mW)

        return -1

    def simulated_distribution(self):
        mu_bg = 100
        mu_sig = 10 # the true signal mean
        poisson = lambda x, m: ((m ** x) / np.math.factorial(x)) * np.exp(-m)
        func = lambda counts: poisson(counts, mu_bg) + poisson(counts, mu_sig+mu_bg)
        counts = range(1,600)
        plot_data = {"x": counts, "y": np.array([func(x) for x in counts])}
        return plot_data

    def simulated_sampling(self):
        mu_bg = 100
        mu_sig = 100  # the true signal mean
        func = lambda counts: poisson(counts, mu_bg) + poisson(counts, mu_sig + mu_bg)
        poisson = lambda x, m: ((m ** x) / np.math.factorial(x)) * np.exp(-m)
        n = 1000  # number of measurements
        domain = [0, 500]
        x1, x2 = domain

        fmax = poisson(mu_bg, mu_bg)  # the maximum
        y_dist = np.empty(n)
        f_dist = np.empty(n)
        x_dist = np.empty(n)  # this is the distribution we want
        j = 0  # dist index
        while j < n:
            x = int((x2 - x1) * np.random.rand() + 0.5)  # rand val on domain of f(x)
            f = func(x)
            y = np.random.rand() * fmax  # rand val on range of f(x)
            if y <= f:
                y_dist[j] = y
                f_dist[j] = f
                x_dist[j] = x  # x vals with approximate gaussian pdf
                j += 1

        function_curve = {"x": x_dist, "y": f_dist}
        ct_curve = {"x": x_dist, "y": y_dist}
        return {"expected": function_curve, "measured": ct_curve}
