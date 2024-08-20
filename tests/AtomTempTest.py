"""
An experiment to test Jake Uribe's release recapture temperature fitting code

Runs a Monte-Carlo simulation to generate simulated data for a given atom temp and hardcoded trap params, then
fits the simulated data to the Monte Carlo simulation to extract a temperature and baseline retention.
"""

from artiq.experiment import *
from numpy import *
from collections import namedtuple

import sys, os
# get the current working directory
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from fitting.run_modeling import *


class AtomExperimentSim(EnvExperiment):
    """
    Atom Temp Test
    """

    def build(self):
        Variable = namedtuple("Variable", "name value value_type kwargs")

        self.setattr_argument("measurements", NumberValue(42, step=1))
        self.setattr_argument("bins", NumberValue(42, step=1))
        self.setattr_argument("bin width", NumberValue(42, step=1))
        self.setattr_argument("atom_temp_uK", NumberValue(50.0, step=0.25))
        self.setattr_argument("atom_temp_guess_uK", NumberValue(50.0, step=0.25))
        self.setattr_argument("atom_retention_guess", NumberValue(0.95, step=0.25))


    def run(self):
        w0 = 2.5e-6  # [m]
        TFORT = 1.5e-3  # [K]
        FORT_lmbda = 8.52e-7 # [m]
        tlist = array([1, 2, 5, 10, 15, 20, 50, 75, 100, 150])  # [us]

        trap_params = {'Tdepth': TFORT, 'wx': w0, 'lmda': FORT_lmbda}

        model = lambda t, T, r: release_recap_retention_at_t(t, T, r, **trap_params)
        simulated_retention = model(tlist, self.atom_temp_uK*1e-6 , 0.95)

        self.set_dataset("tlist", tlist, broadcast=True)

        p0 = [self.atom_temp_guess_uK*1e-6, self.atom_retention_guess]

        trap_params = {'Tdepth': TFORT, 'wx': w0, 'lmda': FORT_lmbda}
        popt, fit_retention = get_release_recap_fit_result(tlist, simulated_retention, p0, retention_at_t_kwargs=trap_params)

        print(f"finished. fit: T={popt[0] * 1e6} [uK], r = {popt[1]}")

        self.set_dataset("optimal_temp_k", popt[0], broadcast=True)
        self.set_dataset("r_opt", popt[1], broadcast=True)
        self.set_dataset("fit_ret", fit_retention, broadcast=True)

        # self.set_dataset("sigma", sigma, broadcast=True)
        self.set_dataset("model_t_data", tlist, broadcast=True)
        self.set_dataset("model_ret_data", simulated_retention, broadcast=True)

    def create_applet(self):
        return -1

    # def model(self, args = None):
    #
    #     return start_modeling("temperature", args)
