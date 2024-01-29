from artiq.experiment import *
from numpy.random import normal
from collections import namedtuple

#from utilities.dropsimulation import *
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
        self.setattr_argument("temp", NumberValue(40.0, step = 0.25))

    def run(self):
        w0 = 2.5e-6  # [m]
        TFORT = 1.5e-3  # [K]
        Tatom = 5e-5
        steps = 100
        lmda =1.064e-6
        tlist = linspace(0, steps, steps)  # time [us]


        reten_1 = retention_at_t(tlist, Tatom, base_retention=.95)


        self.set_dataset("tlist", tlist, broadcast=True)
        self.set_dataset("ret_1", reten_1, broadcast=True)


        if p0 == None:
            p0 = [4.e-5]



        TLIST = linspace(0, 160, 20)  # time [us]
        R1 = retention_at_t(TLIST, 5e-5, base_retention=0.90)
        args = TLIST, R1, p0

        Topt, ropt, modeled_r = self.analyze(args)

        # 72 simulates count above threshold
        sigma = np.ones((len(modeled_r))) / sqrt(72)

        self.set_dataset("optimal_temp_k", Topt, broadcast=True)
        self.set_dataset("r_opt", ropt, broadcast=True)
        self.set_dataset("modeled_r", modeled_r, broadcast=True)
        self.set_dataset("sigma", sigma, broadcast=True)
        self.set_dataset("model_t_data", TLIST, broadcast=True)
        self.set_dataset("model_ret_data", R1, broadcast=True)



    def create_applet(self):
        return -1

    def analyze(self, args):

        return start_modeling("temperature", args)
