"""
Based on AtomLoadingOptimizer.py, but reduced to simply scan one of the coils.
This is intended to find the distribution of atoms loaded vs a given coil voltage

Intended to be run after first (manually) optimizing the coil values for atom loading
"""

from artiq.experiment import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from utilities.BaseExperiment import BaseExperiment


class AtomLoading1DCoilScan(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # overwrite the experiment variables of the same names
        # self.setattr_argument("t_SPCM_exposure", NumberValue(10 * ms, unit='ms'))
        self.setattr_argument("atom_counts_threshold", NumberValue(180))
        self.setattr_argument("t_MOT_loading", NumberValue(500 * ms, unit='ms'))
        self.setattr_argument("n_measurements", NumberValue(400, type='int', scale=1, ndecimals=0, step=1)) # may be necessary to up this if loading rate < 1%
        group = "coil volts"
        self.setattr_argument("dV_sequence", StringValue("np.linspace(-0.1,0.1,20)"), group)


        # this should be close to the mean signal from the atom
        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        self.sampler_buffer = np.zeros(8)
        self.dV_arr = eval(self.dV_sequence)
        self.AZ_bottom_arr = self.AZ_bottom_volts_MOT + self.dV_arr
        self.AZ_top_arr = self.AZ_top_volts_MOT + self.dV_arr
        self.AX_arr = self.AX_volts_MOT + self.dV_arr
        self.AY_arr = self.AY_volts_MOT + self.dV_arr

        # self.minimize_args = eval(self.scipy_minimize_args)
        # self.minimize_options = eval(self.scipy_minimize_options)

        self.coil_volts_default = np.array([self.AZ_bottom_volts_MOT,
                                   self.AZ_top_volts_MOT,
                                   self.AX_volts_MOT,
                                   self.AY_volts_MOT])

        # this returns an array of all rows equal to zeros
        # self.coil_volts_steps = np.zeros((4*len(self.dV_arr),4))
        # coil_volts = np.zeros(4)
        # for coil_i in range(4):
        #     for j,dV in enumerate(self.dV_arr):
        #         coil_volts = self.coil_volts_default
        #         coil_volts[coil_i] += dV
        #         print(coil_volts)
        #         for col in range(4):
        #             self.coil_volts_steps[j,col]=coil_volts[col]

        # print(self.coil_volts_steps)

        # todo: assert that the coil_volts_defaults are within the boundaries

        self.counts_list = np.zeros(self.n_measurements)
        self.set_dataset(self.count_rate_dataset,
                         [0.0],
                         broadcast=True)

        self.cost_dataset = "cost"
        self.set_dataset(self.cost_dataset,
                         [0.0],
                         broadcast=True)


    def run(self):
        self.initialize_hardware()
        self.warm_up()

        # for coil_i in range(4):
        #
        #     for dV in self.dV_arr:

                # this is stupid
                # for coil_j in range(4):
                #     if coil_j == coil_i:
                #         self.coil_volts[coil_j] = self.coil_volts_default[coil_j] + dV
                #     else:
                #         self.coil_volts[coil_j] = self.coil_volts_default[coil_j]

                # # artiq gets this wrong
                # self.coil_volts = self.coil_volts_default
                # self.coil_volts[coil_i] += dV

                # self.optimization_routine(self.coil_volts)

        # coil_volts = np.zeros(4)
        # for coil_i in range(4):
        #     for j in range(len(self.dV_arr)):
        #         coil_volts = self.coil_volts_default
        #         coil_volts[coil_i] += self.dV_arr[j]
        #         self.optimization_routine(coil_volts)

        # for coil_volts in self.coil_volts_steps:
        #
        #     self.optimization_routine(coil_volts)

        """
        none of the above correctly set the coils (ARTIQ fucks up the addition), so just do this the dumbest way
        we can conceive of:
        """
        for i in range(len(self.dV_arr)):
            self.optimization_routine(np.array([self.AZ_bottom_arr[i],
                                               self.AZ_top_volts_MOT,
                                               self.AX_volts_MOT,
                                               self.AY_volts_MOT]))

        for i in range(len(self.dV_arr)):
            self.optimization_routine(np.array([self.AZ_bottom_volts_MOT,
                                                self.AZ_top_arr[i],
                                                self.AX_volts_MOT,
                                                self.AY_volts_MOT]))

        for i in range(len(self.dV_arr)):
            self.optimization_routine(np.array([self.AZ_bottom_volts_MOT,
                                                self.AZ_top_volts_MOT,
                                                self.AX_arr[i],
                                                self.AY_volts_MOT]))

        for i in range(len(self.dV_arr)):
            self.optimization_routine(np.array([self.AZ_bottom_volts_MOT,
                                                self.AZ_top_volts_MOT,
                                                self.AX_volts_MOT,
                                                self.AY_arr[i]]))



    @kernel
    def initialize_hardware(self):
        self.base.initialize_hardware()

    @kernel
    def warm_up(self):
        """hardware init and turn things on"""

        self.core.reset()

        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)

        self.dds_cooling_DP.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_FORT.sw.on()

        # delay for AOMs to thermalize
        delay(2 * s)

        self.core.break_realtime()
        self.laser_stabilizer.run()
        delay(1 * ms)
        self.dds_FORT.sw.on()
        delay(1*s)

    def get_cost(self, data: TArray(TFloat,1)) -> TInt32:
        q = 0 # whether the counts exceeded the atom threshold
        atoms_loaded = 0
        q_last = (data[0] > self.atom_counts_threshold)
        for x in data[1:]:
            q = x > self.atom_counts_threshold
            if q != q_last and q_last:
                atoms_loaded += 1
                self.print_async("optimizer found the single atom signal!")
            q_last = q
        # return -1 if 0 < atoms_loaded < 10 else -1*atoms_loaded # threshold should depend on number of measurements
        return -1 * atoms_loaded

    @kernel
    def optimization_routine(self, coil_values: TArray(TFloat)):# -> TInt32:
        """
        the function that will be passed to the optimizer

        coil_values: array of coil control volt values
        return:
            cost: the cost for the optimizer
        """
        self.core.reset()
        delay(1*ms)

        self.zotino0.set_dac(coil_values, channels=self.coil_channels)
        self.print_async(coil_values)
        delay(1*ms)

        self.laser_stabilizer.run()
        delay(1*ms)
        self.dds_FORT.sw.on()

        delay(self.t_MOT_loading)

        # reset the counts dataset each run so we don't overwhelm the dashboard when plotting
        self.set_dataset(self.count_rate_dataset,[0.0],broadcast=True)

        for i in range(self.n_measurements):
            t_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
            counts_per_s = self.ttl0.count(t_end)/self.t_SPCM_exposure
            delay(1 * ms)
            self.append_to_dataset(self.count_rate_dataset, counts_per_s)
            self.counts_list[i] = counts_per_s*self.t_SPCM_exposure

        cost = self.get_cost(self.counts_list)
        self.append_to_dataset(self.cost_dataset, cost)

        # return cost

    def analyze(self):
        """plot the results here"""

        fig,ax = plt.subplots(nrows=4)
        labels = ['AZ_top','AZ_top','AX','AY']
        colors = ['yellow','cyan','magenta','blue']
        for i,ax in enumerate(ax.flat):
            ax.plot(self.cost_dataset[i],c=colors[i],label=labels[i])





