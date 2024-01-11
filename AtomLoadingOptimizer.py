"""
Simple-minded code for moving the MOT to optimize the atom loading rate.

issues with the current method to be addressed:
- there should be some sort of multidimensional scan, rather than a series of 1D scans
- it might be desirable to fine tune the coil values using an optimization algorithm within the 4-dimensional
  volume defined by the boundaries found with the 1D scans, as parking each coil in the middle of the range
   is unlikely to correspond the global optimum.

better approach:
- loop over coil values to find atom signal, but order the scan arrays to alternate direction
 (i.e. [-2,-1,0,1,2] would be ordered [0,-1,1,-2,2])
- once the signal has been found, 

"""

from artiq.experiment import *
import numpy as np
import scipy as sp

from utilities.BaseExperiment import BaseExperiment


class AtomLoadingOptimizer(EnvExperiment):

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
        self.setattr_argument("n_measurements", NumberValue(100, type='int', scale=1, ndecimals=0, step=1)) # may be necessary to up this if loading rate < 1%
        group = "coil volts boundaries"
        self.setattr_argument("V_AZ_min", NumberValue(3.8,unit="V"), group)
        self.setattr_argument("V_AZ_max", NumberValue(4.5,unit="V"), group)
        self.setattr_argument("V_AX_min", NumberValue(-0.2,unit="V"), group)
        self.setattr_argument("V_AX_max", NumberValue(0.2,unit="V"), group)
        self.setattr_argument("V_AY_min", NumberValue(-0.2, unit="V"), group)
        self.setattr_argument("V_AY_max", NumberValue(0.2, unit="V"), group)

        # this should be close to the mean signal from the atom
        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        self.sampler_buffer = np.zeros(8)

        self.coil_volts_default = [self.AZ_bottom_volts_MOT,
                                   self.AZ_top_volts_MOT,
                                   self.AX_volts_MOT,
                                   self.AY_volts_MOT]

        self.coil_volts_boundaries = [[self.V_AZ_min,self.V_AZ_max],
                                      [self.V_AZ_min, self.V_AZ_max],
                                      [self.V_AX_min,self.V_AX_max],
                                      [self.V_AY_min,self.V_AY_max]]

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

        result = sp.optimize.minimize(fun=self.optimization_routine,
                                      x0=self.coil_volts_default,
                                      bounds=self.coil_volts_boundaries,
                                      method='Nelder-Mead',
                                      options={'xatol':0.003*V # don't try to tune the coils finer than this
                                               } # tolerance in the cost function
                                      )
        print("optimization successful?",result.success)
        print("best coil values (AZ_bottom,AZ_top,AX,AY):",result.x)

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

        atom_events = 0
        for x in data:
            if x > self.atom_counts_threshold:
                atom_events += 1

        if atom_events > 0:
            self.print_async("optimizer found the single atom signal!")

        return -1*atom_events

    @kernel
    def optimization_routine(self, coil_values: TArray(TFloat)) -> TInt32:

        self.core.reset()
        delay(1*ms)

        self.zotino0.set_dac(coil_values, channels=self.coil_channels)
        delay(1*ms)

        self.laser_stabilizer.run()
        delay(1*ms)
        self.dds_FORT.sw.on()

        delay(self.t_MOT_loading)

        for i in range(self.n_measurements):
            t_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
            counts_per_s = self.ttl0.count(t_end)/self.t_SPCM_exposure
            delay(1 * ms)
            self.append_to_dataset(self.count_rate_dataset, counts_per_s)
            self.counts_list[i] = counts_per_s*self.t_SPCM_exposure

        cost = self.get_cost(self.counts_list)
        self.append_to_dataset(self.cost_dataset, cost)

        return cost



