"""
This code allows for monitoring an SPCM whose output is connected to TTL0 ch0. The counts detected per s can be viewed
 with the plot_xyline applet (nicknamed "SPCM count rate" in the Node 1 ARTIQ dashboard). Any Zotino channels and
 Urukul channels that were on before running this code will be left on.
"""

from artiq.experiment import *

import sys, os
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment

import numpy as np

class MonitorSPCMinApplet(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("run_time_minutes", NumberValue(1))
        self.setattr_argument("print_count_rate", BooleanValue(False))

        self.base.set_datasets_from_gui_args()

        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()

        self.n_steps = int(60*self.run_time_minutes/self.t_SPCM_exposure+0.5)
        print(self.n_steps)

        self.set_dataset(self.count_rate_dataset,
                             [0.0],
                             broadcast=True)
        #
        # # testing GRIN1_monitor
        # self.set_dataset("GRIN1_monitor", [0.0], broadcast=True)

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware(turn_off_dds_channels=False, turn_off_zotinos=False)

        delay(10 * ms)

        for i in range(self.n_steps):

            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
            count1 = self.ttl0.count(t_gate_end)
            count_rate_per_s = count1 / self.t_SPCM_exposure
            if self.print_count_rate:
                print(round(count_rate_per_s))
            delay(10 * ms)
            self.append_to_dataset(self.count_rate_dataset, count_rate_per_s)

            # # measure_GRIN1(self)
            # measurement_buf = np.array([0.0] * 8)
            # measurement = 0.0
            # avgs = 50
            # for i in range(avgs):
            #     self.sampler1.sample(measurement_buf)
            #     measurement += measurement_buf[self.GRIN1_sampler_ch]
            #     delay(0.1 * ms)
            # measurement /= avgs
            # self.append_to_dataset("GRIN1_monitor", measurement)


        print("Experiment finished.")
