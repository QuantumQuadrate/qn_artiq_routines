"""
This code allows for monitoring both SPCMs which are connected to ttl0 and ttl1. The counts detected per s can be viewed
 with the plot_xyline applets (nicknamed "SPCM0 count rate" and "SPCM1 count rate" in the Node 1 ARTIQ dashboard).
 Any Zotino channels and Urukul channels that were on before running this code will be left on. For long exposure times
 (in ms regime for example) we need to use _counter.gate_rising rather than ttl.count to avoid overflow.
"""

from artiq.experiment import *

import sys, os
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment

class MonitorSPCMinApplet(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("run_time_minutes", NumberValue(1))
        self.setattr_argument("t_SPCM_exposure", NumberValue(0.05))

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
        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware(turn_off_dds_channels=False, turn_off_zotinos=False)

        self.set_dataset(self.SPCM0_rate_dataset, [0.0], broadcast=True)
        self.set_dataset(self.SPCM1_rate_dataset, [0.0], broadcast=True)

        delay(10 * ms)

        for i in range(self.n_steps):
            with parallel:
                self.ttl_SPCM0_counter.gate_rising(self.t_SPCM_exposure)
                self.ttl_SPCM1_counter.gate_rising(self.t_SPCM_exposure)
            SPCM0_counts = self.ttl_SPCM0_counter.fetch_count()
            SPCM1_counts = self.ttl_SPCM1_counter.fetch_count()

            delay(1 * ms)
            self.append_to_dataset(self.SPCM0_rate_dataset, SPCM0_counts / self.t_SPCM_exposure)
            self.append_to_dataset(self.SPCM1_rate_dataset, SPCM1_counts / self.t_SPCM_exposure)

        print("Experiment finished.")
