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
        self.ttl8.input()
        self.ttl9.input()

        delay(10 * ms)
        self.set_dataset("TTL8_per_s", [0.0], broadcast=True)
        self.set_dataset("TTL9_per_s", [0.0], broadcast=True)

        delay(10 * ms)

        for i in range(self.n_steps):

            t_gate_end = self.ttl8.gate_rising(self.t_SPCM_exposure)
            count8 = self.ttl8.count(t_gate_end)
            count8_per_s = count8 / self.t_SPCM_exposure

            delay(1 * ms)

            t_gate_end = self.ttl9.gate_rising(self.t_SPCM_exposure)
            count9 = self.ttl9.count(t_gate_end)
            count9_per_s = count9 / self.t_SPCM_exposure

            delay(1 * ms)
            self.append_to_dataset("TTL8_per_s", count8_per_s)
            delay(1 * ms)
            self.append_to_dataset("TTL9_per_s", count9_per_s)


        print("Experiment finished.")
