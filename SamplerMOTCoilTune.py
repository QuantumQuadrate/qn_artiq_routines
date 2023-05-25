"""
This code scans the coils to find the optimum coil parameters that put the MOT at the focus of the
parabolic mirror. It saves the photon counts in a file.

"""

from artiq.experiment import *
import csv
from artiq.coredevice import ad53xx # for converting volts to mu for the zotino
from artiq.coredevice.exceptions import RTIOUnderflow
import math # for math
import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt

from BaseExperiment import BaseExperiment

class SamplerMOTCoilTune(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # todo
        # self.setattr_argument("AZ_bottom_V_limits", StringValue('(-3.5, 3.5)'))
        # self.setattr_argument("AZ_top_V_limits", StringValue('(-1.5, 1.5)'))
        # self.setattr_argument("AZ_X_V_limits", StringValue('(-1, 1)'))
        # self.setattr_argument("AZ_Y_V_limits", StringValue('(-0.5, 0.5)'))

        self.setattr_argument("MOT_AOMs_on", BooleanValue(True))
        self.setattr_argument("print_control_volts", BooleanValue(True))
        self.setattr_argument("update_coil_volts_at_finish", BooleanValue(False))
        self.setattr_argument("n_steps", NumberValue(100, type='int', ndecimals=0, scale=1, step=1))

        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        self.base.prepare()
        self.sampler_buffer = [0.0]*8
        self.control_volts_channels = [0,1,2,3] # the sampler channels to read

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()


        if self.MOT_AOMs_on:
            # Turn on AOMs to load the MOT.
            self.dds_cooling_DP.sw.on()
            self.dds_cooling_SP.sw.on()
            self.dds_MOT_RP.sw.on()
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A6.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()

            # wait for AOMs to thermalize
            delay(3000 * ms)

        control_volts = [0.0]*4

        # self.core.break_realtime()

        for i in range(self.n_steps):

            self.sampler1.sample(self.sampler_buffer)
            control_volts = [self.sampler_buffer[ch] for ch in self.control_volts_channels]

            delay(10*ms)
            if self.print_control_volts:
                print(control_volts)

            # set coils based on the sampler values we read
            self.zotino0.set_dac(
                control_volts,
                channels=self.coil_channels)

            delay(50 * ms)

        if self.update_coil_volts_at_finish:
            volt_datasets = ["AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT"]
            for i in range(4):
                self.set_dataset(volt_datasets[i], control_volts[i])

        print("Experiment finished.")
