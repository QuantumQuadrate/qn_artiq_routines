"""
To Monitor the magnetometer connected to Sampler2.

"""

from artiq.experiment import *
import numpy as np
from datetime import datetime as dt

import sys, os
### get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment
from subroutines.experiment_functions import measure_Magnetometer

class MonitorMagnetometer(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(100, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("n_average", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("t_step_ms", NumberValue(10, type='int', scale=1, ndecimals=0, step=1))
        self.setattr_argument("Coils_settings", EnumerationValue(['MOT', 'Optical_pumping', 'Zero_volts', 'None']))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.base.initialize_datasets()
        self.expt()
        print("*************   Experiment finished   *************")

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """

        # measure_Magnetometer(self) ### works, but I want to have more control

        ### x,y, and z axes are connected to Sampler2 Ch2,3, and 4, respectively.

        self.set_dataset("Magnetometer_X", [0.0], broadcast=True)
        self.set_dataset("Magnetometer_Y", [0.0], broadcast=True)
        self.set_dataset("Magnetometer_Z", [0.0], broadcast=True)

        if self.Coils_settings == "MOT":
            ### Set the coils to MOT loading setting
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                channels=self.coil_channels)
            delay(1 * ms)

        if self.Coils_settings == "Optical_pumping":
            ### Set the coils to OP values
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_OP, self.AZ_top_volts_OP, self.AX_volts_OP, self.AY_volts_OP],
                channels=self.coil_channels)
            delay(1 * ms)

        if self.Coils_settings == "Zero_volts":
            ### Turn off all the coils
            self.zotino0.set_dac(
                [0.0, 0.0, 0.0, 0.0],
                channels=self.coil_channels)
            delay(1 * ms)

        if self.Coils_settings == "None":
            ### Does not change the coils values.
            pass
            delay(1 * ms)

        delay(5 * ms)  ### required to avoid a kick in the first measurement

        for n_measrement in range(self.n_measurements):
            measurement_buf = np.array([0.0] * 8)
            MagnetometerX = 0.0
            MagnetometerY = 0.0
            MagnetometerZ = 0.0

            for i in range(self.n_average):
                self.sampler2.sample(measurement_buf)
                MagnetometerX += measurement_buf[self.Magnetometer_X_ch]
                MagnetometerY += measurement_buf[self.Magnetometer_Y_ch]
                MagnetometerZ += measurement_buf[self.Magnetometer_Z_ch]
                delay(0.1 * ms)
            MagnetometerX /= self.n_average
            MagnetometerY /= self.n_average
            MagnetometerZ /= self.n_average
            self.append_to_dataset("Magnetometer_X", MagnetometerY * 350) ### 1V corresponds to 350 mG
            self.append_to_dataset("Magnetometer_Y", MagnetometerX * 350) ### sensor's X axis is coils' Y axis, and vice versa.
            self.append_to_dataset("Magnetometer_Z", MagnetometerZ * 350)
            delay(self.t_step_ms * ms)

        ### finally, in case the worker refuses to die
        self.write_results()