"""
This code loads a MOT and lets it expand for various times t before taking an image
The resulting images can be fit to Gaussians to extract the temperature.

This particular code outputs a trigger to a camera, but the camera setup (spooling, etc)
and image processing are done elsewhere for now.
"""
from artiq.experiment import *

from utilities.BaseExperiment import BaseExperiment

class MOTTemperature(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # experiment variables which are specific to this experiment
        # todo: add tif directory so the script will automatically find the images
        # from the Luca and fit them

        self.setattr_argument("release_times_ms", StringValue(
                '[0.001, 0.005, 0.01, 0.05, 0.1, 1, 2, 5]'))
        # wait for the camera to take the shot
        self.setattr_argument("t_Luca_exposure", NumberValue(50*ms, unit='ms'))

    def prepare(self):
        self.base.prepare()

        self.release_times_ms = eval(self.release_times_ms)
        # convert to ms and make sure everything is a float
        self.release_times_ms = [float(x)*ms for x in self.release_times_ms]
        self.n_steps = len(self.release_times_ms)

    @kernel
    def run(self):
        self.base.initialize_hardware()

        # # URUKUL 0 - MOT and D2 state prep AOMs:

        delay(1*ms)
        self.dds_cooling_DP.sw.on()
        self.dds_cooling_SP.sw.on()
        self.dds_MOT_RP.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        delay(3000*ms)

        # trigger Luca to save an image of the background (i.e. no MOT)
        # todo: this shot shows up way darker than the subsequent ones,
        # probably because of the fluorescence change when the coils turn on
        self.zotino0.write_dac(6, 4.0)
        self.zotino0.load()
        delay(5 * ms)
        self.zotino0.write_dac(6, 0.0)
        self.zotino0.load()
        delay(self.t_Luca_exposure)

        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)
        delay(1*ms)
        t = 0.0

        for i in range(self.n_steps):
            t = self.release_times_ms[i]
            self.dds_cooling_DP.sw.on()
            delay(1000*ms)
            self.dds_cooling_DP.sw.off()
            delay(t)

            # trigger Luca to save an image
            self.zotino0.write_dac(6, 4.0)
            self.zotino0.load()
            delay(5 * ms)
            self.zotino0.write_dac(6, 0.0)
            self.zotino0.load()
            delay(self.t_Luca_exposure)

        # leave the MOT on at the end
        self.dds_cooling_DP.sw.on()
        print("MOT temperature measurement done!")