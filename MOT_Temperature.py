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
        self.setattr_argument("averages", NumberValue(10, type='int', ndecimals=0, scale=1, step=1))

        self.base.set_datasets_from_gui_args()


    def prepare(self):
        self.base.prepare()

        self.release_times_ms = eval(self.release_times_ms)
        # convert to ms and make sure everything is a float
        self.release_times_ms = [float(x)*ms for x in self.release_times_ms]
        self.n_steps = len(self.release_times_ms)

    @kernel
    def run(self):
        self.base.initialize_hardware()

        delay(1*ms)
        self.dds_cooling_DP.sw.on()
        self.dds_D1_pumping_SP.sw.on()
        self.dds_pumping_repump.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        delay(3000*ms)

        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)

        # trigger Luca to save an image - sanity check to make sure there is a MOT
        self.ttl6.pulse(5 * ms)

        delay(1*ms)
        t = 0.0

        for i in range(self.n_steps):
            for n in range(self.averages):
                t = self.release_times_ms[i]

                delay(10*ms)
                self.laser_stabilizer.run()

                delay(3000*ms)
                # # trigger Luca to save an image - sanity check to make sure there is a MOT
                # self.ttl6.pulse(5 * ms)
                # # delay(5*ms)

                # drops the MOT. alternatively we could chop the RP, but the coils may be too slow for this
                self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],
                                     channels=self.coil_channels)
                self.dds_cooling_DP.sw.off()
                self.dds_pumping_repump.sw.off()

                delay(t)

                # turn on imaging beams and trigger Luca to take an image
                self.dds_cooling_DP.sw.on()
                self.dds_pumping_repump.sw.on()
                self.ttl6.pulse(5*ms)

                delay(self.t_Luca_exposure)

                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                    channels=self.coil_channels)
                delay(10*ms)

        print("MOT temperature measurement done!")