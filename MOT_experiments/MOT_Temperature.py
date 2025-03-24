"""
This code loads a MOT and lets it expand for various times t before taking an image
The resulting images can be fit to Gaussians to extract the temperature.

This particular code outputs a trigger to a camera, but the camera setup (spooling, etc)
and image processing are done elsewhere for now.

If using the ThorCam, run
C:\\Networking Experiment\\Camera_Examples\\Python\\save_ims_with_hardware_trigger.py
If using Luca, run Andor Solis, enable spooling and fast external triggering,
with exposure time=1 ms.
"""

from artiq.experiment import *
import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment


class MOTTemperature(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # experiment variables which are specific to this experiment
        # todo: add tif directory so the script will automatically find the images
        # from the Luca and fit them

        self.setattr_argument("release_times_ms", StringValue(
                '[0.001,0.005,0.01,0.05, 0.5]'))
        self.setattr_argument("camera_model", EnumerationValue(['Luca','ThorCam'],default='ThorCam'))
        # wait for the camera to take the shot
        self.setattr_argument("t_exposure", NumberValue(0.5*ms, unit='ms'))
        self.setattr_argument("averages", NumberValue(70, type='int', precision=0, scale=1, step=1))
        self.setattr_argument("do_PGC_in_MOT", BooleanValue(False))

        self.base.set_datasets_from_gui_args()


    def prepare(self):
        self.base.prepare()

        self.release_times_ms = eval(self.release_times_ms)
        # convert to ms and make sure everything is a float
        self.release_times_ms = [float(x)*ms for x in self.release_times_ms]
        self.n_steps = len(self.release_times_ms)

    # @kernel
    # def load_MOT(self, PGC=False):
    #     # set coils and DDS to MOT settings
    #     self.zotino0.set_dac(
    #         [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
    #         channels=self.coil_channels)
    #     self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
    #
    #     delay(self.t_MOT_loading)
    #
    #     if PGC:
    #         # todo: for now, use the same amplitude as the MOT. see my note in BaseExperiment.prepare
    #         self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_MOT)
    #         delay(self.t_PGC_in_MOT)

    @kernel
    def run(self):
        self.base.initialize_hardware()

        delay(1*ms)
        self.dds_cooling_DP.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()
        delay(3000*ms) # wait for AOMs to thermalize

        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)

        # trigger Luca to save an image - sanity check to make sure there is a MOT
        self.ttl6.pulse(5 * ms)

        delay(1*ms)
        t = 0.0

        for i in range(self.n_steps):
            for n in range(self.averages):
                t = self.release_times_ms[i]

                # reset coils and DDS to MOT settings
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                    channels=self.coil_channels)
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
                delay(10 * ms)

                delay(10*ms)
                self.laser_stabilizer.run()

                delay(self.t_MOT_loading)

                # self.load_MOT(PGC=self.do_PGC_in_MOT)

                if self.do_PGC_in_MOT:
                    # todo: for now, use the same amplitude as the MOT. see my note in BaseExperiment.prepare
                    self.dds_cooling_DP.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_MOT)
                    delay(self.t_PGC_in_MOT)

                # # trigger Luca to save an image - sanity check to make sure there is a MOT
                # self.ttl6.pulse(5 * ms)
                # # delay(5*ms)

                # drops the MOT. alternatively we could chop the RP, as the coils may be too slow for this
                self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],
                                     channels=self.coil_channels)
                self.dds_cooling_DP.sw.off()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                self.dds_AOM_A6.sw.off()

                delay(t)

                # turn on imaging beams and trigger Luca to take an image
                self.dds_cooling_DP.sw.on()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                self.dds_AOM_A6.sw.on()
                self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=2*self.ampl_cooling_DP_MOT)

                if self.camera_model == 'ThorCam':
                    # trigger the ThorCam with the Zotino
                    self.zotino0.write_dac(6, 4.0)
                    self.zotino0.load()
                    delay(self.t_exposure)
                    self.zotino0.write_dac(6, 0.0)
                    self.zotino0.load()
                else:# assume self.camera_model == 'Luca'
                    self.ttl6.pulse(5*ms)
                delay(self.t_exposure)
                delay(10*ms)

                # # reset coils and DDS to MOT settings
                # self.zotino0.set_dac(
                #     [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                #     channels=self.coil_channels)
                # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
                # delay(10 * ms)

        print("MOT temperature measurement done!")