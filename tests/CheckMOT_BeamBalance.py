"""
This code turns on and off the on-chip fiber AOMs one by one and sends a triggering signal to save series of images
with Andor_Luca. This can be used, e.g., to measure the on-chip MOT beam imbalance by looking at the fluorescence
from each MOT beam in the chamber.

The option to control the Luca in software currently does not work. It does not seem that spooling setup is
accessible within pylablib.
"""
from artiq.experiment import *
import logging
import numpy as np
import csv
from time import sleep

import matplotlib.pyplot as plt
import pylablib as pll
pll.par["devices/dlls/andor_sdk2"] = "C:\\Program Files\\Andor SOLIS"
from pylablib.devices import Andor

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment

class CheckMOTBalance(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_steps",
                              NumberValue(10, type='int', ndecimals=0, scale=1, step=1))

        self.setattr_argument("loop_delay", NumberValue(0.001 * s, ndecimals=3, unit='s'))
        self.setattr_argument("t_Luca_exposure", NumberValue(0.5 * s, ndecimals=3, unit='s'))

        # the Luca will be setup here using pylablib and the analysis will be carried
        # out in self.analyze
        self.setattr_argument("control_Luca_in_software", BooleanValue(False))
        self.setattr_argument("enable_feedback", BooleanValue(True))

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.base.prepare()
        self.ddsList = [self.dds_AOM_A1, self.dds_AOM_A2, self.dds_AOM_A3, self.dds_AOM_A4]

        n_channels = 8
        self.smp = np.zeros(n_channels, dtype=float)
        self.avg = np.zeros(n_channels, dtype=float)

        self.laser_stabilizer.leave_AOMs_on = False
        self.monitor_only = not self.enable_feedback

        self.i2V_channels = [0,5,3,4] # the sampler channels for i2Vs PD1,2,3,4
        self.mot_beam_voltages = np.zeros(4) # store an average for each of the chip beams

        if self.control_Luca_in_software:
            self.setup_Luca()
            self.frame_list = [] # empty list should be okay off the kernel
        self.shots_per_measurement = len(self.ddsList)+1 # +1 accounts for the dark image

        self.MOT1_i2V_list = [0.0]
        self.MOT2_i2V_list = [0.0]
        self.MOT3_i2V_list = [0.0]
        self.MOT4_i2V_list = [0.0]

    def setup_Luca(self):
        logging.info("found " + str(Andor.get_cameras_number_SDK2()) + " camera(s)")
        self.cam = Andor.AndorSDK2Camera(temperature=-20, fan_mode="full")

        # warm up
        while not self.cam.get_temperature_status() == "stabilized":
            sleep(0.1)
        print("Andor Luca temperature stabilized")
        # self.cam.set_trigger_mode("int") # I succeeded in taking an image with this
        self.cam.set_trigger_mode("ext_start")
        self.cam.set_acquisition_mode("single")
        self.cam.set_exposure(0.001) # units s
        self.cam.set_EMCCD_gain(100)
        self.cam.start_acquisition()

    def get_frame(self) -> TInt32:
        acquired = 0
        # while acquired == 0:
        try:
            frame = np.array(self.cam.grab())
            self.frame_list.append(frame)
            logging.info("frame acquired")
        except Exception as e:
            logging.info(e)

        (acquired, unread, skipped, size) = self.cam.get_frames_status()
        # for stat, lbl in zip((acquired, unread, skipped, size),
        #                      ("acquired", "unread", "skipped", "size")):
        #     print(lbl, ": ", stat)

        return acquired

    @kernel
    def run(self):
        self.base.initialize_hardware()

        acquired = 0

        if not self.control_Luca_in_software:
            print("*****************************  REMEMBER TO START THE CAMERA ACQUISITION  *****************************")
            delay(1000*ms)

        delay(10 * ms)
        self.dds_cooling_DP.sw.on()

        # warm up
        for i in range(10):
            self.laser_stabilizer.run(monitor_only=self.monitor_only)
            delay(1*s)

        delay(1000 * ms)

        n_channels = 8
        iters = 50
        for j in range(self.n_steps):

            print("Loop #: ", j)
            delay(10 * ms)

            ### Sampler reading with averaging:
            for i in range(n_channels):
                self.smp[i] = 0.0
                self.avg[i] = 0.0

            for i in range(iters):
                self.sampler0.sample(self.smp)  # runs sampler and saves to list
                self.avg += self.smp
                delay(1 * ms)
            self.avg /= iters
            self.print_async(self.avg)

            delay(10 * ms)

            ### trigger for Andor_Luca camera by TTL - takes a dark image
            self.ttl6.pulse(5 * ms)

            ### time to wait for camera to take the image
            time2 = now_mu()
            tdelay = self.t_Luca_exposure + 200 * ms # extra delay just to be safe
            tdelay_mu = self.core.seconds_to_mu(tdelay)
            delay(tdelay)  # moves the cursor into the future
            self.core.wait_until_mu(time2 + tdelay_mu)
            delay(10 * ms)

            self.laser_stabilizer.run(monitor_only=self.monitor_only)
            delay(10*ms)

            self.dds_AOM_A1.sw.off()
            self.dds_AOM_A2.sw.off()
            self.dds_AOM_A3.sw.off()
            self.dds_AOM_A4.sw.off()
            self.dds_AOM_A5.sw.off()
            self.dds_AOM_A6.sw.off()

            for i in range(len(self.ddsList)):
                # self.core.reset()

                ### Turn on ith fiber AOM
                self.ddsList[i].sw.on()

                ### this is necessary to have the triggering signal after a certain delay. Otherwise, we do not get a trig signal.
                time1 = now_mu()
                tdelay = self.t_Luca_exposure + 200 * ms # extra delay just to be safe
                tdelay_mu = self.core.seconds_to_mu(tdelay)
                delay(tdelay)  # moves the cursor into the future
                self.core.wait_until_mu(time1 + tdelay_mu)  # wait for the cursor to get there
                delay(1 * ms)

                ### trigger for Andor_Luca camera by TTL:
                self.ttl6.pulse(5 * ms)

                ### time to wait for camera to take the image
                time2 = now_mu()
                tdelay = 300 * ms

                # if self.control_Luca_in_software:
                #     self.core.break_realtime()
                #     delay(1*s)
                #     acquired = self.get_frame()
                #     delay(100*ms)

                self.sampler0.sample(self.smp)
                self.mot_beam_voltages[i] += self.smp[self.i2V_channels[i]]
                tdelay_mu = self.core.seconds_to_mu(tdelay)
                delay(tdelay)  # moves the cursor into the future
                self.core.wait_until_mu(time2 + tdelay_mu)
                delay(1 * ms)

                ### Turn off ith fiber AOM
                self.ddsList[i].sw.off()

                delay(10 * ms)

            time3 = now_mu()
            loop_delay_mu = self.core.seconds_to_mu(self.loop_delay)
            delay(self.loop_delay)
            self.core.wait_until_mu(time3 + loop_delay_mu)
            delay(10 * ms)

        print("MOT beam i2V voltages (1-4):")
        print(self.mot_beam_voltages/self.n_steps)

        print("*****************************  ALL DONE  *****************************")

    def analyze(self):

        if self.control_Luca_in_software:
            """
            It is important to close all camera connections before finishing
            your script. Otherwise, DLL resources might become permanently
            blocked, and the only way to solve it would be to restart the PC.
            - pylablib docs
            """

            try:
                self.frame_array = np.array(self.cam.grab(nframes=self.n_steps*self.shots_per_measurement))
                logging.info("frame acquired")
            except Exception as e:
                logging.info(e)

            (acquired, unread, skipped, size) = self.cam.get_frames_status()
            for stat, lbl in zip((acquired, unread, skipped, size),
                                 ("acquired", "unread", "skipped", "size")):
                print(lbl, ": ", stat)

            self.cam.stop_acquisition()
            self.cam.close()

            print(self.frame_array.shape)
            #
            # # analysis here
            # plt.imshow(self.frame_list[0][0,:,:])
            # plt.show()
            #
            # fig,axes = plt.subplots(ncols = 5)
            # for ax,frame in zip(axes.flat(),self.frames_list):
            #     ax.imshow(frame[0,:,:])
            # plt.imshow()