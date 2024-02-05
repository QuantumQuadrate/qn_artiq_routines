"""
This code turns on and off the fiber AOMs one by one and sends a triggering signal to save series of images with
Andor_Luca. At the same time, record the laser power going into the experiment box, as well as the TA power after
coupling into a fiber. The goal is to analyze the images to see how much each the power in each MOT arm fluctuates
over seconds, minutes, and hours.
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
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')

from utilities.BaseExperiment import BaseExperiment

class CheckMOTBalance(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_steps",
                              NumberValue(10, type='int', ndecimals=0, scale=1, step=1))

        self.setattr_argument("LoopDelay", NumberValue(1.0 * s, ndecimals=3, unit='s')," how often to run the loop in (s)")

        # the Luca will be setup here using pylablib and the analysis will be carried
        # out in self.analyze
        self.setattr_argument("control_Luca_in_software", BooleanValue(True))

        self.setattr_argument("datadir", StringValue('C:\\Users\\QC\\OneDrive - UW-Madison\\Pictures\\'
                                                     'Andor Luca images\\'), "Record counts")

        self.setattr_argument("datafile", StringValue('LaserPowers.csv'), "Record counts")

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.base.prepare()
        self.ddsList = [self.dds_AOM_A1, self.dds_AOM_A2, self.dds_AOM_A3, self.dds_AOM_A4]
        self.datafile = self.datadir + self.datafile

        n_channels = 8
        self.smp = np.zeros(n_channels, dtype=float)
        self.avg = np.zeros(n_channels, dtype=float)


        self.laser_stabilizer.leave_AOMs_on = False

        self.i2V_channels = [7,5,3,4] # the sampler channels for i2Vs PD1,2,3,4
        self.mot_beam_voltages = np.zeros(4) # store an average for each of the chip beams

        if self.control_Luca_in_software:
            self.setup_Luca()
            self.frame_list = [] # empty list should be okay off the kernel

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

    @rpc(flags={"async"})
    def file_setup(self, rowheaders=[]):
        # open file once, then close it at end of the experiment
        self.file_obj = open(self.datafile, 'a', newline='')
        self.csvwriter = csv.writer(self.file_obj)
        if rowheaders != []:
            self.csvwriter.writerow(rowheaders)

    @rpc(flags={"async"})
    def file_write(self, data):
        self.csvwriter.writerow(data)
        self.file_obj.flush()

    @kernel
    def run(self):
        self.base.initialize_hardware()

        acquired = 0

        if not self.control_Luca_in_software:
            print("*****************************  REMEMBER TO START THE CAMERA ACQUISITION  *****************************")
            delay(1000*ms)

        self.file_setup(rowheaders=['cooling_1percent', 'cooling_99percent', '6th_MOT'])

        delay(10 * ms)
        self.dds_cooling_DP.sw.on()

        self.laser_stabilizer.run()
        # self.laser_stabilizer.run(monitor_only=True) # don't actually feedback

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

            delay(1*ms)

            self.file_write([self.avg[7], self.avg[0], self.avg[1]])

            delay(10 * ms)

            ### trigger for Andor_Luca camera by TTL:
            self.ttl6.pulse(5 * ms)


            ### time to wait for camera to take the image
            time2 = now_mu()
            tdelay = 300 * ms
            tdelay_mu = self.core.seconds_to_mu(tdelay)
            delay(tdelay)  # moves the cursor into the future
            self.core.wait_until_mu(time2 + tdelay_mu)
            delay(10 * ms)

            self.laser_stabilizer.run()
            delay(10*ms)

            for i in range(4):
                self.core.reset()

                ### Turn on ith fiber AOM
                self.ddsList[i].sw.on()

                ### this is necessary to have the triggering signal after a certain delay. Otherwise, we do not get a trig signal.
                time1 = now_mu()
                tdelay = 300 * ms
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
                #     acquired = self.get_frame()
                #     delay(1*ms)

                self.sampler0.sample(self.smp)
                self.mot_beam_voltages[i] += self.smp[self.i2V_channels[i]]
                tdelay_mu = self.core.seconds_to_mu(tdelay)
                delay(tdelay)  # moves the cursor into the future
                self.core.wait_until_mu(time2 + tdelay_mu)
                delay(1 * ms)

                ### Turn off ith fiber AOM
                self.ddsList[i].sw.off()

                delay(10 * ms)

            ### Wait several seconds (e.g. 10s) for the next loop... why??
            time3 = now_mu()
            LoopDelay_mu = self.core.seconds_to_mu(self.LoopDelay)
            delay(self.LoopDelay)
            self.core.wait_until_mu(time3 + LoopDelay_mu)
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
            self.cam.stop_acquisition()
            self.cam.close()

            # analysis here
            plt.imshow(self.frame_list[0][0,:,:])
            plt.show()

            fig,axes = plt.subplots(ncols = 5)
            for ax,frame in zip(axes.flat(),self.frames_list):
                ax.imshow(frame[0,:,:])
            plt.imshow()