"""
This code turns on and off the fiber AOMs one by one and sends a triggering signal to save series of images with
Andor_Luca. At the same time, record the laser power going into the experiment box, as well as the TA power after
coupling into a fiber. The goal is to analyze the images to see how much each the power in each MOT arm fluctuates
over seconds, minutes, and hours.
"""
from artiq.experiment import *
import math
import numpy as np
import csv
from subroutines.stabilizer import AOMPowerStabilizer
import DeviceAliases

from BaseExperiment import BaseExperiment

class CheckMOTBallance(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("datadir", StringValue('C:\\Users\\QC\\OneDrive - UW-Madison\\Pictures\\'
                                                     'Andor Luca images\\'), "Record counts")

        self.setattr_argument("datafile", StringValue('LaserPower.csv'), "Record counts")


    def prepare(self):
        self.base.prepare()
        self.ddsList = [self.dds_AOM_A1, self.dds_AOM_A2, self.dds_AOM_A3, self.dds_AOM_A4, self.dds_AOM_A6]
        self.datafile = self.datadir + self.datafile
        self.sampler_buffer = [0.0] * 8


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

        self.file_setup(rowheaders=['cooling_volts'])

        delay(10 * ms)
        self.dds_cooling_DP.sw.on()
        self.dds_cooling_SP.sw.on()
        delay(1000 * ms)

        ### not working!
        # for x in range(2):
        #     for i in range(1,7):
        #         dds_name = "dds_AOM_A" + str(i)
        #         dds_channel = getattr(self, dds_name)
        #         dds_channel.sw.on()
        #         delay(500 * ms)

        ### Working
        # for i in range(1, 6):
        #     self.ddsList[i].sw.on()
        #     delay(500 * ms)


        for j in range(3):
            print("Loop #: ", j)

            self.sampler0.sample(self.sampler_buffer)  # check cooling laser power
            self.file_write([self.sampler_buffer[7]])
            delay(10 * ms)

            ### Dark image
            ### trigger for Andor_Luca camera.
            self.zotino0.write_dac(6, 4.0)
            self.zotino0.load()
            delay(5 * ms)
            self.zotino0.write_dac(6, 0.0)
            self.zotino0.load()

            ### time to wait for camera to take the image
            time2 = now_mu()
            tdelay = 50 * ms
            tdelay_mu = self.core.seconds_to_mu(tdelay)
            delay(tdelay)  # moves the cursor into the future
            self.core.wait_until_mu(time2 + tdelay_mu)
            delay(1 * ms)

            for i in range(5):
                ### Turn on ith fiber AOM
                self.ddsList[i].sw.on()

                ### Delay 100ms
                ### this is necessary to have the triggering signal after a certain delay. Otherwise, we do not get a trig signal.
                time1 = now_mu()
                tdelay = 100 * ms
                tdelay_mu = self.core.seconds_to_mu(tdelay)
                delay(tdelay)  # moves the cursor into the future
                self.core.wait_until_mu(time1 + tdelay_mu)  # wait for the cursor to get there
                delay(1 * ms)

                ### trigger for Andor_Luca camera.
                self.zotino0.write_dac(6, 4.0)
                self.zotino0.load()
                delay(5 * ms)
                self.zotino0.write_dac(6, 0.0)
                self.zotino0.load()

                ### time to wait for camera to take the image
                time2 = now_mu()
                tdelay = 50 * ms
                tdelay_mu = self.core.seconds_to_mu(tdelay)
                delay(tdelay)  # moves the cursor into the future
                self.core.wait_until_mu(time2 + tdelay_mu)
                delay(1 * ms)

                ### Turn off ith fiber AOM
                self.ddsList[i].sw.off()

                delay(10 * ms)

            ### Wait several seconds (e.g. 10s) for the next loop
            time3 = now_mu()
            tdelay = 2 * s
            tdelay_mu = self.core.seconds_to_mu(tdelay)
            delay(tdelay)  # moves the cursor into the future
            self.core.wait_until_mu(time3 + tdelay_mu)
            delay(10 * ms)







        # for x in range(1000):
        #
        #     ### trigger for Andor_Luca camera.
        #     self.zotino0.write_dac(6, 4.0)
        #     self.zotino0.load()
        #     delay(5 * ms)
        #     self.zotino0.write_dac(6, 0.0)
        #     self.zotino0.load()
        #
        #     ### time to wait for camera to take the image
        #     time2 = now_mu()
        #     tdelay = 50 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time2 + tdelay_mu)
        #     delay(1 * ms)
        #
        #     ### Turn on AOM_A1
        #     self.dds_AOM_A1.sw.on()
        #
        #     ### Delay 100ms
        #     ### this is necessary to have the triggering signal after a certain delay. Otherwise, we do not get a trig signal.
        #     time1 = now_mu()
        #     tdelay = 100 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time1 + tdelay_mu)  # wait for the cursor to get there
        #     delay(1 * ms)
        #
        #     ### trigger for Andor_Luca camera.
        #     self.zotino0.write_dac(6, 4.0)
        #     self.zotino0.load()
        #     delay(5 * ms)
        #     self.zotino0.write_dac(6, 0.0)
        #     self.zotino0.load()
        #
        #     ### time to wait for camera to take the image
        #     time2 = now_mu()
        #     tdelay = 50 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time2 + tdelay_mu)
        #     delay(1 * ms)
        #
        #     ### Turn off AOM_A1
        #     self.dds_AOM_A1.sw.off()
        #
        #     delay(10 * ms)
        #     ###############################################################################
        #
        #     ### Turn on AOM_A2
        #     self.dds_AOM_A2.sw.on()
        #
        #     ### Delay 100ms
        #     ### this is necessary to have the triggering signal after a certain delay. Otherwise, we do not get a trig signal.
        #     time1 = now_mu()
        #     tdelay = 100 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time1 + tdelay_mu)  # wait for the cursor to get there
        #     delay(1 * ms)
        #
        #     ### trigger for Andor_Luca camera.
        #     self.zotino0.write_dac(6, 4.0)
        #     self.zotino0.load()
        #     delay(5 * ms)
        #     self.zotino0.write_dac(6, 0.0)
        #     self.zotino0.load()
        #
        #     ### time to wait for camera to take the image
        #     time2 = now_mu()
        #     tdelay = 50 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time2 + tdelay_mu)
        #     delay(1 * ms)
        #
        #     ### Turn off AOM_A2
        #     self.dds_AOM_A2.sw.off()
        #
        #     delay(10 * ms)
        #     ###############################################################################
        #
        #     ### Turn on AOM_A3
        #     self.dds_AOM_A3.sw.on()
        #
        #     ### Delay 100ms
        #     ### this is necessary to have the triggering signal after a certain delay. Otherwise, we do not get a trig signal.
        #     time1 = now_mu()
        #     tdelay = 100 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time1 + tdelay_mu)  # wait for the cursor to get there
        #     delay(1 * ms)
        #
        #     ### trigger for Andor_Luca camera.
        #     self.zotino0.write_dac(6, 4.0)
        #     self.zotino0.load()
        #     delay(5 * ms)
        #     self.zotino0.write_dac(6, 0.0)
        #     self.zotino0.load()
        #
        #     ### time to wait for camera to take the image
        #     time2 = now_mu()
        #     tdelay = 50 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time2 + tdelay_mu)
        #     delay(1 * ms)
        #
        #     ### Turn off AOM_A3
        #     self.dds_AOM_A3.sw.off()
        #
        #     delay(10 * ms)
        #     ###############################################################################
        #
        #     ### Turn on AOM_A4
        #     self.dds_AOM_A4.sw.on()
        #
        #     ### Delay 100ms
        #     ### this is necessary to have the triggering signal after a certain delay. Otherwise, we do not get a trig signal.
        #     time1 = now_mu()
        #     tdelay = 100 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time1 + tdelay_mu)  # wait for the cursor to get there
        #     delay(1 * ms)
        #
        #     ### trigger for Andor_Luca camera.
        #     self.zotino0.write_dac(6, 4.0)
        #     self.zotino0.load()
        #     delay(5 * ms)
        #     self.zotino0.write_dac(6, 0.0)
        #     self.zotino0.load()
        #
        #     ### time to wait for camera to take the image
        #     time2 = now_mu()
        #     tdelay = 50 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time2 + tdelay_mu)
        #     delay(1 * ms)
        #
        #     ### Turn off AOM_A4
        #     self.dds_AOM_A4.sw.off()
        #
        #     delay(10 * ms)
        #     ###############################################################################
        #
        #     ### Turn on AOM_A6
        #     self.dds_AOM_A6.sw.on()
        #
        #     ### Delay 100ms
        #     ### this is necessary to have the triggering signal after a certain delay. Otherwise, we do not get a trig signal.
        #     time1 = now_mu()
        #     tdelay = 100 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time1 + tdelay_mu)  # wait for the cursor to get there
        #     delay(1 * ms)
        #
        #     ### trigger for Andor_Luca camera.
        #     self.zotino0.write_dac(6, 4.0)
        #     self.zotino0.load()
        #     delay(5 * ms)
        #     self.zotino0.write_dac(6, 0.0)
        #     self.zotino0.load()
        #
        #     ### time to wait for camera to take the image
        #     time2 = now_mu()
        #     tdelay = 50 * ms
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time2 + tdelay_mu)
        #     delay(1 * ms)
        #
        #     ### Turn off AOM_A6
        #     self.dds_AOM_A6.sw.off()
        #     ###############################################################################
        #
        #     time3 = now_mu()
        #     tdelay = 10 * s
        #     tdelay_mu = self.core.seconds_to_mu(tdelay)
        #     delay(tdelay)  # moves the cursor into the future
        #     self.core.wait_until_mu(time3 + tdelay_mu)
        #     delay(10 * ms)



        print("*****************************  ALL DONE  *****************************")