"""
This code turns on and off the fiber AOMs one by one and sends a triggering signal to save series of images with
Andor_Luca. At the same time, record the laser power going into the experiment box, as well as the TA power after
coupling into a fiber. The goal is to analyze the images to see how much each the power in each MOT arm fluctuates
over seconds, minutes, and hours.
"""
from artiq.experiment import *
import numpy as np
import csv

import sys, os
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')

from utilities.BaseExperiment import BaseExperiment

class CheckMOTBalance(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_steps",
                              NumberValue(100, type='int', ndecimals=0, scale=1, step=1))

        self.setattr_argument("LoopDelay", NumberValue(10 * s)," how often to run the loop in (s)")

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

        self.file_setup(rowheaders=['cooling_1percent', 'cooling_99percent', '6th_MOT'])

        delay(10 * ms)
        self.dds_cooling_DP.sw.on()
        self.dds_cooling_SP.sw.on()
        self.dds_MOT_RP.sw.off()
        self.dds_AOM_A5.sw.off()

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

        n_channels = 8
        iters = 50
        for j in range(self.n_steps):
            print("Loop #: ", j)
            delay(10 * ms)

            self.dds_AOM_A6.sw.on()
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

            delay(1*ms)

            self.file_write([self.avg[7], self.avg[0], self.avg[1]])



            # ### Sampler reading without averaging:
            # smp = [0.0] * n_channels
            # self.sampler0.sample(smp)  # runs sampler and saves to list
            # self.file_write([smp[7], smp[0], smp[1]])

            delay(10 * ms)

            self.dds_AOM_A6.sw.off()
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
            tdelay = 100 * ms
            tdelay_mu = self.core.seconds_to_mu(tdelay)
            delay(tdelay)  # moves the cursor into the future
            self.core.wait_until_mu(time2 + tdelay_mu)
            delay(10 * ms)

            for i in range(4):
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
                tdelay = 100 * ms
                tdelay_mu = self.core.seconds_to_mu(tdelay)
                delay(tdelay)  # moves the cursor into the future
                self.core.wait_until_mu(time2 + tdelay_mu)
                delay(1 * ms)

                ### Turn off ith fiber AOM
                self.ddsList[i].sw.off()

                delay(10 * ms)

            ### Wait several seconds (e.g. 10s) for the next loop
            time3 = now_mu()
            LoopDelay_mu = self.core.seconds_to_mu(self.LoopDelay)
            delay(self.LoopDelay)
            self.core.wait_until_mu(time3 + LoopDelay_mu)
            delay(10 * ms)

        print("*****************************  ALL DONE  *****************************")