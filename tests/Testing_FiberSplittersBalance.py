"""
This code simply records the voltage of three Sampler channels corresponding to three PD measuring the laser power
after a couple of fiber-splitters.
"""
from artiq.experiment import *
import numpy as np
import csv
from datetime import datetime as dt

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from utilities.BaseExperiment import BaseExperiment

class Test_FiberSplitters(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_steps",
                              NumberValue(100, type='int', ndecimals=0, scale=1, step=1))

        self.setattr_argument("LoopDelay", NumberValue(10 * s)," how often to run the loop in (s)")

        self.setattr_argument("datadir", StringValue("C:\\Networking Experiment\\artiq codes\\artiq-master\\results\\"))

        self.setattr_argument("datafile", StringValue('LaserPowers.csv'), "Record counts")

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.base.prepare()
        self.t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
        self.datafile = self.datadir + self.t_experiment_run + '_' + self.datafile
        self.sampler_buffer = [0.0] * 8

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

        self.file_setup(rowheaders=['P0', 'P1', 'P2', 'P3', 'P4'])

        delay(10 * ms)
        self.dds_cooling_DP.sw.on()

        self.ttl_pumping_repump_switch.on()

        delay(10 * ms)

        n_channels = 8
        iters = 50

        self.core.break_realtime()

        for j in range(self.n_steps):
            print("Loop #: ", j)

            ## Sampler reading with averaging:
            for i in range(n_channels):
                self.smp[i] = 0.0
                self.avg[i] = 0.0

            for i in range(iters):
                self.sampler0.sample(self.smp)  # runs sampler and saves to list
                self.avg += self.smp
                delay(1 * ms)
            self.avg /= iters

            delay(10 * ms)

            self.file_write([self.avg[0], self.avg[1], self.avg[2], self.avg[3], self.avg[4]])

            ### Wait several seconds (e.g. 10s) for the next loop
            time3 = now_mu()
            LoopDelay_mu = self.core.seconds_to_mu(self.LoopDelay)
            delay(self.LoopDelay)
            self.core.wait_until_mu(time3 + LoopDelay_mu)
            delay(10 * ms)

            self.core.reset()
            delay(10 * ms)


        ### The loop for Sampler reading without averaging:
        # for j in range(self.n_steps):
        #     print("Loop #: ", j)
        #
        #     self.sampler0.sample(self.sampler_buffer)
        #     self.file_write([self.sampler_buffer[0], self.sampler_buffer[1], self.sampler_buffer[2], self.sampler_buffer[3]])
        #     delay(10 * ms)
        #
        #     ### Wait several seconds (e.g. 10s) for the next loop
        #     time3 = now_mu()
        #     LoopDelay_mu = self.core.seconds_to_mu(self.LoopDelay)
        #     delay(self.LoopDelay)
        #     self.core.wait_until_mu(time3 + LoopDelay_mu)
        #     delay(10 * ms)

        print("*****************************  ALL DONE  *****************************")