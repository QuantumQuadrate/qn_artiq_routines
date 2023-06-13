"""
This code simply records the voltage of three Sampler channels corresponding to three PD measuring the laser power
after a couple of fiber-splitters.
"""
from artiq.experiment import *
import math
import numpy as np
import csv
from datetime import datetime as dt
from subroutines.stabilizer import AOMPowerStabilizer
import DeviceAliases

from BaseExperiment import BaseExperiment

class Test_FiberSplitters(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_steps",
                              NumberValue(100, type='int', ndecimals=0, scale=1, step=1))

        self.setattr_argument("LoopDelay", NumberValue(10 * s)," how often to run the loop in (s)")

        self.setattr_argument("datadir", StringValue("C:\\Networking Experiment\\artiq codes\\artiq-master\\results\\"))

        self.setattr_argument("datafile", StringValue('LaserPowers.csv'), "Record counts")


    def prepare(self):
        self.base.prepare()
        self.t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
        self.datafile = self.datadir + self.t_experiment_run + '_' + self.datafile
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

        self.file_setup(rowheaders=['P1_1perc_Tap', 'P2_99perc_Tap', 'P3_InTBox'])

        delay(10 * ms)
        self.dds_cooling_DP.sw.on()
        self.dds_cooling_SP.sw.on()
        self.dds_MOT_RP.sw.off()

        delay(1000 * ms)

        for j in range(self.n_steps):
            print("Loop #: ", j)

            self.sampler0.sample(self.sampler_buffer)
            self.file_write([self.sampler_buffer[7], self.sampler_buffer[0], self.sampler_buffer[1]])
            delay(10 * ms)

            ### Wait several seconds (e.g. 10s) for the next loop
            time3 = now_mu()
            LoopDelay_mu = self.core.seconds_to_mu(self.LoopDelay)
            delay(self.LoopDelay)
            self.core.wait_until_mu(time3 + LoopDelay_mu)
            delay(10 * ms)

        print("*****************************  ALL DONE  *****************************")