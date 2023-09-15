"""
Connect the output of the SPCM to TTL0 and the output of Zotino0 ch=6 to an oscilloscope. The voltage seen on the scope
corresponds to counts on the SPCM.

Change Sat1s to change the sensitivity: 10**3 with dt=1s gives 5V signal on the scope for 1000 counts/s.

"""

from artiq.experiment import *
from datetime import datetime as dt
import math
import csv

#### Connect the SPCM to ttl0. This code counts and prints the number of photons received per exptime=50ms, for example.
# class SPCMCount(EnvExperiment):
#
#    def build(self):
#        self.setattr_device("core")
#        self.setattr_device("ttl0")
#
#
#    @kernel
#    def run(self):
#        self.core.reset()
#        self.ttl0.input()
#
#        exptime = 50 * ms
#        delay(1 * us)
#        for x in range(100):
#            tend1 = self.ttl0.gate_rising(exptime)
#            count1 = self.ttl0.count(tend1)
#            print(count1)
#            delay(10 * ms)
#
#        print("code done!")


class SPCMCount(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl0")
        self.setattr_device("zotino0")
        self.setattr_device("sampler0")


        self.setattr_argument("n_steps", NumberValue(100, type='int', ndecimals=0, scale=1, step=1))  # exposure time of the SPCM
        self.setattr_argument("dt_exposure", NumberValue(300*ms))  # saturation limit of the SPCM in counts/s. Can be increased to 10**7 safely, but not higher than 3*10**7.
        self.setattr_argument("sat1s", NumberValue(1*10**5), "# of counts giving 5V output. do not set above 10**7") # saturation limit in counts/dt.
        self.setattr_argument("print_count_rate", BooleanValue(True))
        self.setattr_argument("Calculate_average_rate", BooleanValue(True))
        self.setattr_argument("record_counts", BooleanValue(True),"Record counts")
        self.setattr_argument("datadir", StringValue('C:\\Networking Experiment\\artiq codes\\artiq-master\\results\\'),
                              "Record counts")
        self.setattr_argument("datafile", StringValue('SPCMCounts.csv'), "Record counts")
        self.setattr_argument("prepend_date_to_datafile", BooleanValue(True), "Record counts")


    def prepare(self):
        # where to store the data
        self.t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
        if self.prepend_date_to_datafile:
            self.datafile = self.datadir + self.t_experiment_run + '_' + self.datafile
        else:
            self.datafile = self.datadir + self.datafile

        self.sampler_buffer = [0.0]*8

    @kernel
    def run(self):

        self.file_setup(rowheaders=['counts','cooling_volts'])

        self.core.reset()
        self.ttl0.input()
        self.zotino0.init()

        delay(10 * ms)

        ch = 6
        self.zotino0.write_dac(ch, 0.0)
        self.zotino0.load()

        # dt = 300 * ms #exposure time of the SPCM
        # Sat1s = 1*10**5 #saturation limit of the SPCM in counts/s. Can be increased to 10**7 safely, but not higher than 3*10**7.
        # Satdt = Sat1s * dt # saturation limit in counts/dt.
        # delay(1000 * ms)

        Satdt = self.sat1s * self.dt_exposure  # saturation limit in counts/dt.
        delay(1000 * ms)

        CountList = [0.0] * self.n_steps
        # self.core.break_realtime()

        for x in range(self.n_steps):
            tend1 = self.ttl0.gate_rising(self.dt_exposure)
            count1 = self.ttl0.count(tend1)
            if self.print_count_rate:
                print(round(count1/self.dt_exposure),"Hz")
            delay(10 * ms)
            volt1 = count1 * 5/Satdt # the voltage from zotino0, port 7. Saturation limit corresponds to 5V.
            self.zotino0.write_dac(ch, volt1)
            self.zotino0.load()
            delay(1 * ms)

            if self.record_counts:
                ### all items in the list must be the same type that is why 0.0 is added to counts
                self.sampler0.sample(self.sampler_buffer)  # check cooling laser power
                self.file_write([count1 + 0.0,
                                 self.sampler_buffer[7]])
            delay(10 * ms)

            if self.Calculate_average_rate:
                CountList[x] = round(count1/self.dt_exposure) + 0.0

        self.zotino0.write_dac(ch, 0.0)
        self.zotino0.load()

        if self.Calculate_average_rate:
            ### Calculate the sum:
            CountSum = 0.0
            for cc in CountList:
                CountSum += cc

            ### Calculate the average
            AveCount = CountSum / len(CountList)
            print("Average count = ", AveCount, "Hz")


        print("code done!")

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






