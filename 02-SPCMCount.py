"""
Connect the output of the SPCM to TTL0 and the output of Zotino0 ch=6 to an oscilloscope. The voltage seen on the scope
corresponds to counts on the SPCM.

Change Sat1s to change the sensitivity: 10**3 with dt=1s gives 5V signal on the scope for 1000 counts/s.

"""

from artiq.experiment import *

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

        self.setattr_argument("n_steps", NumberValue(100, type='int', ndecimals=0, scale=1, step=1))  # exposure time of the SPCM
        self.setattr_argument("dt_exposure", NumberValue(300*ms))  # saturation limit of the SPCM in counts/s. Can be increased to 10**7 safely, but not higher than 3*10**7.
        self.setattr_argument("sat1s", NumberValue(1*10**5)) # saturation limit in counts/dt.
        self.setattr_argument("print_counts", BooleanValue(True))

    @kernel
    def run(self):
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

        for x in range(self.n_steps):
            tend1 = self.ttl0.gate_rising(self.dt_exposure)
            count1 = self.ttl0.count(tend1)
            if self.print_counts:
                print(count1)
            delay(10 * ms)
            volt1 = count1 * 5/Satdt # the voltage from zotino0, port 7. Saturation limit corresponds to 5V.
            self.zotino0.write_dac(ch, volt1)
            self.zotino0.load()
            delay(1 * ms)

        self.zotino0.write_dac(ch, 0.0)
        self.zotino0.load()

        print("code done!")






