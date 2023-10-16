"""
To turn on the AOMs of node1 of the networking experiment.

"""
from artiq.experiment import *
import math
# import numpy as np


class Set_AOMs(EnvExperiment):

    def build(self):
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")


        # self.setattr_argument("f_FORT", NumberValue(345.1 * MHz, unit="MHz", ndecimals=1), "AOM1, pumping repumper")
        # self.setattr_argument("p_FORT_loading", NumberValue(-5, unit="dBm", scale=1, ndecimals=1), "AOM1, pumping repumper")
        # self.setattr_argument("FORT_AOM_ON", BooleanValue(default=False), "AOM1, pumping repumper")

        self.setattr_argument("f_FORT", NumberValue(210 * MHz, unit="MHz", ndecimals=1), "AOM1, FORT")
        self.setattr_argument("p_FORT_loading", NumberValue(3, unit="dBm", scale=1, ndecimals=1), "AOM1, FORT")
        self.setattr_argument("FORT_AOM_ON", BooleanValue(default=False), "AOM1, FORT")

        self.setattr_argument("f_cooling_DP_MOT", NumberValue(115.0 * MHz, unit="MHz",ndecimals=1), "AOM2, MOT cooling double pass")
        self.setattr_argument("p_cooling_DP_MOT", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1), "AOM2, MOT cooling double pass")
        self.setattr_argument("Cooling_DP_AOM_ON", BooleanValue(default=False), "AOM2, MOT cooling double pass")

        self.setattr_argument("f_D1_pumping_SP", NumberValue(130.0 * MHz, unit="MHz",ndecimals=1), "AOM3, MOT cooling single pass")
        self.setattr_argument("p_D1_pumping_SP", NumberValue(1, unit="dBm", scale=1, ndecimals=1), "AOM3, MOT cooling single pass")
        self.setattr_argument("D1_pumping_SP_AOM_ON", BooleanValue(default=False), "AOM3, MOT cooling single pass")

        self.setattr_argument("f_pumping_RP", NumberValue(150.5 * MHz, unit="MHz", ndecimals=1), "AOM4, MOT RP/Exc")
        self.setattr_argument("p_pumping_RP", NumberValue(3, unit="dBm", scale=1, ndecimals=1), "AOM4, MOT RP/Exc")
        self.setattr_argument("pumping_RP_AOM_ON", BooleanValue(default=False), "AOM4, MOT RP/Exc")
        # the default power for the fiber AOMs was chosen to give roughly equal diffraction efficiency, empirically
        self.setattr_argument("AOM_A2_freq", NumberValue(78.48 * MHz, unit="MHz", ndecimals=2), "AOM A2")
        self.setattr_argument("AOM_A2_power", NumberValue(-5, unit="dBm", scale=1, ndecimals=1), "AOM A2")
        self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "AOM A2")

        self.setattr_argument("AOM_A3_freq", NumberValue(78.49 * MHz, unit="MHz", ndecimals=2), "AOM A3")
        self.setattr_argument("AOM_A3_power", NumberValue(-3, unit="dBm", scale=1, ndecimals=1), "AOM A3")
        self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "AOM A3")

        self.setattr_argument("AOM_A5_freq", NumberValue(78.5 * MHz, unit="MHz", ndecimals=2), "AOM A5")
        self.setattr_argument("AOM_A5_power", NumberValue(-3, unit="dBm", scale=1, ndecimals=1), "AOM A5")
        self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "AOM A5")

        self.setattr_argument("AOM_A6_freq", NumberValue(78.51 * MHz, unit="MHz", ndecimals=2), "AOM A6")
        self.setattr_argument("AOM_A6_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A6")
        self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "AOM A6")

    def prepare(self):
        # converts RF power in dBm to amplitudes in V
        self.AOM1_ampl = math.sqrt(2*50*10**(self.p_FORT_loading/10-3))
        self.AOM2_ampl = math.sqrt(2*50*10**(self.p_cooling_DP_MOT/10-3))
        self.AOM3_ampl = math.sqrt(2*50*10**(self.p_D1_pumping_SP/10-3))
        self.AOM4_ampl = math.sqrt(2*50*10**(self.p_pumping_RP/10-3))

        self.AOM_A2_ampl = math.sqrt(2*50*10**(self.AOM_A2_power/10-3))
        self.AOM_A3_ampl = math.sqrt(2*50*10**(self.AOM_A3_power/10-3))
        self.AOM_A5_ampl = math.sqrt(2*50*10**(self.AOM_A5_power/10-3))
        self.AOM_A6_ampl = math.sqrt(2*50*10**(self.AOM_A6_power/10-3))


    @kernel
    def run(self):
    #     #initializes the hardware and resets attenuators to 0 dB
        self.core.reset()
        self.urukul0_cpld.init()
        self.urukul1_cpld.init()

        self.urukul0_ch0.init()
        self.urukul0_ch1.init()
        self.urukul0_ch2.init()
        self.urukul0_ch3.init()
        self.urukul1_ch0.init()
        self.urukul1_ch1.init()
        self.urukul1_ch2.init()
        self.urukul1_ch3.init()

        self.urukul0_ch0.set_att(float(0))
        self.urukul0_ch1.set_att(float(0))
        self.urukul0_ch2.set_att(float(0))
        self.urukul0_ch3.set_att(float(0))
        self.urukul1_ch0.set_att(float(0))
        self.urukul1_ch1.set_att(float(0))
        self.urukul1_ch2.set_att(float(0))
        self.urukul1_ch3.set_att(float(0))

        self.core.break_realtime()

        # self.urukul0_cpld.set_profile(0) ### Has to be commented out in artiq7 PC.
        # self.urukul1_cpld.set_profile(0)

        # URUKUL 0 - MOT and D2 state prep AOMs:
        delay(1*ms)
        self.urukul0_ch0.set(frequency=self.f_FORT, amplitude=self.AOM1_ampl)
        if self.FORT_AOM_ON == True:
            self.urukul0_ch0.sw.on()
        else:
            self.urukul0_ch0.sw.off()

        delay(1 * ms)
        self.urukul0_ch1.set(frequency=self.f_cooling_DP_MOT, amplitude=self.AOM2_ampl)
        if self.Cooling_DP_AOM_ON == True:
            self.urukul0_ch1.sw.on()
        else:
            self.urukul0_ch1.sw.off()

        delay(1 * ms)
        self.urukul0_ch2.set(frequency=self.f_D1_pumping_SP, amplitude=self.AOM3_ampl)
        if self.D1_pumping_SP_AOM_ON == True:
            self.urukul0_ch2.sw.on()
        else:
            self.urukul0_ch2.sw.off()

        delay(1 * ms)
        self.urukul0_ch3.set(frequency=self.f_pumping_RP, amplitude=self.AOM4_ampl)
        if self.pumping_RP_AOM_ON == True:
            self.urukul0_ch3.sw.on()
        else:
            self.urukul0_ch3.sw.off()

        # URUKUL 1 - MOT arm fiber AOMs:
        delay(1 * ms)
        self.urukul1_ch0.set(frequency=self.AOM_A2_freq, amplitude=self.AOM_A2_ampl)
        if self.AOM_A2_ON == True:
            self.urukul1_ch0.sw.on()
        else:
            self.urukul1_ch0.sw.off()

        delay(1 * ms)
        self.urukul1_ch1.set(frequency=self.AOM_A3_freq, amplitude=self.AOM_A3_ampl)
        if self.AOM_A3_ON == True:
            self.urukul1_ch1.sw.on()
        else:
            self.urukul1_ch1.sw.off()

        delay(1 * ms)
        self.urukul1_ch2.set(frequency=self.AOM_A6_freq, amplitude=self.AOM_A6_ampl)
        if self.AOM_A6_ON == True:
            self.urukul1_ch2.sw.on()
        else:
            self.urukul1_ch2.sw.off()

        delay(1 * ms)
        self.urukul1_ch3.set(frequency=self.AOM_A5_freq, amplitude=self.AOM_A5_ampl)
        if self.AOM_A5_ON == True:
            self.urukul1_ch3.sw.on()
        else:
            self.urukul1_ch3.sw.off()


        delay(10*ms)
        print("AOMs done!")