"""
This code turns on the MOT AOMs and also the MOT coils.
"""
from artiq.experiment import *
import math
import numpy as np


class AOMsCoils(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")
        self.setattr_device("urukul2_ch0")
        self.setattr_device("urukul2_ch1")
        self.setattr_device("zotino0")  # for controlling coils
        self.setattr_device("ttl6")
        self.setattr_device("ttl1")


        # self.setattr_argument("AOM1_freq", NumberValue(345.1 * MHz, unit="MHz", ndecimals=1), "AOM1, pumping repumper")
        # self.setattr_argument("AOM1_power", NumberValue(-5, unit="dBm", scale=1, ndecimals=1), "AOM1, pumping repumper")
        # self.setattr_argument("AOM1_ON", BooleanValue(default=False), "AOM1, pumping repumper")

        self.setattr_argument("AOM1_freq", NumberValue(210 * MHz, unit="MHz", ndecimals=1), "AOM1, FORT")
        self.setattr_argument("AOM1_power", NumberValue(3, unit="dBm", scale=1, ndecimals=1), "AOM1, FORT")
        self.setattr_argument("AOM1_ON", BooleanValue(default=False), "AOM1, FORT")

        # -0.2 gives the maximum diffraction efficiency. -4 gives about 70% so we can increase to compensate for drift
        self.setattr_argument("AOM2_freq", NumberValue(111.0 * MHz, unit="MHz",ndecimals=1), "AOM2, MOT cooling double pass")
        self.setattr_argument("AOM2_power", NumberValue(-4, unit="dBm", scale=1, ndecimals=1), "AOM2, MOT cooling double pass")
        self.setattr_argument("AOM2_ON", BooleanValue(default=False), "AOM2, MOT cooling double pass")

        self.setattr_argument("AOM3_freq", NumberValue(130.0 * MHz, unit="MHz",ndecimals=1), "AOM3, MOT cooling single pass")
        self.setattr_argument("AOM3_power", NumberValue(1, unit="dBm", scale=1, ndecimals=1), "AOM3, MOT cooling single pass")
        self.setattr_argument("AOM3_ON", BooleanValue(default=False), "AOM3, MOT cooling single pass")

        self.setattr_argument("AOM4_freq", NumberValue(150.5 * MHz, unit="MHz", ndecimals=1), "AOM4, MOT RP/Exc")
        self.setattr_argument("AOM4_power", NumberValue(3, unit="dBm", scale=1, ndecimals=1), "AOM4, MOT RP/Exc")
        self.setattr_argument("AOM4_ON", BooleanValue(default=False), "AOM4, MOT RP/Exc")
        # the default power for the fiber AOMs was chosen to give roughly equal diffraction efficiency, empirically

        self.setattr_argument("AOM_A1_freq", NumberValue(78.51 * MHz, unit="MHz", ndecimals=2), "AOM A1")
        self.setattr_argument("AOM_A1_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A1")
        self.setattr_argument("AOM_A1_ON", BooleanValue(default=False), "AOM A1")

        self.setattr_argument("AOM_A2_freq", NumberValue(78.48 * MHz, unit="MHz", ndecimals=2), "AOM A2")
        self.setattr_argument("AOM_A2_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A2")
        self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "AOM A2")

        self.setattr_argument("AOM_A3_freq", NumberValue(78.49 * MHz, unit="MHz", ndecimals=2), "AOM A3")
        self.setattr_argument("AOM_A3_power", NumberValue(-3, unit="dBm", scale=1, ndecimals=1), "AOM A3")
        self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "AOM A3")

        self.setattr_argument("AOM_A4_freq", NumberValue(78.5 * MHz, unit="MHz", ndecimals=2), "AOM A4")
        self.setattr_argument("AOM_A4_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A4")
        self.setattr_argument("AOM_A4_ON", BooleanValue(default=False), "AOM A4")

        self.setattr_argument("AOM_A5_freq", NumberValue(78.47 * MHz, unit="MHz", ndecimals=2), "AOM A5")
        self.setattr_argument("AOM_A5_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A5")
        self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "AOM A5")

        self.setattr_argument("AOM_A6_freq", NumberValue(78.52 * MHz, unit="MHz", ndecimals=2), "AOM A6")
        self.setattr_argument("AOM_A6_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A6")
        self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "AOM A6")

        self.setattr_argument("AZ_bottom_volts_MOT", NumberValue(0.6, unit="V", ndecimals=3, step=0.01),
                              "A-Z shim/quad bottom coils")
        self.setattr_argument("AZ_top_volts_MOT", NumberValue(-2.5, unit="V", ndecimals=3, step=0.01),
                              "A-Z shim/quad top coil")
        self.setattr_argument("AX_volts_MOT", NumberValue(-0.19, unit="V", ndecimals=3, step=0.01), "A-X shim coils")
        self.setattr_argument("AY_volts_MOT", NumberValue(-0.25, unit="V", ndecimals=3, step=0.01), "A-Y shim coils")
        self.setattr_argument("disable_coils", BooleanValue(default=False))


    def prepare(self):
        # converts RF power in dBm to amplitudes in V
        self.AOM1_ampl = math.sqrt(2*50*10**(self.AOM1_power/10-3))
        self.AOM2_ampl = math.sqrt(2*50*10**(self.AOM2_power/10-3))
        self.AOM3_ampl = math.sqrt(2*50*10**(self.AOM3_power/10-3))
        self.AOM4_ampl = math.sqrt(2*50*10**(self.AOM4_power/10-3))

        self.AOM_A1_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A1_power / 10 - 3))
        self.AOM_A2_ampl = math.sqrt(2*50*10**(self.AOM_A2_power/10-3))
        self.AOM_A3_ampl = math.sqrt(2*50*10**(self.AOM_A3_power/10-3))
        self.AOM_A4_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A4_power / 10 - 3))
        self.AOM_A5_ampl = math.sqrt(2*50*10**(self.AOM_A5_power/10-3))
        self.AOM_A6_ampl = math.sqrt(2*50*10**(self.AOM_A6_power/10-3))


        self.coil_channels = [0, 1, 2, 3]


    @kernel
    def run(self):
    #     #initializes the hardware and resets attenuators to 0 dB
        self.core.reset()
        self.urukul0_cpld.init()
        self.urukul1_cpld.init()
        self.urukul2_cpld.init()

        self.urukul0_ch0.init()
        self.urukul0_ch1.init()
        self.urukul0_ch2.init()
        self.urukul0_ch3.init()
        self.urukul1_ch0.init()
        self.urukul1_ch1.init()
        self.urukul1_ch2.init()
        self.urukul1_ch3.init()
        self.urukul2_ch0.init()
        self.urukul2_ch1.init()

        self.urukul0_ch0.set_att(float(0))
        self.urukul0_ch1.set_att(float(0))
        self.urukul0_ch2.set_att(float(0))
        self.urukul0_ch3.set_att(float(0))
        self.urukul1_ch0.set_att(float(0))
        self.urukul1_ch1.set_att(float(0))
        self.urukul1_ch2.set_att(float(0))
        self.urukul1_ch3.set_att(float(0))

        self.core.break_realtime()

        # self.urukul0_cpld.set_profile(0)  #### We had this in Ubuntu/Garret system. We have to comment out in new Artiq7 PC.
        # self.urukul1_cpld.set_profile(0)

        self.zotino0.init()

        self.ttl6.output()  # for outputting a trigger
        self.ttl1.input()

        # delay(1000 * ms)

        # self.ttl6.pulse(100 * ms)

        # URUKUL 0 - MOT and D2 state prep AOMs:
        delay(1*ms)
        self.urukul0_ch0.set(frequency=self.AOM1_freq, amplitude=self.AOM1_ampl)
        if self.AOM1_ON == True:
            self.urukul0_ch0.sw.on()
        else:
            self.urukul0_ch0.sw.off()

        delay(1 * ms)
        self.urukul0_ch1.set(frequency=self.AOM2_freq, amplitude=self.AOM2_ampl)
        if self.AOM2_ON == True:
            self.urukul0_ch1.sw.on()
        else:
            self.urukul0_ch1.sw.off()

        delay(1 * ms)
        self.urukul0_ch2.set(frequency=self.AOM3_freq, amplitude=self.AOM3_ampl)
        if self.AOM3_ON == True:
            self.urukul0_ch2.sw.on()
        else:
            self.urukul0_ch2.sw.off()

        delay(1 * ms)
        self.urukul0_ch3.set(frequency=self.AOM4_freq, amplitude=self.AOM4_ampl)
        if self.AOM4_ON == True:
            self.urukul0_ch3.sw.on()
        else:
            self.urukul0_ch3.sw.off()

        # URUKUL 1 and 2 - MOT arm fiber AOMs:
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
        self.urukul1_ch2.set(frequency=self.AOM_A1_freq, amplitude=self.AOM_A1_ampl)
        if self.AOM_A1_ON == True:
            self.urukul1_ch2.sw.on()
        else:
            self.urukul1_ch2.sw.off()

        delay(1 * ms)
        self.urukul1_ch3.set(frequency=self.AOM_A6_freq, amplitude=self.AOM_A6_ampl)
        if self.AOM_A6_ON == True:
            self.urukul1_ch3.sw.on()
        else:
            self.urukul1_ch3.sw.off()

        delay(1 * ms)
        self.urukul2_ch0.set(frequency=self.AOM_A4_freq, amplitude=self.AOM_A4_ampl)
        if self.AOM_A4_ON == True:
            self.urukul2_ch0.sw.on()
        else:
            self.urukul2_ch0.sw.off()

        delay(1 * ms)
        self.urukul2_ch1.set(frequency=self.AOM_A5_freq, amplitude=self.AOM_A5_ampl)
        if self.AOM_A5_ON == True:
            self.urukul2_ch1.sw.on()
        else:
            self.urukul2_ch1.sw.off()



        if self.disable_coils:
            self.zotino0.set_dac([0.0,0.0,0.0,0.0],
                             channels=self.coil_channels)
        else:
            self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)

        delay(1*ms)
        print("Coils and AOMs done!")