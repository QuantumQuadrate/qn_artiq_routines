"""
Used for rough measurement of the MOT loading time. It turns on the AOMs, wait for some time, then triggers
the Thorlabs camera with Zotino (the TTL was not triggering the camera, perhaps because it was too fast).

Could be extended to measure the MOT loading rate.
"""

from artiq.experiment import *
import math
import numpy as np
from ExperimentVariables import setattr_variables

class MOT_Load_Time(EnvExperiment):
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
        self.setattr_device("zotino0")  # for controlling coils
        self.setattr_device("ttl7")  # output for experiment trigger
        self.setattr_device("ttl6")
        self.setattr_device("ttl1")

        # import variables by name from datasets created by ExperimentVariables
        self.variables = [
            "f_FORT", "p_FORT_loading",
            "f_cooling_DP_MOT", "p_cooling_DP_MOT",
            "f_cooling_SP", "p_cooling_SP",
            "f_MOT_RP", "p_MOT_RP",
            "AOM_A1_freq", "AOM_A1_power",
            "AOM_A2_freq", "AOM_A2_power",
            "AOM_A3_freq", "AOM_A3_power",
            "AOM_A4_freq", "AOM_A4_power",
            "AOM_A5_freq", "AOM_A5_power",
            "AOM_A6_freq", "AOM_A6_power",
            "AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT",
            "cooling_setpoint_mW"
        ]

        # this adds the variables above as attributes in this experiment and gets their values.
        setattr_variables(self)

        self.setattr_argument("FORT_AOM_ON", BooleanValue(default=False), "AOM1, FORT")
        self.setattr_argument("Cooling_DP_AOM_ON", BooleanValue(default=False), "AOM2, MOT cooling double pass")
        self.setattr_argument("Cooling_SP_AOM_ON", BooleanValue(default=False), "AOM3, MOT cooling single pass")
        self.setattr_argument("MOT_RP_AOM_ON", BooleanValue(default=False), "AOM4, MOT RP/Exc")

        self.setattr_argument("AOM_A1_ON", BooleanValue(default=False), "AOM A2")
        self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "AOM A2")
        self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "AOM A3")
        self.setattr_argument("AOM_A4_ON", BooleanValue(default=False), "AOM A2")
        self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "AOM A5")
        self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "AOM A6")

    def prepare(self):
        # converts RF power in dBm to amplitudes in V
        self.AOM1_ampl = math.sqrt(2*50*10**(self.p_FORT_loading/10-3))
        self.AOM2_ampl = math.sqrt(2*50*10**(self.p_cooling_DP_MOT/10-3))
        self.AOM3_ampl = math.sqrt(2*50*10**(self.p_cooling_SP/10-3))
        self.AOM4_ampl = math.sqrt(2*50*10**(self.p_MOT_RP/10-3))

        self.AOM_A1_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A3_power / 10 - 3))
        self.AOM_A2_ampl = math.sqrt(2*50*10**(self.AOM_A2_power/10-3))
        self.AOM_A3_ampl = math.sqrt(2*50*10**(self.AOM_A3_power/10-3))
        self.AOM_A4_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A3_power / 10 - 3))
        self.AOM_A5_ampl = math.sqrt(2*50*10**(self.AOM_A5_power/10-3))
        self.AOM_A6_ampl = math.sqrt(2*50*10**(self.AOM_A6_power/10-3))

        self.coil_channels = [0, 1, 2, 3]


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

        # self.urukul0_cpld.set_profile(0) #### We had this in Ubuntu/Garret system. We have to comment out in new Artiq7 PC.
        # self.urukul1_cpld.set_profile(0)

        self.zotino0.init()

        self.ttl6.output()  # for outputting a trigger
        self.ttl7.output()  # for outputting a trigger
        self.ttl1.input()

        # delay(1000 * ms)

        # self.ttl6.pulse(100 * ms)

        ### URUKUL 0 - MOT and D2 state prep AOMs:
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
        self.urukul0_ch2.set(frequency=self.f_cooling_SP, amplitude=self.AOM3_ampl)
        if self.Cooling_SP_AOM_ON == True:
            self.urukul0_ch2.sw.on()
        else:
            self.urukul0_ch2.sw.off()

        delay(1 * ms)
        self.urukul0_ch3.set(frequency=self.f_MOT_RP, amplitude=self.AOM4_ampl)
        if self.MOT_RP_AOM_ON == True:
            self.urukul0_ch3.sw.on()
        else:
            self.urukul0_ch3.sw.off()

        ### URUKUL 1 - MOT arm fiber AOMs:
        # todo: channels below are not correct, but this initialization method is
        # about to be deprecated anyway. - PH, 2023.05.23
        delay(1 * ms)
        self.urukul1_ch1.set(frequency=self.AOM_A1_freq, amplitude=self.AOM_A1_ampl)
        if self.AOM_A1_ON == True:
            self.urukul1_ch0.sw.on()
        else:
            self.urukul1_ch0.sw.off()

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
        self.urukul1_ch1.set(frequency=self.AOM_A4_freq, amplitude=self.AOM_A4_ampl)
        if self.AOM_A4_ON == True:
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

        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                         channels=self.coil_channels)


        ### this is necessary to have the triggering signal after a certain delay. Otherwise, we do not get a trig signal.
        time1 = now_mu()
        tdelay = 2000 * ms
        tdelay_mu = self.core.seconds_to_mu(tdelay)
        delay(tdelay) # moves the cursor into the future
        self.core.wait_until_mu(time1 + tdelay_mu) # wait for the cursor to get there
        delay(1 * ms)

        ### trigger for the Thorlabs camera. a ttl pulse seems to rise too quickly
        self.zotino0.write_dac(6, 4.0)
        self.zotino0.load()
        delay(5*ms)
        self.zotino0.write_dac(6, 0.0)
        self.zotino0.load()

        ### time to wait for camera to take the image before shutting off AOMs
        time2 = now_mu()
        tdelay = self.core.seconds_to_mu(500 * ms)
        delay(500 * ms)
        self.core.wait_until_mu(time2 + tdelay)
        delay(1 * ms)

        ### turn off fiber AOMs, but leave the upstream AOMs on for thermal stability
        self.urukul1_ch0.sw.off()
        self.urukul1_ch1.sw.off()
        self.urukul1_ch2.sw.off()
        self.urukul1_ch3.sw.off()

        delay(1*ms)
        print("Coils and AOMs done!")