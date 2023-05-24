"""
This code turns on the MOT AOMs and also the MOT coils.
"""
from artiq.experiment import *
import math
import numpy as np
from subroutines.stabilizer import AOMPowerStabilizer
from ExperimentVariables import setattr_variables
from DeviceAliases import DeviceAliases

class AOMsCoils(EnvExperiment):

    def build(self):

        # import variables by name from datasets created by ExperimentVariables.
        # must come before we set the named devices
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

        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul2_cpld")

        # map the hardware channels. urukul support only for now.
        self.named_devices = DeviceAliases(
            experiment=self,
            device_aliases=[
                'dds_FORT',
                'dds_cooling_SP',
                'dds_cooling_DP',
                'dds_MOT_RP',
                *[f'dds_AOM_A{i+1}' for i in range(6)] # the fiber AOMs
            ]
        )
        print("setattr by device aliases")

        # self.setattr_device("urukul0_ch0")
        # self.setattr_device("dds_cooling_DP")
        # self.setattr_device("dds_cooling_SP")
        # self.setattr_device("dds_MOT_RP")
        # self.setattr_device("dds_AOM_A2")
        # self.setattr_device("dds_AOM_A3")
        # self.setattr_device("dds_AOM_A1")
        # self.setattr_device("dds_AOM_A6")
        # self.setattr_device("dds_AOM_A4")
        # self.setattr_device("dds_AOM_A5")
        self.setattr_device("zotino0")  # for controlling coils
        self.setattr_device("sampler0") # for measuring laser power PD
        self.setattr_device("ttl6")
        self.setattr_device("ttl1")

        # experiment variables which are specific to this experiment
        self.setattr_argument("FORT_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("Cooling_DP_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("Cooling_SP_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("MOT_RP_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("AOM_A1_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A4_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("disable_coils", BooleanValue(default=False))
        self.setattr_argument("enable_laser_feedback", BooleanValue(default=True),"Laser power stabilization")

    def prepare(self):

        self.coil_channels = [0, 1, 2, 3]

        print("prepare - done!")

        # todo: eventually read conversion functions such as this from a config file
        def volts_to_optical_mW(x: TFloat) -> TFloat:
            """
            the conversion of PD voltage to cooling light power at the switchyard MOT 1 path
            """
            x += 0.011  # this accounts for a mismatch between what the Sampler reads and what
            # the multimeter that I used for the fit reads
            return -0.195395 + 17.9214 * x

        self.AOMservo = AOMPowerStabilizer(experiment=self,
                                           dds_names=["dds_cooling_DP"],
                                           sampler_name="sampler0",
                                           sampler_channels=[7],
                                           transfer_functions=[volts_to_optical_mW],
                                           setpoints=[self.cooling_setpoint_mW],  # in mW
                                           proportionals=[0.07],
                                           iters=5,  # if > x you'll underflow the rtio counter
                                           t_meas_delay=20 * ms)

    @kernel
    def run(self):

        # initializes the hardware and resets dds attenuators to 0 dB
        self.named_devices.initialize()

        self.core.break_realtime()

        self.ttl6.output()  # for outputting a trigger
        self.ttl1.input()
        self.sampler0.init()

        # # URUKUL 0 - MOT and D2 state prep AOMs:
        delay(1*ms)
        if self.FORT_AOM_ON == True:
            self.dds_FORT.sw.on()
        else:
            self.dds_FORT.sw.off()
        #
        delay(1 * ms)
        if self.Cooling_DP_AOM_ON == True:
            self.dds_cooling_DP.sw.on()
        else:
            self.dds_cooling_DP.sw.off()

        delay(1 * ms)
        if self.Cooling_SP_AOM_ON == True:
            self.dds_cooling_SP.sw.on()
        else:
            self.dds_cooling_SP.sw.off()

        delay(1 * ms)
        if self.MOT_RP_AOM_ON == True:
            self.dds_MOT_RP.sw.on()
        else:
            self.dds_MOT_RP.sw.off()

        # URUKUL 1 and 2 - MOT arm fiber AOMs:
        delay(1 * ms)
        if self.AOM_A2_ON == True:
            self.dds_AOM_A2.sw.on()
        else:
            self.dds_AOM_A2.sw.off()

        delay(1 * ms)
        if self.AOM_A3_ON == True:
            self.dds_AOM_A3.sw.on()
        else:
            self.dds_AOM_A3.sw.off()

        delay(1 * ms)
        if self.AOM_A1_ON == True:
            self.dds_AOM_A1.sw.on()
        else:
            self.dds_AOM_A1.sw.off()

        delay(1 * ms)
        if self.AOM_A6_ON == True:
            self.dds_AOM_A6.sw.on()
        else:
            self.dds_AOM_A6.sw.off()

        delay(1 * ms)
        if self.AOM_A4_ON == True:
            self.dds_AOM_A4.sw.on()
        else:
            self.dds_AOM_A4.sw.off()

        delay(1 * ms)
        if self.AOM_A5_ON == True:
            self.dds_AOM_A5.sw.on()
        else:
            self.dds_AOM_A5.sw.off()

        if self.disable_coils:
            self.zotino0.set_dac([0.0,0.0,0.0,0.0],
                             channels=self.coil_channels)
        else:
            self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)

        if self.enable_laser_feedback:
            self.AOMservo.get_dds_settings()  # must come after relevant DDS's have been set
            print("waiting for AOMs to thermalize")
            delay(2000 * ms)
            print("running feedback")
            self.AOMservo.run()

            delay(1*ms)
        print("Coils and AOMs done!")