"""
This code turns on the MOT AOMs and also the MOT coils.
"""
from artiq.experiment import *
import math
import numpy as np
from subroutines.stabilizer import AOMPowerStabilizer

from BaseExperiment import BaseExperiment

class AOMsCoils(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

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
        self.base.prepare()

    @kernel
    def run(self):
        self.base.initialize_hardware()

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