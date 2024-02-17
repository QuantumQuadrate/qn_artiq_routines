"""
Used for rough measurement of the MOT loading time. It turns on the AOMs, wait for some time, then triggers
the Thorlabs camera with Zotino (the TTL was not triggering the camera, perhaps because it was too fast).

Could be extended to measure the MOT loading rate.
"""

from artiq.experiment import *

import sys
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\qn_artiq_routines\\')
sys.path.append('C:\\Users\\jakeuribe\\artiq-master\\repository\\')


from utilities.BaseExperiment import BaseExperiment


class MOT_Load_Time(EnvExperiment):
    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("FORT_AOM_ON", BooleanValue(default=False), "AOM1, FORT")
        self.setattr_argument("Cooling_DP_AOM_ON", BooleanValue(default=False), "AOM2, MOT cooling double pass")

        self.setattr_argument("AOM_A1_ON", BooleanValue(default=False), "AOM A2")
        self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "AOM A2")
        self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "AOM A3")
        self.setattr_argument("AOM_A4_ON", BooleanValue(default=False), "AOM A2")
        self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "AOM A5")
        self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "AOM A6")

        self.base.set_datasets_from_gui_args()


    def prepare(self):
        self.base.prepare()

    @kernel
    def run(self):
        self.base.initialize_hardware()

        # # URUKUL 0 - MOT and D2 state prep AOMs:
        delay(1 * ms)
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
        self.dds_AOM_A1.sw.off()
        self.dds_AOM_A2.sw.off()
        self.dds_AOM_A3.sw.off()
        self.dds_AOM_A4.sw.off()
        self.dds_AOM_A5.sw.off()
        self.dds_AOM_A6.sw.off()

        delay(1*ms)
        print("MOT load delay done!")