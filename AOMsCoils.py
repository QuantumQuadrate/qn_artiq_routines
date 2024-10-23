"""
This code turns on the MOT AOMs and also the MOT coils.
"""
from artiq.experiment import *

from utilities.BaseExperiment import BaseExperiment

class AOMsCoils(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # experiment variables which are specific to this experiment
        self.setattr_argument("disable_coils", BooleanValue(default=False))
        self.setattr_argument("FORT_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("Cooling_DP_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("Repump_AOM_ON", BooleanValue(default=True))
        self.setattr_argument("D1_pumping_DP_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("pumping_repump_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("excitation_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("AOM_A1_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A4_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("microwave_dds_ON", BooleanValue(default=False), "Microwaves")
        self.setattr_argument("yes_Im_sure_I_want_the_microwave_dds_ON", BooleanValue(default=False), "Microwaves")
        self.setattr_argument("enable_laser_feedback_loop", BooleanValue(default=True),"Laser power stabilization")
        self.setattr_argument("run_laser_feedback_once", BooleanValue(default=False),"Laser power stabilization")
        self.setattr_argument("t_feedback_period", NumberValue(5*s, unit='s', ndecimals=1, step=1),
                              "Laser power stabilization")

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.base.prepare()

    @kernel
    def turn_on_AOMs(self):
        """
        turns on the AOMs that we have elected to turn on

        if the laser stabilizer is run, this should be called after that to turn on AOMs that were
        shut off for the feedback phase
        :return:
        """

        delay(1 * ms)
        if self.FORT_AOM_ON == True:
            self.dds_FORT.sw.on()
        else:
            self.dds_FORT.sw.off()

        delay(1 * ms)
        if self.Cooling_DP_AOM_ON == True:
            self.dds_cooling_DP.sw.on()
        else:
            self.dds_cooling_DP.sw.off()

        delay(1 * ms)
        if self.Repump_AOM_ON == True:
            self.ttl_repump_switch.off()
        else:
            self.ttl_repump_switch.on()

        delay(1 * ms)
        if self.D1_pumping_DP_AOM_ON == True:
            self.dds_D1_pumping_DP.sw.on()
        else:
            self.dds_D1_pumping_DP.sw.off()

        delay(1 * ms)
        if self.pumping_repump_AOM_ON == True:
            self.dds_pumping_repump.sw.on()
        else:
            self.dds_pumping_repump.sw.off()

        delay(1 * ms)
        if self.excitation_AOM_ON == True:
            self.dds_excitation.sw.on()
        else:
            self.dds_excitation.sw.off()

        # MOT arm fiber AOMs, excitation AOM:
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

        delay(1*ms)
        if self.microwave_dds_ON and self.yes_Im_sure_I_want_the_microwave_dds_ON:
            self.dds_microwaves.sw.on()
            self.ttl_microwave_switch.off()
        else:
            self.dds_microwaves.sw.off()
            self.ttl_microwave_switch.on()
        delay(1*ms)

    @kernel
    def run_feedback(self):
        self.core.reset()

        if self.AOM_A1_ON and self.AOM_A2_ON and self.AOM_A3_ON and self.AOM_A4_ON and self.AOM_A5_ON and self.AOM_A6_ON and self.Cooling_DP_AOM_ON:
            if self.enable_laser_feedback_loop:
                print("Will now run feedback and monitor powers until forcibly stopped")
                delay(100 * ms)

                while True:
                    self.laser_stabilizer.run()
                    delay(1 * ms)
                    self.turn_on_AOMs()
                    delay(1 * ms)
                    delay(self.t_feedback_period)

            elif self.run_laser_feedback_once:
                self.laser_stabilizer.run()  # this tunes the MOT and FORT AOMs
                delay(1 * ms)
                self.turn_on_AOMs()
                delay(1 * ms)
            else:
                # posts one data point for each beam
                delay(10 * ms)
                self.laser_stabilizer.monitor()
                delay(1 * ms)
                self.turn_on_AOMs()

    @kernel
    def aoms_and_coils(self):
        self.base.initialize_hardware()
        self.turn_on_AOMs()

        if self.disable_coils:
            self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],
                                 channels=self.coil_channels)
        else:
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                channels=self.coil_channels)

    def run(self):
        self.aoms_and_coils()
        self.run_feedback()
