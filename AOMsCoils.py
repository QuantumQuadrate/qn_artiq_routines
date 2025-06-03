"""
This code turns on the MOT AOMs and also the MOT coils.
"""
from artiq.experiment import *

from utilities.BaseExperiment import BaseExperiment
from subroutines.k10cr1_functions import *

class AOMsCoils(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # experiment variables which are specific to this experiment
        self.setattr_argument("disable_coils", BooleanValue(default=False))
        self.setattr_argument("FORT_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("Cooling_DP_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("Repump_AOM_switch_ON", BooleanValue(default=True))
        self.setattr_argument("D1_pumping_DP_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("pumping_repump_AOM_ON", BooleanValue(default=False))
        self.setattr_argument("Excitation0_AOM_switch_ON", BooleanValue(default=False))
        self.setattr_argument("GRIN1and2_DDS_ON", BooleanValue(default=False))
        self.setattr_argument("GRIN1_AOM_switch_ON", BooleanValue(default=False))
        self.setattr_argument("GRIN2_AOM_switch_ON", BooleanValue(default=False))
        self.setattr_argument("AOM_A1_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A4_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "Fiber AOMs")
        self.setattr_argument("microwave_dds_ON", BooleanValue(default=False), "Microwaves")
        self.setattr_argument("yes_Im_sure_I_want_the_microwave_dds_ON", BooleanValue(default=False), "Microwaves")
        self.setattr_argument("run_laser_feedback", BooleanValue(default=False), "Laser power stabilization")

        self.setattr_argument("go_to_home_780HWP", BooleanValue(default=False), "K10CR1 780 waveplates")
        self.setattr_argument("go_to_home_780QWP", BooleanValue(default=False), "K10CR1 780 waveplates")
        self.setattr_argument("go_to_target_780HWP", BooleanValue(default=False), "K10CR1 780 waveplates")
        self.setattr_argument("go_to_target_780QWP", BooleanValue(default=False), "K10CR1 780 waveplates")
        self.setattr_argument("move_780HWP_by", BooleanValue(default=False), "K10CR1 780 waveplates")
        self.setattr_argument("move_780QWP_by", BooleanValue(default=False), "K10CR1 780 waveplates")

        self.setattr_argument("go_to_home_852HWP", BooleanValue(default=False), "K10CR1 852 waveplates")
        self.setattr_argument("go_to_home_852QWP", BooleanValue(default=False), "K10CR1 852 waveplates")

        self.setattr_argument("go_to_target_852HWP", BooleanValue(default=False), "K10CR1 852 waveplates")
        self.setattr_argument("go_to_target_852QWP", BooleanValue(default=False), "K10CR1 852 waveplates")
        self.setattr_argument("move_852HWP_by", BooleanValue(default=False), "K10CR1 852 waveplates")
        self.setattr_argument("move_852QWP_by", BooleanValue(default=False), "K10CR1 852 waveplates")

        self.setattr_argument("go_to_optimized_852_settings", BooleanValue(default=False), "K10CR1 852 waveplates")

        self.base.set_datasets_from_gui_args()

    def prepare(self):
        self.base.prepare()

        #todo: put this in BaseExperiment.py
        self.setpoint_datasets = ["best_HWP_to_H","best_QWP_to_H"]
        self.default_setpoints = [getattr(self, dataset) for dataset in self.setpoint_datasets]

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
        if self.Repump_AOM_switch_ON == True:
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
        if self.Excitation0_AOM_switch_ON == True:
            self.ttl_exc0_switch.off()
        else:
            self.ttl_exc0_switch.on()

        delay(1 * ms)
        if self.GRIN1and2_DDS_ON == True:
            self.GRIN1and2_dds.sw.on()
        else:
            self.GRIN1and2_dds.sw.off()

        delay(1 * ms)
        if self.GRIN1_AOM_switch_ON == True:
            self.ttl_GRIN1_switch.off()
        else:
            self.ttl_GRIN1_switch.on()

        delay(1 * ms)
        if self.GRIN2_AOM_switch_ON == True:
            self.ttl_GRIN2_switch.off()
        else:
            self.ttl_GRIN2_switch.on()


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

        # todo: if we redeclare the AOMPowerStabilizer instance in prepare with only the AOMs that we wish to turn on,
        #  we can get rid of this if statement. this will be important for when we have more AOMs that we want to
        #  feed back to, e.g. OP and excitation.
        if self.AOM_A1_ON and self.AOM_A2_ON and self.AOM_A3_ON and self.AOM_A4_ON and self.AOM_A5_ON and self.AOM_A6_ON and self.Cooling_DP_AOM_ON:

            if self.run_laser_feedback:
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

    @kernel
    def k10cr1_operations(self):
        self.core.reset()
        delay(10*ms)
        if self.go_to_home_780HWP:
            go_to_home(self, '780_HWP')
        if self.go_to_home_780QWP:
            go_to_home(self, '780_QWP')

        if self.go_to_target_780HWP:
            move_to_target_deg(self, name="780_HWP", target_deg=self.target_780_HWP)
        if self.go_to_target_780QWP:
            move_to_target_deg(self, name="780_QWP", target_deg=self.target_780_QWP)

        if self.move_780HWP_by:
            move_by_deg(self, name="780_HWP", target_deg=self.move_780_HWP_by)
        if self.move_780QWP_by:
            move_by_deg(self, name="780_QWP", target_deg=self.move_780_QWP_by)

        if self.go_to_home_852HWP:
            go_to_home(self, '852_HWP')
        if self.go_to_home_852QWP:
            go_to_home(self, '852_QWP')

        if self.go_to_target_780HWP:
            move_to_target_deg(self, name="852_HWP", target_deg=self.target_852_HWP)
        if self.go_to_target_780QWP:
            move_to_target_deg(self, name="852_QWP", target_deg=self.target_852_QWP)

        if self.move_852HWP_by:
            move_by_deg(self, name="852_HWP", target_deg=self.move_852_HWP_by)
        if self.move_852QWP_by:
            move_by_deg(self, name="852_QWP", target_deg=self.move_852_QWP_by)


        if self.go_to_optimized_852_settings:
            move_to_target_deg(self, name="852_HWP", target_deg=self.best_852HWP_to_max)
            move_to_target_deg(self, name="852_QWP", target_deg=self.best_852QWP_to_max)



    def run(self):
        self.aoms_and_coils()
        self.run_feedback()
        self.k10cr1_operations()
