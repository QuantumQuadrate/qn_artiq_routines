"""
First turn on the AOMs with the AOMs-coils code, then run this code to change the settings of the coils and
see SPCM counts to put the MOT at the right location.
"""
from artiq.experiment import *

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")
from ExperimentVariables import setattr_variables


class CoilsSPCMCounts(EnvExperiment):
    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl0")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("zotino0")  # for controlling coils
        self.setattr_device("sampler0") # for measuring laser power PD

        self.setattr_argument("n_steps",
                              NumberValue(10, type='int', ndecimals=0, scale=1, step=1))
        # exposure time of the SPCM
        self.setattr_argument("dt_exposure", NumberValue(
            300 * ms))
        self.setattr_argument("print_count_rate", BooleanValue(True))

        # import variables by name from datasets created by ExperimentVariables
        self.variables = [
            "AZ_bottom_volts_MOT",
            "AZ_top_volts_MOT",
            "AX_volts_MOT",
            "AY_volts_MOT",
            "cooling_setpoint_mW"
        ]

        # this adds the variables above as attributes in this experiment and gets their values.
        setattr_variables(self)

        self.setattr_argument("disable_coils", BooleanValue(default=False))
        self.setattr_argument("enable_laser_feedback", BooleanValue(default=True),"Laser power stabilization")


    def prepare(self):

        self.coil_channels = [0, 1, 2, 3]

    @kernel
    def run(self):
        # initializes the hardware and resets attenuators to 0 dB
        self.core.reset()
        self.ttl0.input()

        self.core.break_realtime()

        self.zotino0.init()
        self.sampler0.init()


        if self.disable_coils:
            self.zotino0.set_dac([0.0,0.0,0.0,0.0],
                             channels=self.coil_channels)
        else:
            self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()

        delay(1000 * ms)

        for x in range(self.n_steps):
            tend1 = self.ttl0.gate_rising(self.dt_exposure)
            count1 = self.ttl0.count(tend1)
            if self.print_count_rate:
                print(round(count1/self.dt_exposure))
            delay(10 * ms)

        print("Coils done!")