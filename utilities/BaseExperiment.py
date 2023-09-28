"""
Setup for a "base experiment".

To allow us to be flexible with the type of ARTIQ we run in the future,
we do not inherit from an artiq experiment, but add variables and devices
as attributes to an experiment passed by reference to the functions below.

There are also some quality of life functions here.

intended usage:
----
from BaseExperiment import BaseExperiment


class MyExp(EnvExperiment):

    def build:
        # gets variables and devices
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # define GUI arguments here
        self.setattr_argument("my_argument",NumberValue(5))

        self.base.set_datasets_from_gui_args()

    def prepare:
        self.base.prepare()

    @kernel
    def run:
        self.base.initialize_hardware()
        # now we can do physics

        # reference devices and variables which were initialized in the base experiment
        self.dds_AOM_A1.sw.on()
        self.ttl6.pulse(self.t_pulse)
----
"""
from artiq.experiment import *

import sys, os
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')

from subroutines.aom_feedback import AOMPowerStabilizer
from ExperimentVariables import setattr_variables
from utilities.DeviceAliases import DeviceAliases


def dB_to_V(dB: float) -> float:
    """
    convert power in dB to volts for setting DDS amplitudes
    :return amplitude: float in volts
    """
    return (2 * 50 * 10 ** (dB / 10 - 3)) ** (1 / 2)


class BaseExperiment:

    def __init__(self, experiment: EnvExperiment):
        self.experiment = experiment

    def build(self):
        """
        Put this in your experiment's build method with experiment=self

        :param experiment: your experiment.
        :return:
        """

        # todo: find a way to extract the names of all of the variables
        #  in ExperimentVariables.py, e.g. by reading them in from the
        #  most recent hdf file?
        self.experiment.variables = [
            "f_FORT", "p_FORT_loading", "p_FORT_RO", "p_FORT_PGC",
            "f_cooling_DP_MOT", "p_cooling_DP_MOT",
            "f_cooling_DP_PGC", "p_cooling_DP_PGC",
            "f_cooling_DP_RO", "p_cooling_DP_RO",
            "f_cooling_SP", "p_cooling_SP",
            "f_MOT_RP", "p_MOT_RP",
            "AOM_A1_freq", "AOM_A1_power",
            "AOM_A2_freq", "AOM_A2_power",
            "AOM_A3_freq", "AOM_A3_power",
            "AOM_A4_freq", "AOM_A4_power",
            "AOM_A5_freq", "AOM_A5_power",
            "AOM_A6_freq", "AOM_A6_power",
            "AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT",
            "AZ_bottom_volts_PGC", "AZ_top_volts_PGC", "AX_volts_PGC", "AY_volts_PGC",
            "AZ_bottom_volts_RO", "AZ_top_volts_RO", "AX_volts_RO", "AY_volts_RO",
            "Luca_trigger_for_feedback_verification"
            "enable_laser_feedback",
            "cooling_setpoint_mW",
            "cooling_volts_ch",
            "t_MOT_loading",
            "t_FORT_loading",
            "t_SPCM_exposure",
            "set_point_PD0_AOM_cooling_DP",
            "set_point_fW_AOM_A1",
            "set_point_fW_AOM_A2",
            "set_point_fW_AOM_A3",
            "set_point_fW_AOM_A4",
            "set_point_PD5_AOM_A5",
            "set_point_PD6_AOM_A6",
            'set_point_FORT_MM',
            'MOT_beam_monitor_points',
            'feedback_dds_list',
            'aom_feedback_averages',
            'aom_feedback_iterations'
        ]

        setattr_variables(self.experiment)

        # devices without nicknames. core should come first
        devices_no_alias = ["core",
                            "urukul0_cpld", "urukul1_cpld", "urukul2_cpld",
                            "zotino0",  # for controlling coils
                            "sampler0",  # for measuring laser power PD
                            "sampler1", # for reading in volts in the coil tune experiment
                            *[f"ttl{i}" for i in range(8)]]
        for dev in devices_no_alias:
            # print(f"setting {dev}")
            self.experiment.setattr_device(dev)

        # initialize named channels.
        self.experiment.named_devices = DeviceAliases(
            experiment=self.experiment,
            device_aliases=[
                'dds_FORT',
                'dds_cooling_SP',
                'dds_cooling_DP',
                'dds_MOT_RP',
                *[f'dds_AOM_A{i + 1}' for i in range(6)]  # the fiber AOMs
            ]
        )

        self.experiment.coil_channels = [0, 1, 2, 3]

        # get a list of all attributes of experiment up to this point. if base.build is called in your experiment
        # before any GUI arguments are defined, then this can be used to grab those later by taking a difference
        self.exp_var_names = dir(self.experiment)
        print("base build - done")


    def set_datasets_from_gui_args(self):
        new_exp_var_names = [x for x in dir(self.experiment) if x not in self.exp_var_names]
        for name in new_exp_var_names:
            try:
                self.experiment.set_dataset(name, getattr(self.experiment, name))
            except Exception as e:
                print("ARTIQ complains about this when scanning repository HEAD but then gets over it...")


    def prepare(self):
        """
        Compute DDS amplitudes from powers, instantiate the laser servo, any other math
        that needs to happen before we run stuff on the kernel.
        :return:
        """

        # convert times to machine units
        seconds_to_mu = self.experiment.core.seconds_to_mu
        self.experiment.t_MOT_loading_mu = seconds_to_mu(self.experiment.t_MOT_loading)
        self.experiment.t_FORT_loading_mu = seconds_to_mu(self.experiment.t_FORT_loading)
        self.experiment.t_SPCM_exposure_mu = seconds_to_mu(self.experiment.t_SPCM_exposure)

        # converts RF power in dBm to amplitudes in V
        self.experiment.ampl_FORT_loading = dB_to_V(self.experiment.p_FORT_loading)
        self.experiment.ampl_FORT_RO = dB_to_V(self.experiment.p_FORT_RO)
        self.experiment.ampl_FORT_PGC = dB_to_V(self.experiment.p_FORT_PGC)
        self.experiment.ampl_cooling_DP_MOT = dB_to_V(self.experiment.p_cooling_DP_MOT)
        self.experiment.ampl_cooling_DP_PGC =  dB_to_V(self.experiment.p_cooling_DP_PGC)
        self.experiment.ampl_cooling_DP_RO = dB_to_V(self.experiment.p_cooling_DP_RO)
        self.experiment.ampl_cooling_SP = dB_to_V(self.experiment.p_cooling_SP)
        self.experiment.ampl_MOT_RP = dB_to_V(self.experiment.p_MOT_RP)

        self.experiment.AOM_A1_ampl = dB_to_V(self.experiment.AOM_A1_power)
        self.experiment.AOM_A2_ampl = dB_to_V(self.experiment.AOM_A2_power)
        self.experiment.AOM_A3_ampl = dB_to_V(self.experiment.AOM_A3_power)
        self.experiment.AOM_A4_ampl = dB_to_V(self.experiment.AOM_A4_power)
        self.experiment.AOM_A5_ampl = dB_to_V(self.experiment.AOM_A5_power)
        self.experiment.AOM_A6_ampl = dB_to_V(self.experiment.AOM_A6_power)

        dds_feedback_list = eval(self.experiment.feedback_dds_list)

        self.experiment.laser_stabilizer = AOMPowerStabilizer(experiment=self.experiment,
                                                    dds_names=dds_feedback_list,
                                                    iterations=self.experiment.aom_feedback_iterations,
                                                    averages=self.experiment.aom_feedback_averages,
                                                    leave_AOMs_on=True)

        print("base prepare - done")

    @kernel
    def initialize_hardware(self):
        """
        what it sounds like
        :return:
        """
        self.experiment.named_devices.initialize()

        self.experiment.ttl6.output()  # for outputting a trigger
        self.experiment.ttl1.input()
        self.experiment.sampler0.init() # for reading laser feedback

        print("base initialize_hardware - done")

        self.experiment.core.break_realtime()


# do this so the code above will not actually run when ARTIQ scans the repository
if __name__ == '__main__':
    pass


