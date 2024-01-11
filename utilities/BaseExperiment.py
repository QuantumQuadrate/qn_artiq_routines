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
import logging
import sys, os
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\')

from subroutines.aom_feedback import AOMPowerStabilizer
from ExperimentVariables import setattr_variables
from utilities.DeviceAliases import DeviceAliases
from utilities.write_h5 import write_results

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

        with open('C:\\Networking Experiment\\artiq codes\\artiq-master\\dataset_db.pyon') as f:
            datasets_str = f.read()

        # when the pyon file is saved python True and False are converted to lowercase...
        datasets_str = datasets_str.replace("true", "True")
        datasets_str = datasets_str.replace("false", "False")
        datasets_dict = eval(datasets_str)

        # these are the names of all of the datasets, which we're going to use to create attributes of the same name
        self.experiment.variables = datasets_dict.keys()

        setattr_variables(self.experiment)

        # devices without nicknames. core should come first
        devices_no_alias = ["core",
                            "urukul0_cpld", "urukul1_cpld", "urukul2_cpld",
                            "zotino0",  # for controlling coils
                            "sampler0",  # for measuring laser power PD
                            "sampler1", # for reading in volts in the coil tune experiment
                            "sampler2",
                            *[f"ttl{i}" for i in range(16)]]
        for dev in devices_no_alias:
            # print(f"setting {dev}")
            self.experiment.setattr_device(dev)

        # devices can also be nicknamed here:
        self.experiment.ttl_microwave_switch = self.experiment.ttl4
        self.experiment.ttl_microwave_switch = self.experiment.ttl4
        self.experiment.ttl_repump_switch = self.experiment.ttl5
        self.experiment.ttl_Luca_trigger = self.experiment.ttl6

        # initialize named channels.
        self.experiment.named_devices = DeviceAliases(
            experiment=self.experiment,
            device_aliases=[
                'dds_FORT',
                'dds_D1_pumping_SP',
                'dds_cooling_DP',
                'dds_pumping_repump',
                'dds_excitation',
                'dds_microwaves',
                *[f'dds_AOM_A{i + 1}' for i in range(6)]  # the fiber AOMs
            ]
        )

        # for debugging/logging purposes in experiments
        self.experiment.coil_names = ["AZ bottom","AZ top","AX","AY"]

        self.experiment.AZ_bottom_Zotino_channel = 0
        self.experiment.AZ_top_Zotino_channel = 1
        self.experiment.AX_Zotino_channel = 2
        self.experiment.AY_Zotino_channel = 3

        self.experiment.coil_channels = [self.experiment.AZ_bottom_Zotino_channel,
                                         self.experiment.AZ_top_Zotino_channel,
                                         self.experiment.AX_Zotino_channel,
                                         self.experiment.AY_Zotino_channel]

        # this is an attribute of of the experiment in case we want to access it elsewhere
        self.experiment.all_dds_channels = [getattr(self.experiment,f'urukul{card}_ch{channel}')
                                 for card in range(3) for channel in range(4)]

        # todo: should limit this list to the channels which are outputs
        # self.experiment.all_ttl_channels = [getattr(self.experiment, 'ttl{channel}')
        #                                     for channel in range(16)]

        # dataset names
        self.experiment.count_rate_dataset = 'photocounts_per_s'

        # functions

        @rpc(flags={"async"})
        def print_async(*x):
            """print asynchronously so we don't block the RTIO counter.
            useful for debugging"""
            print(*x)
        self.experiment.print_async = print_async
        self.experiment.write_results = lambda: write_results(experiment=self.experiment)

        # get a list of all attributes of experiment up to this point. if base.build is called in your experiment
        # before any GUI arguments are defined, then this can be used to grab those later by taking a difference
        self.exp_var_names = dir(self.experiment)
        logging.debug("base build - done")


    def set_datasets_from_gui_args(self):
        new_exp_var_names = [x for x in dir(self.experiment) if x not in self.exp_var_names]
        for name in new_exp_var_names:
            try:
                self.experiment.set_dataset(name, getattr(self.experiment, name))
            except Exception as e:
                pass # this is terrible but ARTIQ prints out way too many of these
                # print("ARTIQ complains about this when scanning repository HEAD but then gets over it...")


    def prepare(self):
        """
        Compute DDS amplitudes from powers, instantiate the laser servo, any other math
        that needs to happen before we run stuff on the kernel.
        :return:
        """

        # todo: because we only run feedback strictly <= once per experiment sequence,
        #  different amplitude settings could be expressed as fractions of the amplitudes we feed back to.
        #  e.g., the PGC amplitude for the cooling laser should be a fraction, say, 0.9, that we multiply the
        #  MOT power by. if we want these amplitudes to be able to be tuned completely independently, then we
        #  need to feedback to them individually, but the thermal equilibrium of the AOM is likely to be different
        #  for each setting, so tuning up e.g. the cooling amplitude, may be completely "undone" by subsequently
        #  finding the right amplitude for the PGC.

        # convert times to machine units
        seconds_to_mu = self.experiment.core.seconds_to_mu
        self.experiment.t_MOT_loading_mu = seconds_to_mu(self.experiment.t_MOT_loading)
        self.experiment.t_FORT_loading_mu = seconds_to_mu(self.experiment.t_FORT_loading)
        self.experiment.t_SPCM_exposure_mu = seconds_to_mu(self.experiment.t_SPCM_exposure)

        # converts RF power in dBm to amplitudes in V
        self.experiment.ampl_FORT_loading = dB_to_V(self.experiment.p_FORT_loading)
        self.experiment.ampl_cooling_DP_MOT = dB_to_V(self.experiment.p_cooling_DP_MOT)
        self.experiment.ampl_D1_pumping_SP = dB_to_V(self.experiment.p_D1_pumping_SP)
        self.experiment.ampl_pumping_repump = dB_to_V(self.experiment.p_pumping_repump)
        self.experiment.ampl_D1_pumping_SP = dB_to_V(self.experiment.p_D1_pumping_SP)
        self.experiment.ampl_excitation = dB_to_V(self.experiment.p_excitation)
        self.experiment.AOM_A1_ampl = dB_to_V(self.experiment.AOM_A1_power)
        self.experiment.AOM_A2_ampl = dB_to_V(self.experiment.AOM_A2_power)
        self.experiment.AOM_A3_ampl = dB_to_V(self.experiment.AOM_A3_power)
        self.experiment.AOM_A4_ampl = dB_to_V(self.experiment.AOM_A4_power)
        self.experiment.AOM_A5_ampl = dB_to_V(self.experiment.AOM_A5_power)
        self.experiment.AOM_A6_ampl = dB_to_V(self.experiment.AOM_A6_power)

        # RF powers defined as fractions of the defaults, e.g. the ones we tune during the AOM feedback
        self.experiment.ampl_FORT_RO = self.experiment.ampl_FORT_loading * self.experiment.p_FORT_RO
        self.experiment.ampl_FORT_PGC = self.experiment.ampl_FORT_loading * self.experiment.p_FORT_PGC
        self.experiment.ampl_FORT_blowaway = self.experiment.ampl_FORT_loading * self.experiment.p_FORT_blowaway
        self.experiment.ampl_FORT_OP = self.experiment.ampl_FORT_loading * self.experiment.p_FORT_OP
        self.experiment.ampl_cooling_DP_RO = self.experiment.ampl_cooling_DP_MOT * self.experiment.p_cooling_DP_RO
        self.experiment.ampl_cooling_DP_PGC = self.experiment.ampl_cooling_DP_MOT * self.experiment.p_cooling_DP_PGC

        dds_feedback_list = eval(self.experiment.feedback_dds_list)

        self.experiment.laser_stabilizer = AOMPowerStabilizer(experiment=self.experiment,
                                                    dds_names=dds_feedback_list,
                                                    iterations=self.experiment.aom_feedback_iterations,
                                                    averages=self.experiment.aom_feedback_averages,
                                                    leave_MOT_AOMs_on=True)

        logging.debug("base prepare - done")

    @kernel
    def initialize_hardware(self):
        """
        hardware initialization and setting of ttl switches, and set datasets
        :return:
        """
        self.experiment.named_devices.initialize()

        self.experiment.ttl_microwave_switch.output()
        self.experiment.ttl_repump_switch.output()
        self.experiment.ttl6.output()  # for outputting a trigger
        self.experiment.ttl1.input()
        self.experiment.sampler0.init() # for reading laser feedback
        self.experiment.sampler1.init() # for reading laser feedback
        self.experiment.sampler2.init() # for reading laser feedback

        self.experiment.print_async("base initialize_hardware - done")

        # turn on any switches. this ensures that switches always start in a default state,
        # which might not happen if we abort an experiment in the middle and don't reset it
        delay(1*ms)
        self.experiment.ttl_repump_switch.off() # allow RF to get to the RP AOM
        delay(1*ms)
        self.experiment.ttl_microwave_switch.on() # blocks the microwaves after the mixer
        delay(1*ms)

        # turn off all dds channels
        for ch in self.experiment.all_dds_channels:
            ch.sw.off()

        # todo: turn off all Zotino channels?

        self.experiment.core.break_realtime()


# do this so the code above will not actually run when ARTIQ scans the repository
if __name__ == '__main__':
    pass


