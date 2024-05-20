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
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"

sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")

from subroutines.aom_feedback import AOMPowerStabilizer
from ExperimentVariables import setattr_variables
from utilities.DeviceAliases import DeviceAliases
from utilities.write_h5 import write_results
from utilities.conversions import dB_to_V


class BaseExperiment:

    def __init__(self, experiment: EnvExperiment):
        self.experiment = experiment

    def build(self):
        """
        Put this in your experiment's build method with experiment=self

        Assigning attributes, methods, and devices to the experiment should be done here.

        :param experiment: your experiment.
        :return:
        """

        with open('dataset_db.pyon') as f:
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
                            "core_dma",
                            "scheduler",
                            "urukul0_cpld", "urukul1_cpld", "urukul2_cpld",
                            "zotino0",  # for controlling coils
                            "sampler0",  # for measuring laser power PD
                            "sampler1", # for reading in volts in the coil tune experiment
                            "sampler2",
                            *[f"ttl{i}" for i in range(16)]]
        for dev in devices_no_alias:
            self.experiment.setattr_device(dev)

        # devices can also be nicknamed here:
        self.experiment.ttl_microwave_switch = self.experiment.ttl4
        self.experiment.ttl_repump_switch = self.experiment.ttl5
        self.experiment.ttl_SPCM0 = self.experiment.ttl0
        self.experiment.ttl_scope_trigger = self.experiment.ttl7
        self.experiment.ttl_Luca_trigger = self.experiment.ttl6
        self.experiment.ttl_UV = self.experiment.ttl15
        self.experiment.ttl_SPCM_gate = self.experiment.ttl13

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

        # dataset names
        self.experiment.count_rate_dataset = 'photocounts_per_s'
        self.experiment.scan_var_dataset = "scan_variables"
        self.experiment.scan_sequence1_dataset = "scan_sequence1"
        self.experiment.scan_sequence2_dataset = "scan_sequence2"

        # functions

        @rpc(flags={"async"})
        def print_async(*x):
            """print asynchronously so we don't block the RTIO counter.
            useful for debugging"""
            print(*x)
        self.experiment.print_async = print_async

        def write_results_wrapper(kwargs={}):
            write_results(experiment=self.experiment, **kwargs)
        self.experiment.write_results = write_results_wrapper

        """
        Note that the amplitudes below can be used for setting the urukul channels, but are kernel invariants.
        If you are running the laser_stabilizer in your experiment, and you want to set one of the dds channels we
        feedback to (i.e. it is in one of the dds_feedback_lists in ExperimentVariables), then you should use the
        amplitude attribute of the feedback channel. See subroutines/aom_feedback.py for more details.
        """

        # converts RF power in dBm to amplitudes in V
        self.experiment.ampl_FORT_loading = dB_to_V(self.experiment.p_FORT_loading)
        self.experiment.ampl_cooling_DP_MOT = dB_to_V(self.experiment.p_cooling_DP_MOT)
        self.experiment.ampl_D1_pumping_SP = dB_to_V(self.experiment.p_D1_pumping_SP)
        self.experiment.ampl_pumping_repump = dB_to_V(self.experiment.p_pumping_repump)
        self.experiment.ampl_D1_pumping_SP = dB_to_V(self.experiment.p_D1_pumping_SP)
        self.experiment.ampl_excitation = dB_to_V(self.experiment.p_excitation)
        self.experiment.ampl_microwaves = dB_to_V(self.experiment.p_microwaves)
        self.experiment.ampl_AOM_A1 = dB_to_V(self.experiment.p_AOM_A1)
        self.experiment.ampl_AOM_A2 = dB_to_V(self.experiment.p_AOM_A2)
        self.experiment.ampl_AOM_A3 = dB_to_V(self.experiment.p_AOM_A3)
        self.experiment.ampl_AOM_A4 = dB_to_V(self.experiment.p_AOM_A4)
        self.experiment.ampl_AOM_A5 = dB_to_V(self.experiment.p_AOM_A5)
        self.experiment.ampl_AOM_A6 = dB_to_V(self.experiment.p_AOM_A6)

        # RF powers defined as fractions of the defaults, e.g. the ones we tune during the AOM feedback.
        self.experiment.ampl_FORT_RO = self.experiment.ampl_FORT_loading * self.experiment.p_FORT_RO
        self.experiment.ampl_FORT_PGC = self.experiment.ampl_FORT_loading * self.experiment.p_FORT_PGC
        self.experiment.ampl_FORT_blowaway = self.experiment.ampl_FORT_loading * self.experiment.p_FORT_blowaway
        self.experiment.ampl_FORT_OP = self.experiment.ampl_FORT_loading * self.experiment.p_FORT_OP
        self.experiment.ampl_cooling_DP_RO = self.experiment.ampl_cooling_DP_MOT * self.experiment.p_cooling_DP_RO
        self.experiment.ampl_cooling_DP_PGC = self.experiment.ampl_cooling_DP_MOT * self.experiment.p_cooling_DP_PGC

        # THIS MUST COME LAST IN BASE.BUILD
        # get a list of all attributes of experiment up to this point. if base.build is called in your experiment
        # before any GUI arguments are defined, then this can be used to grab those later by taking a difference
        self.exp_var_names = dir(self.experiment)
        logging.debug("base build - done")


    def set_datasets_from_gui_args(self):
        """
        This should be called at the end of your experiment's build method to archive the GUI arguments.

        For this to work, it is assumed that the line 'self.exp_var_names = dir(self.experiment)' comes
        last in base.build above.
        :return:
        """
        new_exp_var_names = [x for x in dir(self.experiment) if x not in self.exp_var_names]
        for name in new_exp_var_names:
            try:
                self.experiment.set_dataset(name, getattr(self.experiment, name))
            except Exception as e:
                logging.debug(e)

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

        # mainly for cost functions
        try:
            self.experiment.counts_list = [0] * self.experiment.n_measurements
            self.experiment.counts2_list = [0] * self.experiment.n_measurements
        except:
            # if this fails, your experiment probably didn't need it

            logging.warn("experiment does not have variable n_measurements")


        dds_feedback_list = eval(self.experiment.feedback_dds_list)
        slow_feedback_dds_list = eval(self.experiment.slow_feedback_dds_list)
        fast_feedback_dds_list = eval(self.experiment.fast_feedback_dds_list)

        # could implement this but it isn't needed right now
        # self.experiment.slow_laser_stabilizer = AOMPowerStabilizer(experiment=self.experiment,
        #                                                       dds_names=slow_feedback_dds_list,
        #                                                       iterations=self.experiment.aom_feedback_iterations,
        #                                                       averages=self.experiment.aom_feedback_averages,
        #                                                       leave_AOMs_on=True)

        # feedback channels which are fast enough to include both every atom loading attempt.
        # this excludes the on-chip MOT beams because the fW detectors have slow rise time.
        # The external MOT beams and cooling laser could technically be in this list, but
        # why change what isn't broken.
        self.experiment.laser_stabilizer = AOMPowerStabilizer(experiment=self.experiment,
                                                              dds_names=fast_feedback_dds_list,
                                                              iterations=self.experiment.aom_feedback_iterations,
                                                              averages=self.experiment.aom_feedback_averages,
                                                              leave_AOMs_on=True)

        if hasattr(self.experiment, 'n_measurements'):
            self.experiment.set_dataset("n_measurements",self.experiment.n_measurements,broadcast=True)

        logging.debug("base prepare - done")

    @kernel
    def initialize_hardware(self):
        """
        hardware initialization and setting of ttl switches, and set datasets
        :return:
        """

        self.experiment.core.reset()
        self.experiment.set_dataset(self.experiment.scan_var_dataset,'',broadcast=True)
        self.experiment.set_dataset(self.experiment.scan_sequence1_dataset,[0.0],broadcast=True)
        self.experiment.set_dataset(self.experiment.scan_sequence2_dataset,[0.0],broadcast=True)

        self.experiment.named_devices.initialize()

        self.experiment.ttl_microwave_switch.output()
        self.experiment.ttl_repump_switch.output()
        self.experiment.ttl6.output()  # for outputting a trigger
        self.experiment.ttl1.input()

        # channel 3 is configured to read from 14, separated by a switch
        self.experiment.ttl3.input()
        self.experiment.ttl14.output()
        self.experiment.ttl14.on()


        # for diagnostics including checking the performance of fast switches for SPCM gating
        self.experiment.ttl9.output()

        self.experiment.ttl9.off()

        self.experiment.ttl_UV.output()
        self.experiment.ttl_UV.off()

        self.experiment.sampler0.init() # for reading laser feedback
        self.experiment.sampler1.init() # for reading laser feedback
        self.experiment.sampler2.init() # for reading laser feedback


        self.experiment.print_async("base initialize_hardware - done")


        # turn on/off any switches. this ensures that switches always start in a default state,
        # which might not happen if we abort an experiment in the middle and don't reset it
        self.experiment.ttl_repump_switch.off() # allow RF to get to the RP AOM
        delay(1*ms)
        self.experiment.ttl_microwave_switch.on() # blocks the microwaves after the mixer
        delay(1*ms)
        self.experiment.ttl_SPCM_gate.off() # unblocks the SPCM output

        # turn off all dds channels

        for ch in self.experiment.all_dds_channels:
            ch.sw.off()

        # check that the SPCM is plugged in #todo add other SPCM channels as they are added to the experiment
        self.experiment.dds_FORT.sw.on() # we'll get enough Raman scattering to see something
        delay(100*ms)
        t_gate_end = self.experiment.ttl_SPCM0.gate_rising(100*ms)
        counts = self.experiment.ttl_SPCM0.count(t_gate_end)
        delay(10 * ms)
        self.experiment.dds_FORT.sw.off()

        # assert counts > 0, "SPCM0 is likely unplugged"

        # todo: turn off all Zotino channels?
        self.experiment.zotino0.init()
        for zot_ch in range(32):
            self.experiment.zotino0.write_dac(zot_ch, 0.0)
            self.experiment.zotino0.load()
            delay(1*ms)

        self.experiment.print_async("initialize hardware - done")
        self.experiment.core.break_realtime()

        self.experiment.print_async("initialize hardware - done")

# do this so the code above will not actually run when ARTIQ scans the repository
if __name__ == '__main__':
    pass


