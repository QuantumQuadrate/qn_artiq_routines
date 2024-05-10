"""
Map devices defined in a parent experiment to aliases

This class also contains the defaults for devices such as urukul channels.
"""

from artiq.experiment import *


# we assume the dds settings will always start out being those that we
# would use first, e.g. the cooling DDS will default to the MOT settings
# not the readout settings. the variable names for power and frequency must
# exist as a dataset create by ExperimentVariables.
DDS_DEFAULTS = {
    "dds_FORT": {"frequency":"f_FORT", "power":"p_FORT_default"},
    "dds_cooling_DP": {"frequency":"f_cooling_DP_MOT", "power":"p_cooling_DP_MOT"},
    "dds_D1_pumping_SP": {"frequency": "f_D1_pumping_SP", "power": "p_D1_pumping_SP"},
    "dds_pumping_repump": {"frequency": "f_pumping_repump", "power": "p_pumping_repump"},
    "dds_excitation": {"frequency": "f_excitation", "power": "p_excitation"},
    "dds_microwaves": {"frequency": "f_microwaves_dds", "power": "p_microwaves"},
    "dds_AOM_A1": {"frequency": "AOM_A1_freq", "power": "p_AOM_A1"},
    "dds_AOM_A2": {"frequency": "AOM_A2_freq", "power": "p_AOM_A2"},
    "dds_AOM_A3": {"frequency": "AOM_A3_freq", "power": "p_AOM_A3"},
    "dds_AOM_A4": {"frequency": "AOM_A4_freq", "power": "p_AOM_A4"},
    "dds_AOM_A5": {"frequency": "AOM_A5_freq", "power": "p_AOM_A5"},
    "dds_AOM_A6": {"frequency": "AOM_A6_freq", "power": "p_AOM_A6"},
}

ALIAS_MAP = {
        "dds_FORT": "urukul0_ch0",
        "dds_cooling_DP": "urukul0_ch1",
        "dds_D1_pumping_SP": "urukul2_ch2",
        "dds_pumping_repump": "urukul0_ch2",
        "dds_AOM_A2": "urukul1_ch0",
        "dds_AOM_A3": "urukul1_ch1",
        "dds_AOM_A1": "urukul1_ch2",
        "dds_AOM_A6": "urukul1_ch3",
        "dds_AOM_A4": "urukul2_ch0",
        "dds_AOM_A5": "urukul2_ch1",
        "dds_excitation": "urukul0_ch3",
        "dds_microwaves": "urukul2_ch3"
}

@rpc(flags={"async"})  # means this code runs asynchronously; won't block the rtio counter
def print_async(x):
    print(*x)


class DeviceAliases:

    def __init__(self, experiment, device_aliases):

        self.experiment = experiment
        self.dds_list = [] # internal list of references to the dds objects
        self.dds_powers = []
        self.dds_frequencies = []
        self.dds_names_aliases = []

        for alias in device_aliases:
            if alias in ALIAS_MAP.keys():
                try:

                    dev_name = ALIAS_MAP[alias]

                    # setattr for the device, using the device name from device_db.py.
                    # not self.experiment because we are adding the device to the
                    # experiment we passed by reference.
                    experiment.setattr_device(dev_name)

                    # make an attribute named alias which points to the device object, e.g.
                    # dev_ref = gettattr(experiment, "urukul0_ch0")
                    dev_ref = getattr(experiment, dev_name)

                    if dev_name[:6] == 'urukul':
                        self.dds_list.append(dev_ref)
                        self.dds_names_aliases.append((alias,dev_name))

                        # print(dev_name, alias, DDS_DEFAULTS[alias]['power'], DDS_DEFAULTS[alias]['frequency'])
                        self.dds_powers.append(getattr(experiment, DDS_DEFAULTS[alias]['power']))
                        self.dds_frequencies.append(getattr(experiment, DDS_DEFAULTS[alias]['frequency']))
                    else:
                        pass

                    # set an experiment attribute with the name alias that points to the
                    # device object, e.g.
                    # settattr(experiment, "dds_FORT", getattr(experiment,"urukul0_ch0"))
                    setattr(experiment, alias, dev_ref)

                except KeyError:
                    print(f"KeyError: {alias} not defined in alias map. Please define it.")

    @kernel
    def initialize(self):

        # we'll assume that each experiment will want to use these
        self.experiment.core.reset()
        self.experiment.urukul0_cpld.init()
        self.experiment.urukul1_cpld.init()
        self.experiment.urukul2_cpld.init()
        self.experiment.core.break_realtime()

        for i in range(len(self.dds_list)):
            self.dds_list[i].init()
            self.dds_list[i].set_att(float(0))  # set attenuator to 0
            # print_async(self.dds_names_aliases[i])
            ampl = (2*50*10**(self.dds_powers[i]/10 - 3))**(1/2) # convert dBm to volts
            # print_async([self.dds_powers[i], ampl])
            self.dds_list[i].set(frequency=self.dds_frequencies[i],
                    amplitude=ampl)
            delay(1*ms)

    @kernel
    def set_dds_default_settings(self):
        """
        sets the dds frequencies and amplitudes to the defaults defined by experiment variables

        WARNING: this does not get the most up-to-date values from the datasets, which may have
        changed if aom_feedback was running
        :return:
        """

        # for i in range(len(self.dds_list)):
        #     ampl = (2*50*10**(self.dds_powers[i]/10 - 3))**(1/2) # convert dBm to volts
        #     self.dds_list[i].set(frequency=self.dds_frequencies[i],
        #             amplitude=ampl)
        #     delay(1*ms)

        for i in range(len(self.dds_list)):
            ampl = (2*50*10**(self.dds_powers[i]/10 - 3))**(1/2) # convert dBm to volts
            self.dds_list[i].set(frequency=self.dds_frequencies[i],
                    amplitude=ampl)
            delay(1*ms)