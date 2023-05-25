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
    "dds_FORT": {"frequency":"f_FORT", "power":"p_FORT_loading"},
    "dds_cooling_DP": {"frequency":"f_cooling_DP_MOT", "power":"p_cooling_DP_MOT"},
    "dds_cooling_SP": {"frequency": "f_cooling_SP", "power": "p_cooling_SP"},
    "dds_MOT_RP": {"frequency": "f_MOT_RP", "power": "p_MOT_RP"},
    "dds_AOM_A1": {"frequency": "AOM_A1_freq", "power": "AOM_A1_power"},
    "dds_AOM_A2": {"frequency": "AOM_A2_freq", "power": "AOM_A2_power"},
    "dds_AOM_A3": {"frequency": "AOM_A3_freq", "power": "AOM_A3_power"},
    "dds_AOM_A4": {"frequency": "AOM_A4_freq", "power": "AOM_A4_power"},
    "dds_AOM_A5": {"frequency": "AOM_A5_freq", "power": "AOM_A5_power"},
    "dds_AOM_A6": {"frequency": "AOM_A6_freq", "power": "AOM_A6_power"},
}

ALIAS_MAP = {
        "dds_FORT": "urukul0_ch0",
        "dds_cooling_DP": "urukul0_ch1",
        "dds_cooling_SP": "urukul0_ch2",
        "dds_MOT_RP": "urukul0_ch3",
        "dds_AOM_A2": "urukul1_ch0",
        "dds_AOM_A3": "urukul1_ch1",
        "dds_AOM_A1": "urukul1_ch2",
        "dds_AOM_A6": "urukul1_ch3",
        "dds_AOM_A4": "urukul2_ch0",
        "dds_AOM_A5": "urukul2_ch1"
}

class DeviceAliases:

    def __init__(self, experiment, device_aliases):

        self.experiment = experiment
        self.dds_list = [] # internal list of references to the dds objects
        self.dds_powers = []
        self.dds_frequencies = []

        for alias in device_aliases:
            if alias in ALIAS_MAP.keys():
                try:

                    dev_name = ALIAS_MAP[alias]

                    # setattr for the device, using the device name from device_db.py.
                    # not self.experiment because we are adding the device to the
                    # experiment we passed by reference.
                    experiment.setattr_device(dev_name)

                    # make an attribute named alias which points to the device object
                    # print(f"setattr self.{alias} = self.{dev_name}")
                    dev_ref = getattr(experiment, dev_name)

                    if dev_name[:6] == 'urukul':
                        self.dds_list.append(dev_ref)
                        self.dds_powers.append(getattr(experiment, DDS_DEFAULTS[alias]['power']))
                        self.dds_frequencies.append(getattr(experiment, DDS_DEFAULTS[alias]['frequency']))
                    else:
                        pass

                    setattr(experiment, alias, dev_ref)

                except KeyError:
                    print(f"KeyError: {alias} not defined in alias map. Please define it.")

    @kernel
    def initialize(self):

        # we'll assume that each experiment will want to use these
        # todo: could put these as defaults in a keyword arg
        self.experiment.core.reset()
        self.experiment.urukul0_cpld.init()
        self.experiment.urukul1_cpld.init()
        self.experiment.urukul2_cpld.init()
        self.experiment.core.break_realtime()
        self.experiment.zotino0.init()

        for i in range(len(self.dds_list)):
            self.dds_list[i].init()
            self.dds_list[i].set_att(0.0)  # set attenuator to 0
            self.dds_list[i].set(frequency=self.dds_frequencies[i],
                    amplitude=(2*5*10**(self.dds_powers[i]/10 - 3))**(1/2))