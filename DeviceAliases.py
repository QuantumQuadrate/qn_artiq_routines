"""
Map devices defined in a parent experiment to aliases

This class also contains the defaults for devices such as urukul channels.
"""

from artiq.experiment import *
from math import sqrt

# we assume the dds settings will always start out being those that we
# would use first, e.g. the cooling DDS will default to the MOT settings
# not the readout settings. the variable names for power and frequency must
# exist as a dataset create by ExperimentVariables.
dds_defaults = {
    "dds_FORT": {"frequency":"f_FORT", "power":"p_FORT_loading"},
    "dds_cooling_DP": {"frequency":"f_cooling_DP_MOT", "power":"p_cooling_DP_MOT"},
    "dds_cooling_SP": {"frequency": "f_cooling_SP", "power": "p_cooling_SP"},
    "dds_MOT_RP": {"frequency": "f_MOT_RP", "power": "p_MOT_RP"},
    "dds_AOM_A1": {"frequency": "AOM_A1_freq", "power": "AOM_A1_freq"},
    "dds_AOM_A2": {"frequency": "AOM_A2_freq", "power": "AOM_A2_freq"},
    "dds_AOM_A3": {"frequency": "AOM_A3_freq", "power": "AOM_A3_freq"},
    "dds_AOM_A4": {"frequency": "AOM_A4_freq", "power": "AOM_A4_freq"},
    "dds_AOM_A5": {"frequency": "AOM_A5_freq", "power": "AOM_A5_freq"},
    "dds_AOM_A6": {"frequency": "AOM_A6_freq", "power": "AOM_A6_freq"},
}

def init_devices(experiment):

    # we'll assume that each experiment will want to use these
    experiment.core.reset()
    experiment.urukul0_cpld.init()
    experiment.urukul1_cpld.init()
    experiment.urukul2_cpld.init()
    experiment.core.break_realtime()
    experiment.zotino0.init()

    for dds_name in experiment.dds_list:
        dds = getattr(experiment, dds_name)
        dds.init()
        dds.set_att() # set attenuator to 0
        dds.set() # set power/freq to defaults


class UrukulWrapper:

    # todo: do some error handling in case there isn't yet a default in dds_defaults
    def __init__(self, urukul, power, frequency):
        """
        A wrapper for an urukul dds object.

        Instantiation takes care of initialization, and setting the attenuation,
        power, and frequency.

        :param urukul: an urukul dds object
        :param power: power in dB
        :param frequency: frequency in Hz
        """
        self.urukul = urukul
        self._power = power
        self._frequency = power

    @property
    def power(self):
        """Return the power"""
        return self._power

    @property
    def frequency(self):
        """Return the power"""
        return self._frequency

    def init(self):
        self.urukul.init()

    def set_att(self, x=0.0):
        """
        :param x: float value. 0.0 by default
        :return:
        """
        self.urukul.set_att(x)

    def dB_to_V(self, dB):
        """
        convert dB to volt amplitude for the dds
        :param dB:
        :return:
        """
        return sqrt(2 * 50 * 10 ** (dB / 10 - 3))

    def set(self, frequency=None, power=None):
        """
        A set method that allows only setting the parameter we want to change

        Just a wrapper for the urukul's set method. The actual call is:
            urukul.set(frequency=frequency, amplitude=dB_to_V(power))s
        """
        frequency = frequency if frequency != None else self.frequency
        power = power if power != None else self.power

        self.urukul.set(frequency=frequency, amplitude=dB_to_V(power))
        delay(1*ms)

    def on(self):
        self.urukul.sw.on()

    def off(self):
        self.urukul.sw.off()


class DeviceAliases:

    alias_map = {
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

    def __init__(self, experiment, device_aliases):

        # todo: handle the urukul cpld's. what even are they

        experiment.dds_list = [] # aliases of the dds channels

        for alias in device_aliases:
            if alias in self.alias_map.keys():
                try:

                    dev_name = self.alias_map[alias]

                    # setattr for the device, using the device name from device_db.py
                    experiment.setattr_device(dev_name)

                    # make an attribute named alias which points to the device object
                    print(f"setattr self.{alias} = self.{dev_name}")

                    if dev_name[:6] == 'urukul':
                        # setup the urukul. initialize and set the power/freq later
                        dev_ref = UrukulWrapper(dev_name,
                                                power=dds_defaults[alias]['power'],
                                                frequency=dds_defaults[alias]['frequency'])
                        experiment.dds_list.append(alias)
                    else:
                        dev_ref = getattr(experiment, dev_name)

                    setattr(experiment, alias, dev_ref)

                except KeyError:
                    print(f"KeyError: {alias} not defined in alias map. Please define it.")