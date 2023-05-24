"""
Map devices defined in a parent experiment to aliases

This class also contains the defaults for devices such as urukul channels.
"""

from artiq.experiment import *

# we assume the dds settings will always start out being those that we
# would use first, e.g. the cooling DDS will default to the MOT settings
# not the readout settings. the variable names for power and frequency must
# exist as a dataset create by ExperimentVariables.
dds_defaults = {
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


class SwitchWrapper:
    # almost certainly a better way to do this

    def __init__(self, urukul):
        self.urukul = urukul

    @kernel
    def on(self):
        self.urukul.sw.on()

    @kernel
    def off(self):
        self.urukul.sw.off()


class UrukulWrapper:

    # todo: do some error handling in case there isn't yet a default in dds_defaults
    def __init__(self, urukul, power, frequency):
        """
        A wrapper for an urukul dds object.

        The conversion from power in dB to volts is done under the hood.

        :param urukul: an urukul dds object
        :param power: power in dB
        :param frequency: frequency in Hz
        """
        self.urukul = urukul
        self.sw = SwitchWrapper(urukul) # so we can turn the dds on/off as usual
        self.power_default = power
        self.frequency_default = frequency

    @kernel
    def init(self):
        self.urukul.init()

    @kernel
    def set_att(self, x=0.0):
        """
        :param x: float value. 0.0 by default
        :return:
        """
        self.urukul.set_att(x)

    @kernel
    def dB_to_V(self, dB: TFloat) -> TFloat:
        """
        convert dB to volt amplitude for the dds
        :param dB:
        :return:
        """
        return (2 * 50 * 10 ** (dB / 10 - 3))**(1/2) # can't use math.sqrt in a kernel function

    @kernel
    def set(self, frequency=0.0, power=0.0):
        """
        A wrapper for the urukul's set method that allows setting power in dB

        The actual call is:
            urukul.set(frequency=frequency, amplitude=self.dB_to_V(power))

        :param frequency: optional. if not supplied, default from dds_defaults used
        :return power: optional. if not supplied, default from dds_defaults used
        """

        # might be good to remove this, so we just use set the same way
        # we always use set
        frequency = frequency if frequency != 0.0 else self.frequency_default
        power = power if power != 0.0 else self.power_default

        self.urukul.set(frequency=frequency, amplitude=self.dB_to_V(power))


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

        self.experiment = experiment
        self.dds_list = [] # internal list of references to the dds objects

        for alias in device_aliases:
            if alias in self.alias_map.keys():
                try:

                    dev_name = self.alias_map[alias]

                    # setattr for the device, using the device name from device_db.py.
                    # not self.experiment because we are adding the device to the
                    # experiment we passed by reference.
                    experiment.setattr_device(dev_name)

                    # make an attribute named alias which points to the device object
                    print(f"setattr self.{alias} = self.{dev_name}")

                    if dev_name[:6] == 'urukul':
                        # setup the urukul. initialize and set the power/freq later
                        dev_ref = UrukulWrapper(getattr(experiment, dev_name),
                                                power=getattr(experiment, dds_defaults[alias]['power']),
                                                frequency=getattr(experiment, dds_defaults[alias]['frequency']))
                        self.dds_list.append(dev_ref)
                    else:
                        dev_ref = getattr(experiment, dev_name)

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
        self.experiment.zotino0.init()

        for dds in self.dds_list:
            dds.init()
            dds.set_att()  # set attenuator to 0
            dds.set()  # set power/freq to defaults