"""
Map devices defined in a parent experiment to aliases

This class also contains the defaults for devices such as urukul channels.
"""

from artiq.experiment import *


class DeviceAliases:

    def __init__(self, experiment, device_aliases):

        self.experiment = experiment
        self.dds_list = [] # internal list of references to the dds objects
        self.dds_powers = []
        self.dds_frequencies = []
        self.dds_names_aliases = []

        for alias in device_aliases:
            if alias in self.experiment.alias_map.keys():
                try:

                    dev_name = self.experiment.alias_map[alias]

                    # setattr for the device, using the device name from device_db.py.
                    # not self.experiment because we are adding the device to the
                    # experiment we passed by reference.
                    experiment.setattr_device(dev_name)

                    # make an attribute named alias which points to the device object, e.g.
                    # dev_ref = gettattr(experiment, "urukul0_ch0")
                    dev_ref = getattr(experiment, dev_name)

                    if dev_name[:6] == 'urukul':
                        self.dds_list.append(dev_ref)
                        self.dds_names_aliases.append((alias, dev_name))

                        self.dds_powers.append(getattr(experiment, self.experiment.dds_defaults[alias]['power']))
                        self.dds_frequencies.append(getattr(experiment, self.experiment.dds_defaults[alias]['frequency']))
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
            ampl = (2*50*10**(self.dds_powers[i]/10 - 3))**(1/2) # convert dBm to volts
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