"""
This is a simple code to read tune_sync_delay and tune_io_update_delay for each urukul channel individually.
The results can be used in device_db for each channel for calibrating the dds.

tune_sync_delay outputs like (15, 2) only the first number will be used. The 2nd number is a width or something.
tune_io_update_delay outputs a number, often either 0 or 1.

These values can also be obtained when running "artiq_sinara_tester -o urukuls". However, you will need to enable the EEPROM
in device_db for each channel, like the following:

device_db["urukul2_ch0"] = {
    "type": "local",
    "module": "artiq.coredevice.ad9910",
    "class": "AD9910",
    "arguments": {
        "pll_en": 1,
        "pll_n": 32,
        "chip_select": 4,
        "cpld_device": "urukul2_cpld",
        "sw_device": "ttl_urukul2_sw0",
        "sync_delay_seed": "eeprom_urukul2:64",
        "io_update_delay": "eeprom_urukul2:64",
        # "sync_delay_seed": 15,
        # "io_update_delay": 0
    }
}


"""

from artiq.experiment import *
from artiq.coredevice.ad9910 import (
    PHASE_MODE_ABSOLUTE,
    PHASE_MODE_CONTINUOUS,
    PHASE_MODE_TRACKING,
)

import numpy as np

class UrukulSyncCalibrate(EnvExperiment):

    def build(self):
       self.setattr_device("core")
       self.dds1 = self.get_device("urukul2_ch3")

    @kernel
    def run(self):
        self.core.reset()
        self.dds1.cpld.init()
        self.dds1.init()
        self.dds1.set_att(float(0))

        delay(1 * ms)

        delay_value=self.dds1.tune_sync_delay()
        delay(10*ms)
        print(delay_value)
        delay(10*ms)

        io_value = self.dds1.tune_io_update_delay()
        delay(10 * ms)
        print(io_value)
        delay(10 * ms)

        print("code done!")