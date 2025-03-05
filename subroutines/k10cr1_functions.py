"""
2025/03/04 EO

Wrapper functions to control K10CR1 Thorlabs waveplate rotator

* Why we need this?
if you want to call NDSP functions from the kernel and return things,
you need to define wrapper functions. if you call the NDSP functions
directly, you'll get an error saying that the expression of NoneType
can't be unified with float (or whatever your return type is).

"""


from artiq.experiment import *
import numpy as np
import os, sys
# cwd = os.getcwd() + "\\"
#
# sys.path.append(cwd)
# sys.path.append(cwd+"\\repository\\qn_artiq_routines")

@kernel
def working_example(self):
    """
    ratator_test function to record the Sampler value while rotating the waveplate

    :param self:
    :return:
    """

    self.core.reset()
    delay(2*s)
    # get position and move the 780 QWP.
    qwp780_pos = get_rotator_position(self,'780_QWP')
    hwp780_pos = get_rotator_position(self,'780_HWP')

    self.print_async('780_QWP initially at ', qwp780_pos/self.deg_to_pos, ' deg')
    self.print_async('780_HWP initially at ', hwp780_pos/self.deg_to_pos, ' deg')
    delay(1 * s)

    # home the devices
    go_to_home(self, '780_QWP')
    go_to_home(self, '780_HWP')

    self.print_async("home done")
    delay(2*s)

    qwp780_pos = get_rotator_position(self, '780_QWP')
    hwp780_pos = get_rotator_position(self, '780_HWP')

    self.print_async('780_QWP homed at', qwp780_pos/self.deg_to_pos, ' deg')
    self.print_async('780_HWP homed at', hwp780_pos/self.deg_to_pos, ' deg')

    move_to_target_deg(self, name="780_QWP", target_deg=self.target_qwp_deg)
    move_to_target_deg(self, name="780_HWP", target_deg=self.target_hwp_deg)

    delay(10 * ms)

def wr_is_homed(self, name: TStr) -> TBool:
    """
    There is a built_in function that can be called by:
            self.k10cr1_ndsp.is_homed(name)
    However, this function will always return True, if the device has been homed at least once.
    """
    homed = self.k10cr1_ndsp.is_homed(name)
    return homed

def is_rotator_moving(self, name: TStr) -> TBool:
    """wrapper function"""
    is_moving = self.k10cr1_ndsp.is_moving(name)
    return is_moving


def go_to_home(self, name: TStr):
    """
    k10cr1_ndsp.home() includes wait_for_home() inside - defined in k10cr1_driver.py
    """
    self.k10cr1_ndsp.home(name)


def get_rotator_position(self, name: TStr) -> TFloat:
    """
    returns the actual position of the rotator in device units.
    """
    positions = self.k10cr1_ndsp.get_position(name)
    return positions

def get_rotator_deg(self, name: TStr) -> TFloat:
    """
    returns the actual position of the rotator in device units.
    """
    positions = self.k10cr1_ndsp.get_position(name)
    degs = positions / self.deg_to_pos
    return degs


def move_to_target_deg(self, name: TStr, target_deg):
    """
    move to target position
    k10cr1_ndsp.move_to() includes wait_move() inside - defined in k10cr1_driver.py
    """

    deg_to_pos = self.deg_to_pos # sw_position = 136533
    target_pos = int(target_deg * deg_to_pos)
    self.k10cr1_ndsp.move_to(target_pos, name)

    # after_move_pos = self.k10cr1_ndsp.get_position(name)
    # after_move_deg = int(after_move_pos / deg_to_pos)
    # print("now ", name, " at : deg = ", after_move_deg)


def move_by_deg(self, name: TStr, target_deg):
    """
    move to target position
    k10cr1_ndsp.move_to() includes wait_move() inside - defined in k10cr1_driver.py
    """

    deg_to_pos = self.deg_to_pos  # sw_position = 136533
    target_pos = int(target_deg * deg_to_pos)

    self.k10cr1_ndsp.move_by(target_pos, name)


def wait_move(self, name: TStr):
    self.k10cr1_ndsp.wait_move(name)


@kernel
def record_PDA_power(self):
    """
    move by increment and measure from the sampler

    :return: power
    """

    # Sampler2 ch1
    test_PDA_ch = 1

    measurement_buf = np.array([0.0] * 8)
    measurement1 = 0.0  # 1

    avgs = 50

    # APD
    for i in range(avgs):
        self.sampler2.sample(measurement_buf)
        delay(0.1 * ms)
        measurement1 += measurement_buf[test_PDA_ch]

    measurement1 /= avgs

    self.append_to_dataset("test_PDA_monitor", measurement1)


    delay(0.1 * ms)