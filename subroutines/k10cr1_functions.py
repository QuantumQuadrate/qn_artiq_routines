"""
2025/03/04 EO

Wrapper functions to control K10CR1 Thorlabs waveplate rotator

* Why we need this?
if you want to call NDSP functions from the kernel and return things,
you need to define wrapper functions. if you call the NDSP functions
directly, you'll get an error saying that the expression of NoneType
can't be unified with float (or whatever your return type is).

* Notes on @rpc decorator
    - @rpc decorator is required when calling host functions from a kernel
    - @rpc functions runs on the host PC, while the rest of the run() runs on the FPGA

* Notes on @rpc(flags={"async"})
    - makes RPCs asynchronous
    - allows the kernel to continue execution without waiting for the host function to return


Few ways to improve the code:
1. implement "asyncio"
- this will allow rotating two waveplates at the same time.

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

@rpc
def wr_is_homed(self, name: TStr) -> TBool:
    """
    There is a built_in function that can be called by:
            self.k10cr1_ndsp.is_homed(name)
    However, this function will always return True, if the device has been homed at least once.
    """
    homed = self.k10cr1_ndsp.is_homed(name)
    return homed

@rpc
def is_rotator_moving(self, name: TStr) -> TBool:
    """wrapper function"""
    is_moving = self.k10cr1_ndsp.is_moving(name)
    return is_moving

@rpc
def go_to_home(self, name: TStr):
    """
    k10cr1_ndsp.home() includes wait_for_home() inside - defined in k10cr1_driver.py

    * problem of this built_in function is that, it doesn't know the fastest way to home.
        ex) if the initial position is at 350, it rotates -350 to go back to home rather than +10.
    * good thing is that we can check the position of the rotator and then let the device know which is the fastest way to go.
    * this way, we can save maximum of 18s.
    * oh no... this gives underflow error
    """
    # # delay(2*s)
    # # delay(4*s)
    #
    # pos_in_deg = self.get_rotator_deg(self, name)
    # pos_in_deg %= 360   # this will return 0~359.99999
    #
    # delay(10 * s)
    #
    # if pos_in_deg < 180:
    #     move_by_deg(self, name, -1 * pos_in_deg)
    # else:
    #     move_by_deg(self, name, (360 - pos_in_deg))

    self.k10cr1_ndsp.home(name)

@rpc
def get_rotator_position(self, name: TStr) -> TFloat:
    """
    returns the actual position of the rotator in device units.
    """
    positions = self.k10cr1_ndsp.get_position(name)
    return positions

@rpc
def get_rotator_deg(self, name: TStr) -> TFloat:
    """
    returns the actual position of the rotator in device units.
    """
    #todo: deg % 360
    positions = self.k10cr1_ndsp.get_position(name)
    degs = positions / self.deg_to_pos
    return degs

@rpc
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

@rpc
def move_by_deg(self, name: TStr, target_deg):
    """
    move to target position
    k10cr1_ndsp.move_to() includes wait_move() inside - defined in k10cr1_driver.py
    """

    deg_to_pos = self.deg_to_pos  # sw_position = 136533
    target_pos = int(target_deg * deg_to_pos)

    self.k10cr1_ndsp.move_by(target_pos, name)

@rpc
def wait_move(self, name: TStr):
    self.k10cr1_ndsp.wait_move(name)


@kernel
def record_PDA_power(self):
    """
    This is for the test setup.
    :return: power
    """

    #todo: define this in BaseExperiment.
    # this function was for the test setup.
    # change the name to record_APD_power before merging to the actual setup.
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

    return measurement1

@kernel
def record_FORT_MM_power(self):
    """
    same thing as "measure_FORT_MM_fiber(self)" in experiment_functions.py

    move by increment and measure from the sampler

    :return: power
    """

    measurement_buf = np.array([0.0] * 8)
    measurement1 = 0.0  # 1

    avgs = 50

    # APD
    for i in range(avgs):
        self.sampler1.sample(measurement_buf)
        delay(0.1 * ms)
        measurement1 += measurement_buf[self.FORT_MM_sampler_ch]

    measurement1 /= avgs

    self.append_to_dataset("FORT_MM_monitor", measurement1)

    delay(0.1 * ms)

    return measurement1

@kernel(flags={"fast-math"})
def time_to_rotate_in_ms(self, deg:TFloat) -> TInt32:
    rotator_ave_speed = 10  # deg/s

    if deg < 0:
        deg = -deg

    delay_time = int(deg / rotator_ave_speed * 1000)

    return delay_time