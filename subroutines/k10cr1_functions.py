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

