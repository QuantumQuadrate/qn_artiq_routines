


"""
Control two K10CR1 units with only one NDSP server, where the calls to the
NDSP object are made from a kernel.

1. connect thorlabs K10CR1 rotators to your machine with USB, and make sure
the serial numbers in the device_db match the rotators you connected.

2. open a terminal in the same directory as this script,
enter your artiq virtual environment, then run the
NDSP server for this test with

    python launcher_multi_rotor.py

3. In another terminal in this directory in the artiq environment,
run this script with artiq:

    artiq_run example_multi_k10cr1.py --device-db=.\device_db.py

If you wanted to connect to multiple K10CR1 units, at least if we stick to the
way this example is structured using the K10CR1_NDSP_Driver class, you would need
to launch a separate server for every Thorlabs rotator, which is probably not a
good solution. A better solution would have only one server for all the K10CR1 units,
and this is what the K10CR1_Multi_NDSP_Driver tries to accomplish. See the multi
examples."""

# here's how you do ndsp stuff in your artiq scripts
from artiq.experiment import *
from time import sleep
from pylablib.devices import Thorlabs  # for Kinesis instrument control



import logging
import numpy as np
import sys, os

# cwd = os.getcwd() + "\\"
# sys.path.append(cwd)
# sys.path.append(cwd+"\\repository\\qn_artiq_routines")
# from utilities.BaseExperiment import BaseExperiment



class K10CR1Example_EO(EnvExperiment):

    def build(self):

        self.setattr_device("core")
        self.setattr_device("led0")

        try:
            self.setattr_device("k10cr1_ndsp")
        except Exception as e:
            print(f"Error connecting to device {e}")

        self.setattr_argument("initiallize_to_home", BooleanValue(True))
        self.setattr_argument("target_qwp_deg", NumberValue(0.0, ndecimals=1, step=1))
        self.setattr_argument("target_hwp_deg", NumberValue(0.0, ndecimals=1, step=1))

    # if you want to call NDSP functions from the kernel and return things,
    # you need to define wrapper functions. if you call the NDSP functions
    # directly, you'll get an error saying that the expression of NoneType
    # can't be unified with float (or whatever your return type is).
    def get_rotator_position(self, name: TStr) -> TFloat:
        """wrapper function"""
        positions = self.k10cr1_ndsp.get_position(name)
        return positions

    def is_rotator_moving(self, name: TStr) -> TBool:
        """wrapper function"""
        is_moving = self.k10cr1_ndsp.is_moving(name)
        return is_moving

    def wr_is_homed(self, name: str) -> TBool:
        """wrapper function"""
        homed = self.k10cr1_ndsp.is_homed(name)
        return homed


    def move_to_target_deg(self, name, target_deg):
        """ move to target position """
        deg_to_pos = 136533 # sw_position = 136533
        target_pos = target_deg * deg_to_pos
        self.k10cr1_ndsp.move_to(target_pos, name)
        after_move_pos = self.k10cr1_ndsp.get_position(name)
        after_move_deg = int(after_move_pos / deg_to_pos)
        print("now ", name, " at : deg = ", after_move_deg, "position = ", after_move_pos)



    @kernel
    def experiment_function(self):
        """
        In a real use case, this would be an experiment that lives
        in experiment_functions.py
        """

        self.core.reset()

        # get position and move the 780 QWP.
        qwp780_pos = self.get_rotator_position('780_QWP')
        hwp780_pos = self.get_rotator_position('780_HWP')

        print('780_QWP initially at', qwp780_pos)
        print('780_HWP initially at', hwp780_pos)
        delay(10 * ms)

        #
        # # we can do things conditional on the position, for example move the rotator
        # # or change some output from a Sinara card:

        if qwp780_pos != 0:
            self.k10cr1_ndsp.move_to(0.0, '780_QWP')
            print("780_QWP reset to 0")
        if hwp780_pos != 0:
            self.k10cr1_ndsp.move_to(0.0, '780_HWP')
            print("780_HWP reset to 0")

        # get position and move the 780 QWP.
        qwp780_pos = self.get_rotator_position('780_QWP')
        hwp780_pos = self.get_rotator_position('780_HWP')

        print('780_QWP now at', qwp780_pos)
        print('780_HWP now at', hwp780_pos)
        delay(10 * ms)

        # we can do things conditional on the position, for example move the rotator
        # or change some output from a Sinara card:

        angle = 0

        while angle < 361:
            angle = angle + 10

            self.k10cr1_ndsp.move_to(angle, '780_QWP')
            self.k10cr1_ndsp.move_to(angle, '780_HWP')

            # get position and move the 780 QWP.
            qwp780_pos = self.get_rotator_position('780_QWP')
            hwp780_pos = self.get_rotator_position('780_HWP')

            print("angle at ", angle)
            print('780_QWP at', qwp780_pos)
            print('780_HWP at', hwp780_pos)

        # wait for the rotators to stop moving before proceeding with your experiment
        i = 0
        while self.is_rotator_moving('780_QWP') and self.is_rotator_moving('780_HWP'):
            delay(0.1 * s)
            i += 1
        print(i)
        delay(10 * ms)

    @kernel
    def how_many_trials_to_initialization(self):
        """
        In a real use case, this would be an experiment that lives
        in experiment_functions.py
        """

        self.core.reset()

        # get position and move the 780 QWP & 780 HWP
        qwp780_pos = self.get_rotator_position('780_QWP')
        hwp780_pos = self.get_rotator_position('780_HWP')
        print('780_QWP initially at', qwp780_pos)
        print('780_HWP initially at', hwp780_pos)


        if self.initiallize_to_home:
            print("initialize the devices to home")
            self.k10cr1_ndsp.home('780_QWP')
            self.k10cr1_ndsp.home('780_HWP')

            qwp780_pos = self.get_rotator_position('780_QWP')
            hwp780_pos = self.get_rotator_position('780_HWP')
            print('780_QWP homed at', qwp780_pos)
            print('780_HWP homed at', hwp780_pos)

        #todo: step = 682667; it seems like 1 deg ~ 682667/2/3.14=108704 but not exact.
        # target_qwp = int(self.target_qwp_deg * 682667/2/3.14)
        # target_hwp = int(self.target_hwp_deg * 682667/2/3.14)
        # sw_position = 136533
        # target_qwp = int(self.target_qwp_deg * 136533)


        print('780_QWP target in deg: ', self.target_qwp_deg)
        print('780_HWP target in deg: ', self.target_hwp_deg)

        delay(10 * ms)

        self.move_to_target_deg(name="780_QWP", target_deg=self.target_qwp_deg)
        self.move_to_target_deg(name="780_HWP", target_deg=self.target_hwp_deg)




    @kernel
    def simple_test(self):
        self.core.reset()

        # get position and move the 780 QWP
        qwp780_pos = self.get_rotator_position('780_QWP')
        # hwp780_pos = self.get_rotator_position('780_HWP')


        #todo: Make this as a function; probably put it in initialize_hardware;
        #todo: check how much time it needs for homing the device
        ### Check if the waveplates are homed
        if self.wr_is_homed('780_QWP') != True:
            print("780_QWP is not homed. lets home it before we start exp")
            self.k10cr1_ndsp.home('780_QWP')

        if self.wr_is_homed('780_QWP'):
            print("780_QWP is homed! Can proceed!")



    def initialize_hardware(self):
        print("setting up hardware!")


    def initialize_datasets(self):
        print("setting up datasets!")

    # run is not on the kernel but it calls methods that are
    def run(self):
        """
        hint hint: does this remind you of GeneralVariableScan?
        """

        self.initialize_datasets()
        self.initialize_hardware()

        # our experiment on the kernel
        # self.experiment_function()
        self.how_many_trials_to_initialization()

        # self.simple_test()