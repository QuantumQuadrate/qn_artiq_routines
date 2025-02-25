"""control two K10CR1 units with only one NDSP server"""

from time import time, sleep
from pylablib.devices import Thorlabs  # for Kinesis instrument control
from artiq.experiment import *


class K10CR1_NDSP_Driver():
    """control two K10CR1 units with only one NDSP server"""

    def __init__(self, devices):
        self.motors = {}
        for name, kwargs in devices:

            self.motors[name] = Thorlabs.KinesisMotor(**kwargs)
            print(f"initialized {name} with {kwargs}")

            self.get_jog_parameters(name)   # prints the jog parameters of the device



            self.get_scale(name)            # prints the scale of the device
            self.get_scale_units(name)      # prints the scale unit of the device
            self.get_stage(name)            # prints the name of the stage

            # device initialized to home
            self.home(name)

            # prints the status of the device
            self.get_status(name)

            print(self.motors[name].get_settings())

            print("------------")
        print("Initialization Done!")

    def get_jog_parameters(self, name: str):
        jog_parameters = self.motors[name].get_jog_parameters()
        print("jog: ", jog_parameters)

        # TJogParams(mode='step', step_size=682667, min_velocity=0.0, acceleration=2047962.2766257639,
        #           max_velocity=2048270.5476900148, stop_mode='profiled')

    def get_scale(self, name: str):
        scale = self.motors[name].get_scale()
        print("scale: (position_scale, velocity_scale, acceleration_scale) = ", scale)

    def get_scale_units(self, name: str):
        scale_units = self.motors[name].get_scale_units()
        print("scale units: ", scale_units)

        if scale_units == "step":
            print("** Autodetected a driver but not detected step scale **")

    def get_stage(self, name: str):
        """Return the name of the stage which was supplied by the usr or autodetected."""
        print("name of the stage: ", self.motors[name].get_stage())

    def get_status(self, name: str):
        print("status of the device: ", self.motors[name].get_status())

    def wait_for_status(self, name: str, status: str):
        self.motors[name].wait_for_status(status = status, enabled = True)


    def get_position(self, name: str) -> float:
        """
        Get current position
        If scale == True, return value in physical units
        else,                          in device internal units (steps)
        """
        position = self.motors[name].get_position()
        return float(position)

    def is_moving(self, name: str) -> bool:
        moving = self.motors[name].is_moving()
        return moving

    #todo: what is the difference btw move_by and move_to?
    def move_by(self, degrees, name: str):
        """Move piezo-motor by a given distance (positive or negative)"""
        self.motors[name].move_by(degrees)
        self.motors[name].wait_move()
        print("K10CR1 moving!")

    def move_to(self, degrees, name: str):
        """Move piezo-motor to position (positive or negative)"""
        self.motors[name].move_to(degrees)
        self.motors[name].wait_move()

    def wait_move(self, name: str):
        self.motors[name].wait_move()

    def home(self, name: str):
        """Home the device and wait until it is done
        force = True: force must be set to True
        once the device is homed at least once, it thinks that it is homed all the time.


        * sync = True: wait until homing is done(with the given timeout)
        * force = False: only home if the device isn't homed already
        """
        self.motors[name].home(sync=True, force=True)
        self.motors[name].wait_for_home()

        # sync = True: wait until homing is done(with the given timeout)
        # force = False: only home if the device isn't homed already



    def is_homed(self, name: str) -> bool:
        """Check if the device is homed"""
        homed = self.motors[name].is_homed()
        return homed





# >>>>>>> 39fa172 (don't hardcode sn's in __init__. two rotators hardware example not working)
# >>>>>>> ndsp_simplify
