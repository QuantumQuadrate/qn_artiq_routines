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

            self.get_scale(name)
            self.get_scale_units(name)
            self.get_stage(name)
            self.get_status(name)

            print("------------")
        print("Initialization Done!")


    def get_scale(self, name: str):
        scale = self.motors[name].get_scale()
        print("scale: (position_scale, velocity_scale, acceleration_scale) = ", scale)

    def get_scale_units(self, name: str):
        scale_units = self.motors[name].get_scale_units()
        print("scale units: ", scale_units)

        if scale_units == "step":
            print("** Autodetected a driver but not detected step scale **")

    def get_stage(self, name: str):
        print("name of the stage: ", self.motors[name].get_stage())

    def get_status(self, name: str):
        print("status of the device: ", self.motors[name].get_status())



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
        self.motors[name].move_by(degrees)
        print("K10CR1 moving!")

    def move_to(self, degrees, name: str):
        self.motors[name].move_to(degrees)

    def wait_move(self, name: str):
        self.motors[name].wait_move()

    def home(self, name: str):
        """Home the device"""
        self.motors[name].home(sync=True, force=False)
        # sync = True: wait until homing is done(with the given timeout)
        # force = False: only home if the device isn't homed already

    def is_homed(self, name: str) -> bool:
        """Check if homing is in progress"""
        homed = self.motors[name].is_homed()
        return homed

    def wait_for_home(self, name: str):
        """Wait until the device is homed"""
        self.motors[name].wait_for_home()


# >>>>>>> 39fa172 (don't hardcode sn's in __init__. two rotators hardware example not working)
# >>>>>>> ndsp_simplify
