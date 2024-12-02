"""control two K10CR1 units with only one NDSP server"""

from time import time, sleep
from pylablib.devices import Thorlabs  # for Kinesis instrument control
from artiq.experiment import *


class K10CR1_NDSP_Driver():
    """control two K10CR1 units with only one NDSP server"""

<<<<<<< HEAD
    def __init__(self, devices):
        self.motors = {}
        for name, kwargs in devices:
                self.motors[name] = Thorlabs.KinesisMotor(**kwargs)
                print(f"initialized {name} with {kwargs}")
=======
<<<<<<< HEAD
    def __init__(self):
        self.motor = {}
        self.motor['780_QWP'] = Thorlabs.KinesisMotor(conn=55000741, scale='K10CR1')
        self.motor['780_HWP'] = Thorlabs.KinesisMotor(conn=55105674, scale='K10CR1')
>>>>>>> ndsp_simplify

    def get_position(self, name: str) -> float:
        position = self.motors[name].get_position()
        return float(position)

    def is_moving(self, name: str) -> bool:
        moving = self.motors[name].is_moving()
        return moving

    def move_by(self, degrees, name: str):
        self.motors[name].move_by(degrees)
        print("K10CR1 moving!")

    def move_to(self, degrees, name: str):
        self.motors[name].move_to(degrees)

    # todo: what does this do again? block until movement is done?
    def wait_move(self):
<<<<<<< HEAD
        self.motors.wait_move()
=======
        self.motor.wait_move()
=======
    def __init__(self, devices):
        self.motors = {}
        for name, kwargs in devices:
                self.motors[name] = Thorlabs.KinesisMotor(**kwargs)
                print(f"initialized {name} with {kwargs}")

    def get_position(self, name: str) -> float:
        position = self.motors[name].get_position()
        return float(position)

    def is_moving(self, name: str) -> bool:
        moving = self.motors[name].is_moving()
        return moving

    def move_by(self, degrees, name: str):
        self.motors[name].move_by(degrees)
        print("K10CR1 moving!")

    def move_to(self, degrees, name: str):
        self.motors[name].move_to(degrees)

    # todo: what does this do again? block until movement is done?
    def wait_move(self):
        self.motors.wait_move()
>>>>>>> 39fa172 (don't hardcode sn's in __init__. two rotators hardware example not working)
>>>>>>> ndsp_simplify
