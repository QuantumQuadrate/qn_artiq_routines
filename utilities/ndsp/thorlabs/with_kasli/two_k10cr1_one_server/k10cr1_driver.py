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
