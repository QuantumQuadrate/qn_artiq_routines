"""control two K10CR1 units with only one NDSP server"""

from time import time, sleep
from pylablib.devices import Thorlabs  # for Kinesis instrument control
from artiq.experiment import *


class K10CR1_NDSP_Driver():
    """control two K10CR1 units with only one NDSP server"""

    def __init__(self):
        self.motor = {}
        self.motor['780_QWP'] = Thorlabs.KinesisMotor(conn=55000741, scale='K10CR1')
        self.motor['780_HWP'] = Thorlabs.KinesisMotor(conn=55105674, scale='K10CR1')

    def get_position(self, name: str) -> float:
        position = self.motor[name].get_position()
        return position

    def is_moving(self, name: str) -> bool:
        moving = self.motor[name].is_moving()
        return moving

    def move_by(self, degrees, name: str):
        self.motor[name].move_by(degrees)
        print("K10CR1 moving!")

    def move_to(self, degrees, name: str):
        self.motor[name].move_to(degrees)

    # todo: what does this do again? block until movement is done?
    def wait_move(self):
        self.motor.wait_move()
