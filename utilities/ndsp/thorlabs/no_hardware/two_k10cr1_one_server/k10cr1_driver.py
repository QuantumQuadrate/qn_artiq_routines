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

    @rpc
    def get_position(self, name) -> TFloat:
        position = self.motor[name].get_position()
        return position

    @rpc
    def is_moving(self, name) -> TBool:
        moving = self.motor[name].is_moving()
        return moving

    @rpc(flags={'async'})
    def move_by(self, degrees, name):
        self.motor[name].move_by(degrees)
        print("K10CR1 moving!")

    @rpc(flags={'async'})
    def move_to(self, degrees):
        self.motor.move_to(degrees)

    @rpc(flags={'async'})
    def wait_move(self):
        self.motor.wait_move()

    # additional functions for diagnostics
    @rpc
    def get_position_query_time(self) -> TFloat:
        _ = self.motor.get_position()
        return time()
