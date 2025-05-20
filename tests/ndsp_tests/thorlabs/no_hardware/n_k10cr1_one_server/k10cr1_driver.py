"""control two K10CR1 units with only one NDSP server"""

from time import time, sleep
from pylablib.devices import Thorlabs  # for Kinesis instrument control
from artiq.experiment import *


class K10CR1_NDSP_Driver():
    """control two K10CR1 units with only one NDSP server"""

    def __init__(self, devices):
        self.motor = {}
        for name, kwargs in devices:
                self.motor[name] = Thorlabs.KinesisMotor(**kwargs)
                print(f"initialized {name} with {kwargs}")
    @rpc
    def get_position(self, names: TList(TStr)) -> TList(TFloat):
        positions = []
        for name in names:
            position = self.motor[name].get_position()
            positions.append(position)
        return positions    

    @rpc
    def is_moving(self, name) -> TBool:
        moving = self.motor[name].is_moving()
        return moving

    @rpc(flags={'async'})
    def move_by(self, degrees: TList(TFloat), names: TList(TStr)):
        for name,deg in zip(names, degrees):
            self.motor[name].move_by(deg)
            print(f"K10CR1 {name} moving by {deg}!")

    # todo
    @rpc(flags={'async'})
    def move_to(self, degrees):
        self.motor.move_to(degrees)

    @rpc(flags={'async'})
    def wait_move(self):
        self.motor.wait_move()
