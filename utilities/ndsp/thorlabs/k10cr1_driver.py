from time import time, sleep
from pylablib.devices import Thorlabs  # for Kinesis instrument control
from artiq.experiment import *

# driver class. We'll expose the functions of this class
# after it's been instantiated

class K10CR1_NDSP_Driver():

    def __init__(self, **kwargs):
        self.motor = Thorlabs.KinesisMotor(**kwargs)

    @rpc
    def get_position(self) -> TFloat:
        position = self.motor.get_position()
        return position

    @rpc
    def is_moving(self) -> TBool:
        moving = self.motor.is_moving()
        return moving

    @rpc(flags={'async'})
    def move_by(self, degrees):
        self.motor.move_by(degrees)
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

