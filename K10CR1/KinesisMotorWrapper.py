from time import time
from pylablib.devices import Thorlabs  # for Kinesis instrument control
from artiq.experiment import *


class KinesisMotorWrapper:
    """
    wrapper class to explicitly declare return types and add some functions for diagnostics

    Where applicable, the methods here shadow the names of the methods belonging to Thorlabs.KinesisMotor

    For now, the hardware settings can be configured in the APT User program with 'persist to hardware' checked.
    See ./setup
    """

    def __init__(self, *args, **kwargs):
        self.motor = Thorlabs.KinesisMotor(*args, **kwargs)

    @rpc
    def get_position(self) -> TFloat:
        position = self.motor.get_position()
        return position

    @rpc
    def is_moving(self) -> TBool:
        moving = self.motor.is_moving()
        return moving

    @rpc
    def move_by(self, degrees):
        self.motor.move_by(degrees)

    @rpc
    def move_to(self, degrees):
        self.motor.move_to(degrees)

    @rpc
    def wait_move(self):
        self.motor.wait_move()

    # additional functions for diagnostics
    @rpc
    def get_position_query_time(self) -> TFloat:
        _ = self.motor.get_position()
        return time()
    