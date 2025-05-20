import artiq *
import pylablib *

class Thorlabs(EnvExperiment):

    def build(self):


    def prepare(self):
        self.motor = ThorlabsMotor()
        # setup ...
    rpc
    def move(self):
        self.motor = ThorlabsMotor()

        self.motor.move_by(x)

        self.motor.move_to(x)


        self.motor.close()

    @kernel
    def run(self):

        self.dds_FORT.sw.on()

        self.move()

        self.dds_FORT.sw.off()