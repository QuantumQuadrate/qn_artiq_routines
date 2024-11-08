# driver class. We'll expose the functions of this class
# after it's been instantiated

from artiq.experiment import *

class ExampleDriver():

    def __init__(self, info):
        self.info = info

    def do_something(self):
        print("I'm alive!")

    def return_a_number(self):
        return 3.14

    def close(self):
        print("bye!")