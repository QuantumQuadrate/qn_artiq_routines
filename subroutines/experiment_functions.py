from artiq.experiment import *
import numpy as np

@kernel
def test(self):
    x = self.t_blowaway
    self.print_async(x)