from artiq.experiment import *

def dB_to_V(dB: TFloat) -> TFloat:
    """
    convert power in dB to volts for setting DDS amplitudes
    :return amplitude: float in volts
    """
    return (2 * 50 * 10 ** (dB / 10 - 3)) ** (1 / 2)

@kernel
def dB_to_V_kernel(dB: TFloat) -> TFloat:
    """
    convert power in dB to volts for setting DDS amplitudes
    :return amplitude: float in volts
    """
    return (2 * 50 * 10 ** (dB / 10 - 3)) ** (1 / 2)