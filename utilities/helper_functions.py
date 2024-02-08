from artiq.experiment import *

@rpc(flags={"async"})
def print_async(*x):
    """print asynchronously so we don't block the RTIO counter.
    useful for debugging"""
    print(*x)