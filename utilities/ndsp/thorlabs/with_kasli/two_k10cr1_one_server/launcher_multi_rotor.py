#!/usr/bin/env python3

# you can install sipyco via nix, or from the github using pip
from sipyco.pc_rpc import simple_server_loop
from pylablib.devices import Thorlabs  # for Kinesis instrument control
# import sys, os

# # local imports
# cwd = os.getcwd() + "\\"
# sys.path.append(cwd)
# sys.path.append(cwd+"\\repository\\qn_artiq_routines")
import k10cr1_driver as k10cr1_driver
from device_db import device_db

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def main():
    logger.info("Running unit tests")
    # network settings
    device_name = "k10cr1_ndsp"
    host = device_db[device_name]["host"]
    port = device_db[device_name]["port"]

    sn_list = device_db[device_name]["sn_list"]
    nicknames = device_db[device_name]["nickname_list"]
    devices = [(name, {'conn':sn}) for name,sn in zip(nicknames,sn_list)]
    # devices = [
    #     ("780_QWP", {"conn": 55000741}),
    #     ("780_HWP", {"conn": 55105674})
    # ]

    # instantiate driver object
    driver = k10cr1_driver.K10CR1_NDSP_Driver(devices)

    try:
        simple_server_loop(
            {device_name: driver},
            host,
            port
        )
    except Exception as e:
        print("Insert Error handling here")
    finally:
        driver.close()

if __name__ == "__main__":
     main()