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

def main():
    # network settings
    device_name = "k10cr1_ndsp"
    host = device_db[device_name]["host"]
    port = device_db[device_name]["port"]

    # could get the serial number of a device this way:
    # sn = device_db[device_name]["sn_list"][0]

    # for now, we'll just try to connect to the first available device
    sn, name = Thorlabs.list_kinesis_devices()[0]

    # instantiate driver object
    driver = k10cr1_driver.K10CR1_NDSP_Driver()

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