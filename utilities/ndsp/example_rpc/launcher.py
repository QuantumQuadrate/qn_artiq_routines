# you can install sipyco via nix, or from the github using pip
from sipyco.pc_rpc import simple_server_loop

# local imports
import example_driver
from device_db import device_db

def main():
    # network settings
    device_name = "example_ndsp"
    host = device_db[device_name]["host"]  # run it locally
    port = device_db[device_name]["port"]

    # instantiate driver object
    driver = example_driver.ExampleDriver("I'm a little example :3")

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