Controlling Thorlabs devices with a NDSP

Single Thorlabs K10CR1 example:

1. connect a thorlabs K10CR1 rotator to your machine with USB

2. open a terminal in the same directory as this file, 
enter your artiq virtual environment, then run the 
NDSP server for this test with

    python launcher_single_rotor.py

3. In another terminal in this directory in the artiq environment,
run this script with artiq:

    artiq_run ndsp_k10cr1_example.py --device-db=.\device_db.py