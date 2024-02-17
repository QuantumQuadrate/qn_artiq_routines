"""
Script for moving a Thorlabs K10CR1 and reading an analog signal.

Preston Huft, Fall 2022

This script is for 
1. aligning a polarizer held by the K10CR1 to the polarization of an incoming beam using slope detection.
2. aligning polarization of an input beam to a reference angle (found in step 1) about which the K10CR1 will rotate 
a polarizer. That is, the K10CR1 will rotate to -/+ 45 deg wrt a refereence angle that
we want to define the polarization.

Note that the APT Thorlabs drivers must be installed for this to work. If the connection works 
but the stage fails to move, use a shorter cable or avoid using a USB hub.

References:
https://github.com/AlexShkarin/pyLabLib/blob/main/docs/devices/Thorlabs_kinesis.rst#id25
https://pylablib.readthedocs.io/en/latest/.apidoc/pylablib.devices.Thorlabs.html#pylablib.devices.Thorlabs.kinesis.KinesisMotor.get_velocity_parameters
"""

from time import sleep, time
from pylablib.devices import Thorlabs # for Kinesis instrument control
import nidaqmx as daq
import nidaqmx.constants as daq_constants
from nidaqmx.errors import DaqError, DaqWarning
from nidaqmx.error_codes import DAQmxErrors, DAQmxWarnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


## global "constants" - set these before running

DRY_RUN = False # will run the loop without moving the rotator or taking data to test plotting

# if true, align a polarizer held by the K10CR1 to a reference polarized beam. if false, just rotate the polarizer
# by +/-45 degrees about a reference angle and plot the difference in photodetector values at the end points.
# the only difference between these modes is whether the rotator angle recieves feedback.
ALIGN_ROTATOR = False


## generic functions

def figax():
    """
    plot setup
    """
    plt.ion() # this line ensures that the program doesn't hang here until the figure is closed
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax.set_ylim((-1,1))
    ax2.set_ylim((-0.1,0.1))
    ax2.set_xlabel('iter')
    ax.set_ylabel('I(45) - I(-45)')
    fig.add_axes(ax)
    return fig, (ax, ax2)

## setup a connection to a Thorlabs rotator

ser = "55000915" #'55000741'
# devices = Thorlabs.list_kinesis_devices() # can use this to print out the detected devices
stage = Thorlabs.KinesisMotor(ser, scale='K10CR1')
scl = True # whether to use physical units - this apparently has no effect
if not stage.is_homed():
    stage.home()

def print_pos():
    print("position =",stage.get_position(scale=scl))

## create an NIDAQmx task to read from an Analog input device

task = daq.Task()
channel='Dev1/ai0' # the formula is device_name/channel_name where device_name can be found in NI MAX
acq_type = getattr(daq_constants.AcquisitionType, 'FINITE')
terminal_cfg = getattr(daq_constants.TerminalConfiguration,'NRSE')
sample_rate = 50 # per second
samples = 100 # per channel
task.ai_channels.add_ai_voltage_chan(physical_channel=channel)
task.timing.cfg_samp_clk_timing(sample_rate, sample_mode=acq_type, samps_per_chan=samples)

## the main part of the code. move the rotator, take data.

# for measuring polarization with slope detection, we want to rotate a polarizer +/- 45 degrees about
# the desired polarization angle. we will record a photodetector voltage at each of these endpoints 
# several times and take the difference of the average of each. This should be zero when the polarizer
# is aligned to the 0 of the polarizer.

ref_degs = 61.33 # the angle of the polarizer about which to move +/-45 degs
# print_pos()
# stage.move_to(ref_degs)
# stage.wait_move() # don't do anything until the motion is complete
print_pos()
iters = 20 # the maximum number of attempts to align the polarizer
data = np.empty(iters, float)
fig, axes = figax()
ax,ax2 = axes
ax.set_xlim((0, iters))
ax2.set_xlim((0, iters))
plt.show()
line, = ax.plot(np.zeros(iters), np.zeros(iters), color='crimson')
line2, = ax2.plot(np.zeros(iters), np.zeros(iters), color='crimson')

proportional = 20 # P in our feedback loop
tolerance_deg = 0.01 # the amount of tolerable error; break out of the loop if we get there
err_degs = 360

 # for the dry run simulation
if DRY_RUN:

    if ALIGN_ROTATOR:
        target_degs = 10 # the angle of the polarization we are aligning to
    else:
        target_degs = ref_degs
    
    # mimic the laser noise
    fast_noise = lambda x: 0.2*np.sin(2*np.pi*20*x - 0.01*(np.random.rand()-0.5))
    slow_noise = lambda x: 0.5*np.sin(2*np.pi*0.1*x - 0.01*(np.random.rand()-0.5))
    noise = lambda x: slow_noise(x) + fast_noise(x)
     
for i in range(iters):
    
    # move the stage and take data
    if not DRY_RUN:
        
        stage.move_to(ref_degs - 45)
        stage.wait_move()
        print_pos()
        data[i] = np.mean(task.read(number_of_samples_per_channel=samples)) # take data
        
        stage.move_to(ref_degs + 45)
        stage.wait_move()
        print_pos()
        data[i] -= np.mean(task.read(number_of_samples_per_channel=samples)) # take data
        
        
        
        if ALIGN_ROTATOR:
            ref_degs -= proportional*data[i]
            
            # if abs(data[i] - data[i-1]) < deviation:
                # print(f"reached target in {i} iterations")
                # print(f"ref angle = {ref_degs} degs")
                # break
            if abs(err_degs) < tolerance_deg:
                print(f"reached target error in {i+1} iterations")
                print(f"ref angle = {ref_degs} degs")
                break
    
    # simulate moving the stage and taking data
    else: 
        
        error_degs = target_degs-ref_degs
        data[i] = np.cos(np.pi/4+error_degs*np.pi/180)**2 + noise(i)
        data[i] -= np.cos(np.pi/4-error_degs*np.pi/180)**2 + noise(i)
        print("faking rotator movement delay. are we fooling you?")
        sleep(1)
        
        if ALIGN_ROTATOR:
            # update the rotator position using the error
            ref_degs -= proportional*data[i]
            
            # if abs(data[i] - data[i-1]) < deviation:
                # print(f"reached target in {i} iterations")
                # print(f"ref angle = {ref_degs} degs")
                # print(f"error = {target_degs - ref_degs} deg")
                # break
                
            if abs(err_degs) < tolerance_deg:
                print(f"reached target error in {i+1} iterations")
                break
    
    # estimated error (valid for small error)
    err_degs = -0.5*data[i]*180/np.pi # from Taylor expanding and solving for error_degs
    
    # update the plot
    err_title = f"estimated error = {err_degs} deg"
    ax.set_title(err_title)
    line.set_data(range(i+1), data[:i+1])
    line2.set_data(range(i+1), data[:i+1])
    fig.canvas.draw()
    fig.canvas.flush_events()
    
    # let user know it's safe to adjust the polarizer now because we're between measurements
    if not ALIGN_ROTATOR:
        fig.set_facecolor('green') # okay to change GRIN tube polarizer while
        ax.set_title(f"Safe to adjust the polarizer. err = {err_degs} deg")
        fig.canvas.draw()
        fig.canvas.flush_events()
        sleep(10)
        fig.set_facecolor('white')
        ax.set_title(err_title)
        fig.canvas.draw()
        fig.canvas.flush_events()
    
print(f"final error = {err_degs} degs")
print(f"ref angle = {ref_degs} degs")
    
# it's important to clean up when we're done
plt.show() # keep the plot open
# plt.close()
stage.close()
task.stop()
task.close()