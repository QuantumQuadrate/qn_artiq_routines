"""
Trigger the Luca camera to take images of MOT5,6 while the other is pulsed off
and capture the MOT using the ThorCam which images through the JenOptik objective.

ThorCam code adapted from Thorlabs' Camera_Examples\Python\grab_single_frame.py
"""

from artiq.experiment import *

import numpy as np
import os
import cv2
from PIL import Image
from thorlabs_tsi_sdk.tl_camera import TLCameraSDK, OPERATION_MODE
import matplotlib.pyplot as plt
from datetime import datetime as dt

import sys, os
# get the current working directory
current_working_directory = os.getcwd()
cwd = os.getcwd() + "\\"
sys.path.append(cwd)
sys.path.append(cwd+"\\repository\\qn_artiq_routines")
from utilities.BaseExperiment import BaseExperiment

NUM_FRAMES = 1  # adjust to the desired number of frames
dll_parent_dir = '/qn_artiq_routines/third_party/thorlabs'

os.add_dll_directory(dll_parent_dir + "\\dlls")
os.environ['PATH'] = dll_parent_dir + "\\dlls\\" + os.pathsep + os.environ['PATH']


class MonitorMOTandExternalBeamPositions(EnvExperiment):

    def build(self):
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("run_time_minutes", NumberValue(60))
        self.setattr_argument("t_ThorCam_exposure_ms", NumberValue(100))

        # take images every image_frequency iterations
        self.setattr_argument("image_frequency", NumberValue(10))
        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        self.AOM_feedback_period_cycles = 10
        self.last_photocounts_to_plot = 200 # how many datapoints will be broadcast

        self.n_steps = int(60 * self.run_time_minutes / (
                self.t_MOT_loading + self.t_FORT_loading + 0.5 / self.AOM_feedback_period_cycles) + 1)

        self.t_Luca_exposure = 10*ms

        # where to store the data
        self.t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
        self.data_dir = os.getcwd()
        self.file_suffix = 'MonitorMOTandExternalBeamPositions'
        self.file_prefix = os.path.join(self.data_dir, self.t_experiment_run + '_' + self.file_suffix)

    @kernel
    def run(self):
        self.base.initialize_hardware()

        self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
                             channels=self.coil_channels)

        self.set_dataset("photocount_bins", [50], broadcast=True)
        self.set_dataset("photocounts", [0])
        self.set_dataset("photocounts_current_iteration", [0], broadcast=True, archive=False)
        self.set_dataset("FORT_TTI_volts", [0.0])

        # turn off AOMs we aren't using, in case they were on previously
        self.dds_D1_pumping_DP.sw.off()
        self.dds_excitation.sw.off()

        # turn on cooling MOT AOMs and warm up FORT AOM
        self.dds_cooling_DP.sw.on()  # cooling double pass
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()

        delay(2000 * ms)  # wait for AOMS to thermalize in case they have been off.

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()
        delay(1 * ms)

        counts = 0

        # ttl channel for toggling between what we send to the VCA: the zotino voltage or external modulation signal
        self.ttl12.output()
        self.ttl12.off()

        ############################
        # take a background image with the shims offset so the MOT doesn't form
        ############################
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT-0.5, self.AY_volts_MOT],
            channels=self.coil_channels)
        # delay(2 * ms)
        self.dds_cooling_DP.sw.on()

        # trigger the ThorCam with the Zotino
        self.zotino0.write_dac(6, 4.0)
        self.zotino0.load()
        delay(self.t_ThorCam_exposure_ms * ms + 1 * ms)
        self.zotino0.write_dac(6, 0.0)
        self.zotino0.load()

        delay(1*ms)

        # experiment loop
        for measurement in range(self.n_steps):

            # make sure any beams that can mess up the fW measurement (e.g. OP, excitation) are turned off
            if self.enable_laser_feedback:
                if measurement % 10 == 0:
                    self.laser_stabilizer.run()
                    delay(1 * ms)
                self.dds_FORT.sw.on()
                self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)

            self.ttl7.pulse(self.t_exp_trigger)  # in case we want to look at signals on an oscilloscope

            ############################
            # load the MOT
            ############################
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                channels=self.coil_channels)
            # delay(2 * ms)
            self.dds_cooling_DP.sw.on()

            # trigger the Luca to capture scattering from MOT5,6
            # then capture the MOT with ThorCam

            if measurement % self.image_frequency == 0:
                self.dds_AOM_A5.sw.off()
                delay(1 * ms)
                self.ttl_Luca_trigger.pulse(5 * ms)
                delay(self.t_Luca_exposure+100*ms)
                self.dds_AOM_A5.sw.on()

                # pulse off MOT6, take an image
                self.dds_AOM_A6.sw.off()
                delay(1 * ms)
                self.ttl_Luca_trigger.pulse(5 * ms)
                delay(self.t_Luca_exposure)
                self.dds_AOM_A6.sw.on()

            # wait for the MOT to load
            delay_mu(self.t_MOT_loading_mu)

            if measurement % self.image_frequency == 0:
                # trigger the ThorCam with the Zotino
                self.zotino0.write_dac(6, 4.0)
                self.zotino0.load()
                delay(self.t_ThorCam_exposure_ms * ms + 1*ms)
                self.zotino0.write_dac(6, 0.0)
                self.zotino0.load()

            # turn on the dipole trap and wait to load atoms
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
            delay_mu(self.t_FORT_loading_mu)

            # turn off the coils
            self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0], channels=self.coil_channels)

            delay(3 * ms)  # should wait for the MOT to dissipate

            # set the cooling DP AOM to the readout settings
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_MOT)

            ############################
            # take the first shot
            ############################
            self.dds_cooling_DP.sw.on()
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_first_shot)
            counts = self.ttl0.count(t_gate_end)
            delay(1 * ms)
            self.dds_cooling_DP.sw.off()

            delay(10 * ms)

            ############################
            # take the second shot
            ############################

            delay(self.t_SPCM_second_shot + 2*ms)

            # effectively turn the FORT AOM off
            self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)
            # set the cooling DP AOM to the MOT settings
            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

            # update the datasets
            self.append_to_dataset('photocounts', counts)
            self.append_to_dataset('photocounts_current_iteration', counts)

            # so we don't store too much data in artiq-master
            if measurement % self.last_photocounts_to_plot == 0:
                self.set_dataset("photocounts_current_iteration", [0], broadcast=True, archive=False)

        print("Experiment finished")
