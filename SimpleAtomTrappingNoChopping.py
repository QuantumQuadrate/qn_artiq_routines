"""
A simple experiment to look for a trapped single atom signal

1. Turn on cooling and RP AOMs
Experiment cycle (repeats n times)
2. Turn on magnetic fields (Zotino)
3. Turn on fiber AOMs (Urukul), wait some time to load the MOT
4. Turn on the dipole trap AOM
5. Turn off the quadrupole fields (PGC phase)
6. Turn off the fiber AOMs
7. Turn on the fiber AOMs and read from a single SPCM (TTL) for a certain exposure time
8. Store the number of counts registered by the SPCM in an array
End of experiment
9. Save the array of SPCM0_RO1 to a file


* Tested that all the delays and timings are functioning as expected by pulsing TTLs at different points and monitoring
    with an oscilloscope.

"""

from artiq.experiment import *
import csv
import numpy as np
from datetime import datetime as dt

from utilities.BaseExperiment import BaseExperiment


class SimpleAtomTrapNoChop(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        self.setattr_argument("n_measurements", NumberValue(10, ndecimals=0, step=1))
        self.setattr_argument("print_measurement_number", BooleanValue(False))
        self.setattr_argument("dt_exposure", NumberValue(10 * ms, unit='ms'))
        self.setattr_argument("iteration_delay", NumberValue(0 * s, unit='s'))
        self.setattr_argument("drop_MOT", BooleanValue(False))
        ### When selected, the MOT is dropped after atom loading.

        self.setattr_argument("bins", NumberValue(50, ndecimals=0, step=1), "Histogram setup (set bins=0 for auto)")
        self.setattr_argument("print_counts", BooleanValue(False))
        self.setattr_argument("enable_laser_feedback", BooleanValue(default=True),"Laser power stabilization")
        self.setattr_argument("control_experiment", BooleanValue(False))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """
        self.base.prepare()
        self.t_SPCM_exposure = self.dt_exposure # overwrite the exposure time for this experiment

        self.t_exp_trigger = 1*ms

        self.sampler_buffer = np.full(8, 0.0)
        self.cooling_volts_ch = 7

        # self.hist_bins = np.zeros(self.bins, dtype=int)
        # self.SPCM0_RO1 = np.full(self.n_measurements, 0.0)

        print("prepare - done")

    @kernel
    def run(self):
        self.base.initialize_hardware()
        self.expt()
        print("******************** Experiment finished *****************")

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """
        self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
                             channels=self.coil_channels)

        self.set_dataset("SPCM0_counts_per_s", [0], broadcast=True)
        self.set_dataset("photocount_bins", [self.bins], broadcast=True)

        # turn on cooling/RP AOMs
        self.dds_cooling_DP.sw.on() # cooling double pass

        delay(2000*ms) # wait for AOMS to thermalize in case they have been off.

        if self.enable_laser_feedback:
            self.laser_stabilizer.run()

        # loop the experiment sequence
        for measurement in range(self.n_measurements):

            if self.enable_laser_feedback:
                if measurement % 10 == 0:
                    self.laser_stabilizer.run()
                    delay(1 * ms)
                self.dds_FORT.sw.on()
                self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)

            # self.ttl7.pulse(self.t_exp_trigger) # in case we want to look at signals on an oscilloscope

            if self.control_experiment and measurement % 2 == 0:
                # change the Y MOT coil enough to lose the MOT but not enough to significantly change the fluorescence
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT-0.8],
                    channels=self.coil_channels)
            else:
                # Set magnetic fields for MOT loading
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                    channels=self.coil_channels)
            delay(1 * ms)  # avoid RTIOSequence error

            # Set and turn on fiber AOMs to load the MOT. The MOT AOMs upstream are assumed to be on.
            self.dds_AOM_A2.sw.on()
            self.dds_AOM_A3.sw.on()
            self.dds_AOM_A1.sw.on()
            self.dds_AOM_A6.sw.on()
            self.dds_AOM_A4.sw.on()
            self.dds_AOM_A5.sw.on()
            delay(1 * ms)

            # ### turn on the dipole trap with the MOT beams
            # self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

            ### wait for the MOT to load
            delay_mu(self.t_MOT_loading_mu)

            ### turn on the dipole trap and wait to load atoms
            # self.dds_FORT.sw.on()
            self.dds_FORT.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)
            delay_mu(self.t_FORT_loading_mu)

            if self.drop_MOT:
                ### Turn off the coils
                self.zotino0.set_dac(
                    [self.AZ_bottom_volts_PGC, -self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                    channels=self.coil_channels)
                delay(1 * ms)

                ### turning off the MOT beams
                self.dds_AOM_A2.sw.off()
                self.dds_AOM_A3.sw.off()
                self.dds_AOM_A1.sw.off()
                self.dds_AOM_A6.sw.off()
                self.dds_AOM_A4.sw.off()
                self.dds_AOM_A5.sw.off()
                delay(30 * ms)

                ### turning on the MOT beams for imaging
                # self.dds_cooling_DP.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)
                self.dds_AOM_A2.sw.on()
                self.dds_AOM_A3.sw.on()
                self.dds_AOM_A1.sw.on()
                self.dds_AOM_A6.sw.on()
                self.dds_AOM_A4.sw.on()
                self.dds_AOM_A5.sw.on()
                # delay(1 * ms)

            ### take the shot
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
            SPCM0_RO1 = self.ttl0.count(t_gate_end)
            delay(1*ms)
            if self.print_counts:
                print(SPCM0_RO1)
            delay(10 * ms)

            self.ttl6.pulse(5*ms)
            delay(50*ms)

            # reset AOMs and coils
            # self.dds_AOM_A2.sw.off()  # fiber AOMs off
            # self.dds_AOM_A3.sw.off()
            # self.dds_AOM_A1.sw.off()
            # self.dds_AOM_A6.sw.off()
            # self.dds_AOM_A4.sw.off()
            # self.dds_AOM_A5.sw.off()
            # self.dds_FORT.sw.off()  # FORT AOM off

            self.dds_cooling_DP.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

            self.append_to_dataset('SPCM0_counts_per_s', SPCM0_RO1/self.t_SPCM_exposure)

            if self.print_measurement_number:
                print("measurement", measurement)
            delay(10*ms)

            delay(self.iteration_delay)

        delay(1*ms)
        # leave MOT on at end of experiment, but turn off the FORT
        self.dds_FORT.sw.off()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A6.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)
        self.dds_FORT.set(frequency=self.f_FORT - 30 * MHz, amplitude=self.stabilizer_FORT.amplitude)