"""
Simple-minded code for moving the MOT to optimize the atom loading rate.
"""

from artiq.experiment import *
import numpy as np

from utilities.BaseExperiment import BaseExperiment

class CoilScanChannel:
    """object for grouping a coil channel and the scan steps"""
    def __init__(self, Zotino_channel_index, V_steps_array, V_default):
        self.Zotino_channel_index = Zotino_channel_index
        self.V_steps_array = V_steps_array
        self.V_default = V_default


class MOTAutoPositioner(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # this should be close to the mean signal from the atom
        self.setattr_argument("trapped_atom_count_rate", NumberValue(10000))

        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        self.dt_signal_exposure = 10*ms
        self.dt_background_exposure = 500*ms

        self.sampler_buffer = np.zeros(8)

        self.dVz_bottom_steps_one_direction = np.linspace(0, 0.2, 10)
        self.dVZ_top_steps_one_direction = np.linspace(0, 0.2, 10)
        self.dVX_steps_one_direction = np.linspace(0,0.2,10)
        self.dVY_steps_one_direction = np.linspace(0,0.2,10)

        self.AZ_bottom_scan_ch = CoilScanChannel(self.AZ_bottom_Zotino_channel,
                                              self.dVZ_bottom_steps_one_direction,
                                              self.AZ_bottom_volts_MOT)
        self.AZ_top_scan_ch = CoilScanChannel(self.AZ_top_Zotino_channel,
                                              self.dVZ_top_steps_one_direction,
                                              self.AZ_top_volts_MOT)
        self.AX_scan_ch = CoilScanChannel(self.AX_Zotino_channel,
                                              self.dVX_steps_one_direction,
                                              self.AX_volts_MOT)
        self.AY_scan_ch = CoilScanChannel(self.AY_Zotino_channel,
                                              self.dVY_steps_one_direction,
                                              self.AY_volts_MOT)

        coil_scan_channels = [self.AZ_bottom_scan_ch, self.AZ_top_scan_ch, self.AX_scan_ch, self.AY_scan_ch]

        self.mean_counts = 0.0
        self.n_exposures = 10
        self.counts_list = [0.0]*self.n_exposures

        self.max_signal_search_attempts = 1

    @kernel
    def run(self):
        self.base.initialize_hardware()

        # load a MOT with the default settings
        self.dds_cooling_DP.sw.on()
        self.dds_AOM_A1.sw.on()
        self.dds_AOM_A2.sw.on()
        self.dds_AOM_A3.sw.on()
        self.dds_AOM_A4.sw.on()
        self.dds_AOM_A5.sw.on()
        self.dds_AOM_A6.sw.on()

        # delay for AOMs to thermalize
        delay(1*s)

        self.core.break_realtime()
        self.laser_stabilizer.run()
        delay(100 * ms)

        # offset the MOT enough to not have an atom signal
        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT - 0.3, self.AY_volts_MOT],
                             channels=self.coil_channels)

        # turn on the FORT
        self.dds_FORT.sw.on()
        delay(1*s) # delay for AOM to thermalize

        #  measure the counts from the FORT and background fluorescence and scattering
        t_end = self.ttl0.gate_rising(self.dt_background_exposure)
        background_counts = self.ttl0.count(t_end)
        background_counts_per_s = background_counts / self.dt_background_exposure
        delay(1*ms)

        # load a MOT with the default settings
        self.zotino0.set_dac(
            [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT - 0.3, self.AY_volts_MOT],
            channels=self.coil_channels)
        delay(self.t_MOT_loading)

        n_optimized_coils = 0

        """
        - while coils_for_which_we_latched < 4 and attempts < max_attempts:
            for each [coils that didn't latch], 
            
            advance the coil in one direction
            expose the SPCM several times sequentially and see if any of the counts
            is above the threshold. we can set the threshold high so that ambiguous events
            are ignored. if two events above the threshold are recorded for one coil step,
            we can "latch" and record the coil value as the entry boundary of the loading region.
            do this for a finite number of steps, breaking when the entry boundary is found.
            if the signal latches on the zeroth step, when the coil has not yet been advanced,
            flag this as an auto-latch and continue to the exit boundary loop.
            if entry boundary was not found, reset the coil, wait 1s to load the MOT, repeat the 
            loop above but with the coil advancing in the other direction.
            
            if an entry boundary was found, repeat the sequence of coil advance and serial exposures, 
            checking each time for at least one event that was above threshold. when a serial exposure 
            finds no counts above the threshold, record the coil value as the exit boundary of the
            loading region. break out of the loop when the exit boundary is found.
            
            if the signal auto-latched in the entry loop, reset the coil value, wait 1s, then repeat the
            loop above but with the coil advancing in the other direction. when the signal unlatches, 
            flag this coil value as the entry boundary.
            
            set the coil value to the mean of the entry and exit boundary values.
            else, flag this coil as not latched. the reason is that we may find the signal with the next
            coil, and then we should make sure to come back and optimize the signal with for this coil.
        """

        # this whole block can go inside a bigger loop to have repeated iterations
        n_attempts_to_find_signal = 0
        optimized_coil_channels = []
        while n_optimized_coils < 4 and n_attempts_to_find_signal < self.max_signal_search_attempts:

            # loop over the unoptimized coils
            for coil_ch in self.coil_scan_channels:

                


        self.dds_FORT.sw.off()

        print("Experiment finished.")


