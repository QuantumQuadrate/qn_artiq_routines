"""
Simple-minded code for moving the MOT to optimize the atom loading rate.
"""

from artiq.experiment import *
import numpy as np
import sys
sys.path.append('C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\')


from utilities.BaseExperiment import BaseExperiment

class CoilScanChannel:
    """object for grouping a coil channel and the scan steps"""
    def __init__(self, name, Zotino_ch_index, dV_steps_array_one_direction, V_default):
        self.name = name
        self.Zotino_ch_index = Zotino_ch_index
        self.dV_steps_array_one_direction = dV_steps_array_one_direction
        self.dV = dV_steps_array_one_direction[1] - dV_steps_array_one_direction[0]
        self.V_default = V_default
        self.V_steps_array_forward = dV_steps_array_one_direction + V_default
        self.V_steps_array_backward = -1 * dV_steps_array_one_direction + V_default
        self.optimized=False
        self.V_signal_boundary1 = 0.0
        self.V_signal_boundary2 = 0.0

class MOTAutoPositioner(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """
        self.base = BaseExperiment(experiment=self)
        self.base.build()

        # this should be close to the mean signal from the atom
        self.setattr_argument("trapped_atom_counts_per_s", NumberValue(10000))
        self.base.set_datasets_from_gui_args()
        print("build - done")

    def prepare(self):
        self.base.prepare()

        self.dt_signal_exposure = 10*ms
        self.dt_background_exposure = 500*ms

        self.sampler_buffer = np.zeros(8)

        # the zeroth step should be zero so that we first check the counts with the current coil settings
        self.dVZ_bottom_steps_one_direction = np.linspace(0, 0.1, 40)
        self.dVZ_top_steps_one_direction = np.linspace(0, 0.1, 40)
        self.dVX_steps_one_direction = np.linspace(0,0.1,40)
        self.dVY_steps_one_direction = np.linspace(0,0.1,40)

        self.AZ_bottom_scan_ch = CoilScanChannel(
            "AZ_bottom",
            self.AZ_bottom_Zotino_channel,
            self.dVZ_bottom_steps_one_direction,
            self.AZ_bottom_volts_MOT)
        self.AZ_top_scan_ch = CoilScanChannel(
            "AZ_top",
            self.AZ_top_Zotino_channel,
            self.dVZ_top_steps_one_direction,
            self.AZ_top_volts_MOT)
        self.AX_scan_ch = CoilScanChannel(
            "AX",
            self.AX_Zotino_channel,
            self.dVX_steps_one_direction,
            self.AX_volts_MOT)
        self.AY_scan_ch = CoilScanChannel(
            "AY",
            self.AY_Zotino_channel,
            self.dVY_steps_one_direction,
            self.AY_volts_MOT)

        self.coil_scan_channels = [self.AZ_bottom_scan_ch,
                                   self.AZ_top_scan_ch,
                                   self.AX_scan_ch,
                                   self.AY_scan_ch]

        self.mean_counts = 0.0
        self.n_exposures = 10
        self.counts_list = [0.0]*self.n_exposures

        self.max_signal_search_attempts = 1

        self.set_dataset(self.count_rate_dataset,
                         [0.0],
                         broadcast=True)

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

        self.dds_FORT.sw.on()
        delay(1*s) # delay for AOM to thermalize

        #  measure the counts from the FORT and background fluorescence and scattering
        t_end = self.ttl0.gate_rising(self.dt_background_exposure)
        background_counts = self.ttl0.count(t_end)
        background_counts_per_s = background_counts / self.dt_background_exposure
        counts_per_s_threshold = self.trapped_atom_counts_per_s + background_counts_per_s
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

        # while n_optimized_coils < 4: #and n_attempts_to_find_signal < self.max_signal_search_attempts:

        # loop over the unoptimized coils.
        for coil_ch in self.coil_scan_channels:
            if not coil_ch.optimized:

                self.print_async("optimizing", coil_ch.name)

                self.laser_stabilizer.run()
                delay(1*ms)
                self.dds_FORT.sw.on()

                boundary1_latched = False
                boundary2_latched = False
                needs_reverse_scan = False
                for V in coil_ch.V_steps_array_forward:

                    self.zotino0.write_dac(coil_ch.Zotino_ch_index, V)
                    self.zotino0.load()

                    # expose the SPCM
                    there_is_still_a_signal = False
                    for i in range(self.n_exposures):
                        t_end = self.ttl0.gate_rising(self.dt_signal_exposure)
                        counts_per_s = self.ttl0.count(t_end)/self.dt_signal_exposure
                        delay(1 * ms)
                        self.append_to_dataset(self.count_rate_dataset, counts_per_s)

                        if counts_per_s > counts_per_s_threshold:
                            if not boundary1_latched:
                                # if we find any atom signal, we count that as having found the signal
                                coil_ch.V_signal_boundary1 = V
                                boundary1_latched = True
                                self.print_async("found the single atom signal!")

                                # technically, boundary 1 may be "behind" us, so we still need to scan
                                # in the other direction later, so we flag this
                                if V == coil_ch.V_default:
                                    needs_reverse_scan = True
                            else:
                                there_is_still_a_signal = True

                            # either we've just found boundary1, or we haven't reached boundary2, so we break
                            break

                    # if there was no signal found during that exposure loop, we've found boundary2
                    if boundary1_latched and not there_is_still_a_signal:
                        boundary2_latched = True
                        coil_ch.V_signal_boundary2 = V
                        V_signal_center = (coil_ch.V_signal_boundary2 + coil_ch.V_signal_boundary1)/2
                        self.print_async(coil_ch.name,"set to the center of the signal range")
                        self.zotino0.write_dac(coil_ch.Zotino_ch_index, V_signal_center)
                        self.zotino0.load()
                        break

                if needs_reverse_scan:
                    boundary1_latched = False

                    self.zotino0.write_dac(coil_ch.Zotino_ch_index, coil_ch.V_default)
                    self.zotino0.load()
                    delay(self.t_MOT_loading)

                    for V in coil_ch.V_steps_array_backward:

                        self.zotino0.write_dac(coil_ch.Zotino_ch_index, V)
                        self.zotino0.load()

                        # expose the SPCM
                        there_is_still_a_signal = False
                        for i in range(self.n_exposures):
                            t_end = self.ttl0.gate_rising(self.dt_signal_exposure)
                            counts_per_s = self.ttl0.count(t_end) / self.dt_signal_exposure
                            delay(1*ms)
                            self.append_to_dataset(self.count_rate_dataset, counts_per_s)

                            if counts_per_s > counts_per_s_threshold:
                                there_is_still_a_signal = True
                                break

                        if not there_is_still_a_signal:
                            boundary1_latched = True
                            break

                    if boundary1_latched and boundary2_latched:
                        coil_ch.optimized = True
                        n_optimized_coils += 1
                        V_signal_center = (coil_ch.V_signal_boundary2 + coil_ch.V_signal_boundary1) / 2
                        self.print_async(coil_ch.name, "set to the center of the signal range")
                        self.zotino0.write_dac(coil_ch.Zotino_ch_index, V_signal_center)
                        self.zotino0.load()


        self.dds_FORT.sw.off()

        print("Experiment finished.")


