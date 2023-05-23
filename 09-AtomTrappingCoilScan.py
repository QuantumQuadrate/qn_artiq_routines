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
9. Save the array of counts to a file


* Tested that all the delays and timings are functioning as expected by pulsing TTLs at different points and monitoring
    with an oscilloscope.

"""

from artiq.experiment import *
import csv
import math # for math
import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt
from subroutines.stabilizer import AOMPowerStabilizer


class SimpleAtomTrapping(EnvExperiment):

    def build(self):
        """
        declare hardware and user-configurable independent variables
        """

        # declare the hardwire devices/channels we will use
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul2_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")
        self.setattr_device("urukul2_ch0")
        self.setattr_device("urukul2_ch1")
        self.setattr_device("sampler0")  # for laser feedback
        self.setattr_device("zotino0")  # for controlling coils
        self.setattr_device("ttl0")  # input for counting SPCM clicks
        self.setattr_device("ttl7")  # output for experiment trigger

        # MOT coil settings
        self.setattr_argument("AZ_bottom_volts_MOT", NumberValue(0.92, unit="V", ndecimals=3, step=0.025),
                              "MOT coil settings")
        self.setattr_argument("AZ_top_volts_MOT", NumberValue(-2.76, unit="V", ndecimals=3, step=0.025),
                              "MOT coil settings")
        self.setattr_argument("AX_volts_MOT", NumberValue(-0.135, unit="V", ndecimals=3, step=0.025),
                              "MOT coil settings")
        self.setattr_argument("AY_volts_MOT", NumberValue(-0.021, unit="V", ndecimals=3, step=0.025),
                              "MOT coil settings")

        # dipole trap loading/PGC coil settings
        self.setattr_argument("AZ_bottom_volts_PGC", NumberValue(0, unit="V", ndecimals=3, step=0.025),
                              "PGC coil settings")
        self.setattr_argument("AZ_top_volts_PGC", NumberValue(0, unit="V", ndecimals=3, step=0.025),
                              "PGC coil settings")
        self.setattr_argument("AX_volts_PGC", NumberValue(-0.135, unit="V", ndecimals=3, step=0.025),
                              "PGC coil settings")
        self.setattr_argument("AY_volts_PGC", NumberValue(-0.021, unit="V", ndecimals=3, step=0.025),
                              "PGC coil settings")

        # imaging coil settings
        self.setattr_argument("AZ_bottom_volts_imaging", NumberValue(0, unit="V", ndecimals=3, step=0.025),
                              "Imaging coil settings")
        self.setattr_argument("AZ_top_volts_imaging", NumberValue(0, unit="V", ndecimals=3, step=0.025),
                              "Imaging coil settings")
        self.setattr_argument("AX_volts_imaging", NumberValue(-0.135, unit="V", ndecimals=3, step=0.025),
                              "Imaging coil settings")
        self.setattr_argument("AY_volts_imaging", NumberValue(-0.021, unit="V", ndecimals=3, step=0.025),
                              "Imaging coil settings")

        self.setattr_argument("f_FORT", NumberValue(210.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM1, FORT switching AOM")
        self.setattr_argument("p_FORT_loading", NumberValue(3, unit="dBm", scale=1, ndecimals=1),
                              "AOM1, FORT switching AOM")
        self.setattr_argument("p_FORT_imaging", NumberValue(3, unit="dBm", scale=1, ndecimals=1),
                              "AOM1, FORT switching AOM")
        self.setattr_argument("p_FORT_PGC", NumberValue(3, unit="dBm", scale=1, ndecimals=1),
                              "AOM1, FORT switching AOM")
        self.setattr_argument("f_780DP_MOT", NumberValue(115.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("f_780DP_PGC", NumberValue(115.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("f_780DP_imaging", NumberValue(115.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_780DP_MOT", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_780DP_PGC", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_780DP_imaging", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")

        self.setattr_argument("AOM3_freq", NumberValue(130.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM3, MOT cooling single pass")
        self.setattr_argument("AOM3_power", NumberValue(1, unit="dBm", scale=1, ndecimals=1),
                              "AOM3, MOT cooling single pass")

        self.setattr_argument("AOM4_freq", NumberValue(150.5 * MHz, unit="MHz", ndecimals=1), "AOM4, MOT RP/Exc")
        self.setattr_argument("AOM4_power", NumberValue(3, unit="dBm", scale=1, ndecimals=1), "AOM4, MOT RP/Exc")

        # the default power for the fiber AOMs was chosen to give roughly equal diffraction efficiency
        self.setattr_argument("AOM_A1_freq", NumberValue(78.51 * MHz, unit="MHz", ndecimals=2), "AOM A1")
        self.setattr_argument("AOM_A1_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A1")
        # self.setattr_argument("AOM_A1_ON", BooleanValue(default=False), "AOM A1")

        self.setattr_argument("AOM_A2_freq", NumberValue(78.48 * MHz, unit="MHz", ndecimals=2), "AOM A2")
        self.setattr_argument("AOM_A2_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A2")
        # self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "AOM A2")

        self.setattr_argument("AOM_A3_freq", NumberValue(78.49 * MHz, unit="MHz", ndecimals=2), "AOM A3")
        self.setattr_argument("AOM_A3_power", NumberValue(-3, unit="dBm", scale=1, ndecimals=1), "AOM A3")
        # self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "AOM A3")

        self.setattr_argument("AOM_A4_freq", NumberValue(78.5 * MHz, unit="MHz", ndecimals=2), "AOM A4")
        self.setattr_argument("AOM_A4_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A4")
        # self.setattr_argument("AOM_A4_ON", BooleanValue(default=False), "AOM A4")

        self.setattr_argument("AOM_A5_freq", NumberValue(78.47 * MHz, unit="MHz", ndecimals=2), "AOM A5")
        self.setattr_argument("AOM_A5_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A5")
        # self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "AOM A5")

        self.setattr_argument("AOM_A6_freq", NumberValue(78.52 * MHz, unit="MHz", ndecimals=2), "AOM A6")
        self.setattr_argument("AOM_A6_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A6")
        # self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "AOM A6")

        self.setattr_argument("n_measurements", NumberValue(10, ndecimals=0, step=1))
        self.setattr_argument("t_MOT_loading", NumberValue(350 * ms, unit="ms", ndecimals=0, step=10 * ms))
        self.setattr_argument("t_FORT_loading", NumberValue(100 * ms, unit="ms", ndecimals=1, step=10 * ms))
        self.setattr_argument("t_SPCM_exposure", NumberValue(50 * ms, unit="ms", ndecimals=1, step=5 * ms))

        self.setattr_argument("cooling_setpoint_mW", NumberValue(0.7), "Laser power stabilization")

        self.setattr_argument("datadir",
                              StringValue('C:\\Networking Experiment\\artiq codes\\artiq-master\\results\\'),"File to save data")
        self.setattr_argument("datafile", StringValue('atom_loading_counts.csv'),"File to save data")
        self.setattr_argument("prepend_date_to_datafile", BooleanValue(True),"File to save data")
        self.setattr_argument("print_measurement_number", BooleanValue(False), "Developer options")
        self.setattr_argument("print_meas_result", BooleanValue(False), "Developer options")
        self.setattr_argument("save_data", BooleanValue(True), "Developer options")

        self.setattr_argument("bins", NumberValue(50, ndecimals=0, step=1), "Histogram setup")
        self.setattr_argument("counts_per_bin", NumberValue(10, ndecimals=0, step=1), "Histogram setup")
        self.setattr_argument("print_counts", BooleanValue(True))

        print("build - done")

    def prepare(self):
        """
        performs initial calculations and sets parameter values before
        running the experiment. also sets data filename for now.

        any conversions from human-readable units to machine units (mu) are done here
        """

        # where to store the data
        self.t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
        if self.prepend_date_to_datafile:
            self.datafile = self.datadir + self.t_experiment_run + '_' + self.datafile
        else:
            self.datafile = self.datadir + self.datafile

        # experiment trigger pulse width
        self.t_exp_trigger = 1*ms

        # convert times to machine units
        self.t_MOT_loading_mu = self.core.seconds_to_mu(self.t_MOT_loading)
        self.t_FORT_loading_mu = self.core.seconds_to_mu(self.t_FORT_loading)
        self.t_SPCM_exposure_mu = self.core.seconds_to_mu(self.t_SPCM_exposure)

        # converts RF power in dBm to amplitudes in V
        self.ampl_FORT_loading = math.sqrt(2 * 50 * 10 ** (self.p_FORT_loading / 10 - 3))
        self.ampl_FORT_imaging = math.sqrt(2 * 50 * 10 ** (self.p_FORT_imaging / 10 - 3))
        self.ampl_FORT_PGC = math.sqrt(2 * 50 * 10 ** (self.p_FORT_PGC / 10 - 3))
        self.ampl_780DP_MOT = math.sqrt(2 * 50 * 10 ** (self.p_780DP_MOT / 10 - 3))
        self.ampl_780DP_PGC = math.sqrt(2 * 50 * 10 ** (self.p_780DP_PGC / 10 - 3))
        self.ampl_780DP_imaging = math.sqrt(2 * 50 * 10 ** (self.p_780DP_imaging / 10 - 3))
        self.AOM3_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM3_power / 10 - 3))
        self.AOM4_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM4_power / 10 - 3))

        self.AOM_A1_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A1_power / 10 - 3))
        self.AOM_A2_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A2_power / 10 - 3))
        self.AOM_A3_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A3_power / 10 - 3))
        self.AOM_A4_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A4_power / 10 - 3))
        self.AOM_A5_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A5_power / 10 - 3))
        self.AOM_A6_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A6_power / 10 - 3))

        # setup stabilization for the cooling laser power

        # todo: eventually read conversion functions such as this from a config file
        def volts_to_optical_mW(x: TFloat) -> TFloat:
            """
            the conversion of PD voltage to cooling light power at the switchyard MOT 1 path
            """
            x += 0.011  # this accounts for a mismatch between what the Sampler reads and what
            # the multimeter that I used for the fit reads
            return -0.195395 + 17.9214 * x

        self.sampler_buffer = [0.0] * 8
        self.cooling_volts_ch = 7
        self.AOMservo = AOMPowerStabilizer(experiment=self,
                                           dds_names=["urukul0_ch1"],
                                           sampler_name="sampler0",
                                           sampler_channels=[self.cooling_volts_ch],
                                           transfer_functions=[volts_to_optical_mW],
                                           setpoints=[self.cooling_setpoint_mW],  # in mW
                                           proportionals=[0.04],
                                           iters=5,  # keep iters/t_meas_delay small or rtio underflow
                                           t_meas_delay=20 * ms)

        self.coil_channels = [0, 1, 2, 3]

        self.hist_bins = np.zeros(self.bins, dtype=int)
        print("prepare - done")

    @kernel
    def run(self):
        self.init_hardware()
        # self.core.break_realtime()
        self.expt()
        print("Experiment finished.")

    @rpc(flags={"async"}) # means this code runs asynchronously; won't block the rtio counter
    def file_setup(self, rowheaders=[]):
        with open(self.datafile, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(rowheaders)
            f.close()

    @rpc(flags={"async"})
    def file_write(self, data):
        with open(self.datafile, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(data)
            f.close()

    # plot with the applet plot_hist instead
    # @rpc(flags={"async"})
    # def plot_data(self):
    #     """assumes one numeric datum per row"""
    #     with open(self.datafile, 'r', newline='') as f:
    #         reader = csv.reader(f)
    #         reader.__next__() # skip the header
    #         data = [int(row[0]) for row in reader]
    #         f.close()
    #     # xpts = range(len(data))
    #     plt.hist(data) # to-do set up some binning?
    #     plt.xlabel("Measurement index")
    #     plt.ylabel("Counts")
    #     plt.show()

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """
        self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
                             channels=self.coil_channels)

        self.set_dataset("photocounts", self.hist_bins, broadcast=True)
        # self.set_dataset("mot_photocounts", self.hist_bins, broadcast=True)

        self.file_setup(rowheaders=['counts'])

        # turn on cooling/RP AOMs
        self.urukul0_ch1.sw.on() # cooling double pass
        self.urukul0_ch2.sw.on()  # cooling single pass
        self.urukul0_ch3.sw.on()  # MOT repump

        delay(2000*ms) # wait for AOMS to thermalize in case they have been off.

        self.AOMservo.run()

        # loop the experiment sequence
        for measurement in range(self.n_measurements):

            self.ttl7.pulse(self.t_exp_trigger) # in case we want to look at signals on an oscilloscope

            # Set magnetic fields for MOT loading
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                channels=self.coil_channels)
            delay(1 * ms)  # avoid RTIOSequence error

            # Set and turn on fiber AOMs to load the MOT. The MOT AOMs upstream are assumed to be on.
            self.urukul1_ch0.sw.on()
            self.urukul1_ch1.sw.on()
            self.urukul1_ch2.sw.on()
            self.urukul1_ch3.sw.on()
            self.urukul2_ch0.sw.on()
            self.urukul2_ch1.sw.on()
            delay(1 * ms)

            # wait for the MOT to load
            delay_mu(self.t_MOT_loading_mu)

            # change the magnetic fields for loading the dipole trap
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_PGC, self.AZ_top_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                channels=self.coil_channels)

            # change double pass power and frequency to PGC settings
            # self.urukul1_ch1.set(frequency=self.f_780DP_PGC, amplitude=self.ampl_780DP_PGC)

            # turn on the dipole trap and wait to load atoms
            self.urukul0_ch0.sw.on()
            delay_mu(self.t_FORT_loading_mu)

            # change AOMs to "imaging" settings
            # self.urukul0_ch0.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_imaging)
            # self.urukul0_ch1.set(frequency=self.f_780DP_imaging, amplitude=self.ampl_780DP_imaging)

            # change the magnetic fields for imaging
            t_gate_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
            self.zotino0.set_dac(
                [self.AZ_bottom_volts_imaging, self.AZ_top_volts_imaging, self.AX_volts_imaging, self.AY_volts_imaging],
                channels=self.coil_channels)
            # delay(1*ms)

            # # take the shot
            # t_gate_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
            counts = self.ttl0.count(t_gate_end)
            if self.print_counts:
                print(counts)
            delay(10 * ms)

            # reset parameters
            # self.urukul1_ch0.sw.off()  # fiber AOMs off
            # self.urukul1_ch1.sw.off()
            # self.urukul1_ch2.sw.off()
            # self.urukul1_ch3.sw.off()
            # self.urukul2_ch0.sw.off()
            # self.urukul2_ch1.sw.off()
            self.urukul0_ch0.sw.off()  # FORT AOM off
            # self.urukul1_ch1.set(frequency=self.f_780DP_MOT, amplitude=self.ampl_780DP_MOT)
            # self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0],  # voltages must be floats or ARTIQ complains
            #                      channels=self.coil_channels)

            bin_idx = int(counts / self.counts_per_bin)
            if bin_idx < self.bins:
                self.hist_bins[bin_idx] += 1
                self.mutate_dataset("photocounts", bin_idx, self.hist_bins[bin_idx])

            if self.print_measurement_number:
                print("measurement", measurement)
            if self.print_meas_result:
                print("counts", counts)
            if self.save_data:
                self.file_write([counts])

        delay(1*ms)
        # leave MOT on at end of experiment, but turn off the FORT
        self.urukul0_ch0.sw.off()
        self.urukul1_ch0.sw.on()
        self.urukul1_ch1.sw.on()
        self.urukul1_ch2.sw.on()
        self.urukul1_ch3.sw.on()
        self.urukul2_ch0.sw.on()
        self.urukul2_ch1.sw.on()
        self.zotino0.set_dac([self.AZ_bottom_volts_MOT, self.AZ_top_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=self.coil_channels)

    @kernel
    def init_hardware(self):
        """
        Sets amplitudes and frequencies for the Urukul0 channels
        :return:
        """

        self.core.reset()
        self.ttl0.input()  # for reading pulses from SPCM
        self.ttl7.output()  # for outputting a trigger each cycle

        self.urukul0_cpld.init()
        self.urukul1_cpld.init()
        self.urukul2_cpld.init()
        self.urukul0_ch0.init()
        self.urukul0_ch1.init()
        self.urukul0_ch2.init()
        self.urukul0_ch3.init()
        self.urukul1_ch0.init()
        self.urukul1_ch1.init()
        self.urukul1_ch2.init()
        self.urukul1_ch3.init()
        self.urukul2_ch0.init()
        self.urukul2_ch1.init()

        self.urukul0_ch0.set_att(float(0))
        self.urukul0_ch1.set_att(float(0))
        self.urukul0_ch2.set_att(float(0))
        self.urukul0_ch3.set_att(float(0))
        self.urukul1_ch0.set_att(float(0))
        self.urukul1_ch1.set_att(float(0))
        self.urukul1_ch2.set_att(float(0))
        self.urukul1_ch3.set_att(float(0))
        self.urukul2_ch0.set_att(float(0))
        self.urukul2_ch1.set_att(float(0))
        self.zotino0.init()

        self.core.break_realtime()

        # URUKUL 0 - FORT, MOT and D2 state prep AOMs:
        delay(1*ms)
        self.urukul0_ch0.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_loading)
        delay(1*ms)
        self.urukul0_ch1.set(frequency=self.f_780DP_MOT, amplitude=self.ampl_780DP_MOT)
        delay(1*ms)
        self.urukul0_ch2.set(frequency=self.AOM3_freq, amplitude=self.AOM3_ampl)
        delay(1*ms)
        self.urukul0_ch3.set(frequency=self.AOM4_freq, amplitude=self.AOM4_ampl)

        # URUKUL 1 - MOT arm fiber AOMs:
        delay(1*ms)
        self.urukul1_ch0.set(frequency=self.AOM_A2_freq, amplitude=self.AOM_A2_ampl)
        self.urukul1_ch1.set(frequency=self.AOM_A3_freq, amplitude=self.AOM_A3_ampl)
        self.urukul1_ch2.set(frequency=self.AOM_A1_freq, amplitude=self.AOM_A1_ampl)
        self.urukul1_ch3.set(frequency=self.AOM_A6_freq, amplitude=self.AOM_A6_ampl)
        self.urukul2_ch0.set(frequency=self.AOM_A4_freq, amplitude=self.AOM_A4_ampl)
        self.urukul2_ch1.set(frequency=self.AOM_A5_freq, amplitude=self.AOM_A5_ampl)

        self.AOMservo.get_dds_settings()  # must come after relevant DDS's have been set
        delay(100*ms)

        print("init_hardware - done")