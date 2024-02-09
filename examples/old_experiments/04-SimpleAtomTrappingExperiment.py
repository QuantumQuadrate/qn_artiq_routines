"""
A simple experiment to look for signs of trapped atoms

1. Turn on cooling and RP AOMs
Experiment cycle (repeats n times)
2. Turn on magnetic fields (Zotino)
3. Turn on fiber AOMs (Urukul), wait some time to load the MOT
4. Turn on the dipole trap AOM
5. Turn off the magnetic fields (PGC phase)
6. Turn off the fiber AOMs
7. Turn on the fiber AOMs and read from a single SPCM (TTL) for a certain exposure time
8. Store the number of counts registered by the SPCM in an array
End of experiment
9. Save the array of counts to a file

To-do and notes:
 [*] ttl output on channel 4:
    [*] find out how to configure TTL card channel for output
    [*] initialize output to low
 [] zotino for coil control
    [*] check if Akbar had an impedance conversion
    [*] c/p his code for setup
    [] convert voltages to machine units? this is what Garrett was doing
    [] get Zotino calibration code from Kais in the future
 [*] urukul for AOMs
    [*] update init_hardware with the MOT frequency/power for the cooling double pass AOM


"""

from artiq.experiment import *
import csv
from artiq.coredevice import ad53xx # for converting volts to mu for the zotino
from artiq.coredevice.exceptions import RTIOUnderflow
import math # for math
import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt

class SimpleAtomTrapping(EnvExperiment):


    def build(self):
        """
        declare hardware and user-configurable independent variables
        """

        # declare the hardwire devices/channels we will use
        self.setattr_device("core")
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul0_ch0")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")
        self.setattr_device("urukul1_ch0")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")
        self.setattr_device("zotino0") # for controlling coils
        self.setattr_device("ttl0")  # input for counting SPCM clicks
        self.setattr_device("ttl7")  # output for experiment trigger

        self.setattr_argument("AZ_top_volts_MOT", NumberValue(1.8*(-1.64/2.5), unit="V", ndecimals=3, step=0.025), "A-Z shim/quad top coil")
        self.setattr_argument("AZ_top_volts_PGC", NumberValue(0 * (-1.64 / 2.5), unit="V", ndecimals=3, step=0.025), "A-Z shim/quad top coil")
        self.setattr_argument("AZ_top_volts_RO", NumberValue(0 * (-1.64 / 2.5), unit="V", ndecimals=3, step=0.025), "A-Z shim/quad top coil")
        self.setattr_argument("AZ_bottom_volts_MOT", NumberValue(-1.75*(-1.64/2.5), unit="V", ndecimals=3, step=0.025), "A-Z shim/quad bottom coils")
        self.setattr_argument("AZ_bottom_volts_PGC", NumberValue(0*(-1.64/2.5), unit="V", ndecimals=3, step=0.025), "A-Z shim/quad bottom coils")
        self.setattr_argument("AZ_bottom_volts_RO", NumberValue(0*(-1.64/2.5), unit="V", ndecimals=3, step=0.025), "A-Z shim/quad bottom coils")
        self.setattr_argument("AX_volts_MOT", NumberValue(0.5*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-X shim coils")
        self.setattr_argument("AX_volts_PGC", NumberValue(0.5*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-X shim coils")
        self.setattr_argument("AX_volts_RO", NumberValue(0.5*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-X shim coils")
        self.setattr_argument("AY_volts_MOT", NumberValue(0.4*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-Y shim coils")
        self.setattr_argument("AY_volts_PGC", NumberValue(0.4*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-Y shim coils")
        self.setattr_argument("AY_volts_RO", NumberValue(0.4*(-1.64/2.5), unit="V", ndecimals=3, step=0.025),"A-Y shim coils")

        # define user-configurable independent variables in human-readable units
        # these will show up in the GUI
        self.setattr_argument("f_FORT", NumberValue(100.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM1, FORT switching AOM")
        self.setattr_argument("p_FORT_loading", NumberValue(0, unit="dBm", scale=1, ndecimals=1),
                              "AOM1, FORT switching AOM")
        self.setattr_argument("p_FORT_RO", NumberValue(0, unit="dBm", scale=1, ndecimals=1),
                              "AOM1, FORT switching AOM")
        self.setattr_argument("p_FORT_PGC", NumberValue(0, unit="dBm", scale=1, ndecimals=1),
                              "AOM1, FORT switching AOM")
        self.setattr_argument("f_cooling_DP_MOT", NumberValue(115.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("f_cooling_DP_PGC", NumberValue(115.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("f_cooling_DP_RO", NumberValue(115.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_cooling_DP_MOT", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_cooling_DP_PGC", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")
        self.setattr_argument("p_cooling_DP_RO", NumberValue(-0.2, unit="dBm", scale=1, ndecimals=1),
                              "AOM2, MOT cooling double pass")
        # self.setattr_argument("Cooling_DP_AOM_ON", BooleanValue(default=False), "AOM2, MOT cooling double pass")

        self.setattr_argument("f_D1_pumping_SP", NumberValue(130.0 * MHz, unit="MHz", ndecimals=1),
                              "AOM3, MOT cooling single pass")
        self.setattr_argument("p_D1_pumping_SP", NumberValue(1, unit="dBm", scale=1, ndecimals=1),
                              "AOM3, MOT cooling single pass")
        # self.setattr_argument("D1_pumping_SP_AOM_ON", BooleanValue(default=False), "AOM3, MOT cooling single pass")

        self.setattr_argument("f_pumping_repump", NumberValue(150.5 * MHz, unit="MHz", ndecimals=1), "AOM4, MOT RP/Exc")
        self.setattr_argument("p_pumping_repump", NumberValue(3, unit="dBm", scale=1, ndecimals=1), "AOM4, MOT RP/Exc")
        # self.setattr_argument("pumping_repump_AOM_ON", BooleanValue(default=False), "AOM4, MOT RP/Exc")

        # the default power for the fiber AOMs was chosen to give roughly equal diffraction efficiency, empirically
        self.setattr_argument("AOM_A2_freq", NumberValue(78.48 * MHz, unit="MHz", ndecimals=2), "AOM A2")
        self.setattr_argument("AOM_A2_power", NumberValue(-5, unit="dBm", scale=1, ndecimals=1), "AOM A2")
        # self.setattr_argument("AOM_A2_ON", BooleanValue(default=False), "AOM A2")

        self.setattr_argument("AOM_A3_freq", NumberValue(78.49 * MHz, unit="MHz", ndecimals=2), "AOM A3")
        self.setattr_argument("AOM_A3_power", NumberValue(-3, unit="dBm", scale=1, ndecimals=1), "AOM A3")
        # self.setattr_argument("AOM_A3_ON", BooleanValue(default=False), "AOM A3")

        self.setattr_argument("AOM_A5_freq", NumberValue(78.5 * MHz, unit="MHz", ndecimals=2), "AOM A5")
        self.setattr_argument("AOM_A5_power", NumberValue(-3, unit="dBm", scale=1, ndecimals=1), "AOM A5")
        # self.setattr_argument("AOM_A5_ON", BooleanValue(default=False), "AOM A5")

        self.setattr_argument("AOM_A6_freq", NumberValue(78.51 * MHz, unit="MHz", ndecimals=2), "AOM A6")
        self.setattr_argument("AOM_A6_power", NumberValue(0, unit="dBm", scale=1, ndecimals=1), "AOM A6")
        # self.setattr_argument("AOM_A6_ON", BooleanValue(default=False), "AOM A6")

        self.setattr_argument("n_measurements", NumberValue(10, ndecimals=0, step=1))
        self.setattr_argument("t_MOT_loading", NumberValue(350 * ms, unit="ms", ndecimals=0, step=10 * ms))
        self.setattr_argument("t_FORT_loading", NumberValue(100 * ms, unit="ms", ndecimals=1, step=10 * ms))
        self.setattr_argument("t_SPCM_exposure", NumberValue(50 * ms, unit="ms", ndecimals=1, step=5 * ms))
        self.setattr_argument("datadir", StringValue('/home/fiber-cavity/NetworkExperimentData/'),"File to save data")
        self.setattr_argument("datafile", StringValue('counts_test.csv'),"File to save data")
        self.setattr_argument("prepend_date_to_datafile", BooleanValue(True),"File to save data")
        self.setattr_argument("print_measurement_number", BooleanValue(False), "Developer options")
        self.setattr_argument("print_meas_result", BooleanValue(False), "Developer options")
        self.setattr_argument("save_data", BooleanValue(True), "Developer options")

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
        # self.AOM1_ampl = math.sqrt(2 * 50 * 10 ** (self.p_FORT_loading / 10 - 3))
        # self.AOM2_ampl = math.sqrt(2 * 50 * 10 ** (self.p_cooling_DP_MOT / 10 - 3))
        self.stabilizer_FORT.amplitude = math.sqrt(2 * 50 * 10 ** (self.p_FORT_loading / 10 - 3))
        self.ampl_FORT_RO = math.sqrt(2 * 50 * 10 ** (self.p_FORT_RO / 10 - 3))
        self.ampl_FORT_PGC = math.sqrt(2 * 50 * 10 ** (self.p_FORT_PGC / 10 - 3))
        self.ampl_cooling_DP_MOT = math.sqrt(2 * 50 * 10 ** (self.p_cooling_DP_MOT / 10 - 3))
        self.ampl_cooling_DP_PGC = math.sqrt(2 * 50 * 10 ** (self.p_cooling_DP_PGC / 10 - 3))
        self.ampl_cooling_DP_RO = math.sqrt(2 * 50 * 10 ** (self.p_cooling_DP_RO / 10 - 3))
        self.AOM3_ampl = math.sqrt(2 * 50 * 10 ** (self.p_D1_pumping_SP / 10 - 3))
        self.AOM4_ampl = math.sqrt(2 * 50 * 10 ** (self.p_pumping_repump / 10 - 3))

        self.AOM_A2_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A2_power / 10 - 3))
        self.AOM_A3_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A3_power / 10 - 3))
        self.AOM_A5_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A5_power / 10 - 3))
        self.AOM_A6_ampl = math.sqrt(2 * 50 * 10 ** (self.AOM_A6_power / 10 - 3))

        # convert frequencies to frequency-tuning-words... when is this necessary?

        # convert voltage for coils to machine units?
        # self.zero_volts = NumberValue(0, unit="V") #ad53xx.voltage_to_mu(0)
        # print(self.zero_volts)

    @kernel
    def run(self):
        self.init_hardware()
        self.expt()
        print("Experiment finished.")

    @kernel
    def mot_and_shot(self) -> TInt32: # this ARTIQ type declaration means the function returns a float
        """
        Load a MOT, then dipole trap, then take a single shot

        In more complicated experiment, this function should be decomposed into
        the constituent experiment phases: load the MOT, etc.

        :return:
        """

        # Turn on the magnetic fields
        self.zotino0.set_dac([self.AZ_top_volts_MOT, self.AZ_bottom_volts_MOT, self.AX_volts_MOT, self.AY_volts_MOT],
                             channels=[0, 1, 2, 3])
        delay(1*ms) # avoid RTIOSequence error

        # Set and turn on fiber AOMs to load the MOT. The MOT AOMs upstream are assumed to be on.
        self.urukul1_ch0.sw.on()
        self.urukul1_ch1.sw.on()
        self.urukul1_ch2.sw.on()
        self.urukul1_ch3.sw.on()
        delay(1*ms)

        # wait for the MOT to load
        delay_mu(self.t_MOT_loading_mu)

        # change the magnetic fields for loading the dipole trap
        self.zotino0.set_dac([self.AZ_top_volts_PGC, self.AZ_bottom_volts_PGC, self.AX_volts_PGC, self.AY_volts_PGC],
                             channels=[0, 1, 2, 3])

        # change double pass power and frequency to PGC settings
        self.urukul1_ch1.set(frequency=self.f_cooling_DP_PGC, amplitude=self.ampl_cooling_DP_PGC)

        # turn on the dipole trap and wait to load atoms
        # self.urukul0_ch0.sw.on() # todo: this channel is currently connected to the pumping RP AOM.
        delay_mu(self.t_FORT_loading_mu)

        # change AOMs to imaging settings
        self.urukul0_ch0.set(frequency=self.f_FORT, amplitude=self.ampl_FORT_RO)
        self.urukul0_ch1.set(frequency=self.f_cooling_DP_RO, amplitude=self.ampl_cooling_DP_RO)


        # change the magnetic fields for imaging
        self.zotino0.set_dac([self.AZ_top_volts_RO, self.AZ_bottom_volts_RO, self.AX_volts_RO, self.AY_volts_RO],
                             channels=[0, 1, 2, 3])

        # take the shot
        t_gate_end = self.ttl0.gate_rising(self.t_SPCM_exposure)
        counts = self.ttl0.count(t_gate_end)
        # print(counts) # To-do print value to file instead. is it possible to have a return?
        delay(10*ms)

        # reset parameters
        self.urukul1_ch0.sw.off() # fiber AOMs off
        self.urukul1_ch1.sw.off()
        self.urukul1_ch2.sw.off()
        self.urukul1_ch3.sw.off()
        self.urukul0_ch0.sw.off() # FORT AOM off
        self.urukul1_ch1.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)
        self.zotino0.set_dac([0.0, 0.0, 0.0, 0.0], # voltages must be floats or ARTIQ complains
                             channels=[0, 1, 2, 3])

        return counts

    @rpc(flags={"async"}) # means this code runs asynchronously; won't block the rtio counter
    def file_setup(self,rowheaders=[]):
        with open(self.datafile,'w',newline='') as f:
            writer = csv.writer(f)
            writer.writerow(rowheaders)
            f.close()

    @rpc(flags={"async"})
    def file_write(self, data):
        with open(self.datafile, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(data)
            f.close()

    @rpc(flags={"async"})
    def plot_data(self):
        """assumes one numeric datum per row"""
        with open(self.datafile, 'r', newline='') as f:
            reader = csv.reader(f)
            reader.__next__() # skip the header
            data = [int(row[0]) for row in reader]
            f.close()
        # xpts = range(len(data))
        plt.hist(data) # to-do set up some binning?
        plt.xlabel("Measurement index")
        plt.ylabel("Counts")
        plt.show()

    @kernel
    def expt(self):
        """
        The experiment loop.

        :return:
        """

        # file setup
        # datafile = '/home/fiber-cavity/NetworkExperimentData/counts_test.csv'
        self.file_setup(rowheaders=['counts'])

        # turn on cooling/RP AOMs
        self.urukul0_ch1.sw.on() # cooling double pass
        self.urukul0_ch2.sw.on()  # cooling single pass
        self.urukul0_ch3.sw.on()  # MOT repump

        delay(2000*ms) # wait for AOMS to thermalize in case they have been off.

        for measurement in range(self.n_measurements):

            # do the experiment sequence
            self.ttl7.pulse(self.t_exp_trigger)
            data = self.mot_and_shot()
            if self.print_measurement_number:
                print("measurement", measurement)
            if self.print_meas_result:
                print("counts", data)
            if self.save_data:
                self.file_write([data])

        self.plot_data()


    @kernel
    def init_hardware(self):
        """
        Sets amplitudes and frequencies for the Urukul0 channels
        :return:
        """

        self.core.reset()
        self.ttl0.input()  # for reading pulses from SPCM
        self.ttl7.output()  # for outputting a trigger each cycle

        self.urukul0_ch0.init()
        self.urukul0_ch1.init()
        self.urukul0_ch2.init()
        self.urukul0_ch3.init()
        self.urukul1_ch0.init()
        self.urukul1_ch1.init()
        self.urukul1_ch2.init()
        self.urukul1_ch3.init()

        self.urukul0_ch0.set_att(float(0))
        self.urukul0_ch1.set_att(float(0))
        self.urukul0_ch2.set_att(float(0))
        self.urukul0_ch3.set_att(float(0))
        self.urukul1_ch0.set_att(float(0))
        self.urukul1_ch1.set_att(float(0))
        self.urukul1_ch2.set_att(float(0))
        self.urukul1_ch3.set_att(float(0))
        self.zotino0.init()
        self.urukul0_cpld.init()
        self.urukul1_cpld.init()

        self.core.break_realtime()

        self.urukul0_cpld.set_profile(0)
        self.urukul1_cpld.set_profile(0)

        self.core.break_realtime()

        # URUKUL 0 - FORT, MOT and D2 state prep AOMs:
        delay(1 * ms)
        self.urukul0_ch0.set(frequency=self.f_FORT, amplitude=self.stabilizer_FORT.amplitude)

        delay(1 * ms)
        self.urukul0_ch1.set(frequency=self.f_cooling_DP_MOT, amplitude=self.ampl_cooling_DP_MOT)

        delay(1 * ms)
        self.urukul0_ch2.set(frequency=self.f_D1_pumping_SP, amplitude=self.AOM3_ampl)

        delay(1 * ms)
        self.urukul0_ch3.set(frequency=self.f_pumping_repump, amplitude=self.AOM4_ampl)

        # URUKUL 1 - MOT arm fiber AOMs:
        delay(1 * ms)
        self.urukul1_ch0.set(frequency=self.AOM_A2_freq, amplitude=self.AOM_A2_ampl)

        delay(1 * ms)
        self.urukul1_ch1.set(frequency=self.AOM_A3_freq, amplitude=self.AOM_A3_ampl)

        delay(1 * ms)
        self.urukul1_ch2.set(frequency=self.AOM_A6_freq, amplitude=self.AOM_A6_ampl)

        delay(1 * ms)
        self.urukul1_ch3.set(frequency=self.AOM_A5_freq, amplitude=self.AOM_A5_ampl)