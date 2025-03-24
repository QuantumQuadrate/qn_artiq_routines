from artiq.experiment import *
import numpy as np
from collections import namedtuple


def setattr_variables(experiment, exclude_list=[], exclude_keywords=[]):
    """
    Set variables in your experiment from existing datasets of the same names.

    Usage: in the build method of your experiment, declare a list of variables
    then call this method (having imported it at the top) with experiment=self.

    :param experiment: an ARTIQ Experiment object with an attribute variables which is a list
        of strings which are the variables names we want. It is assumed that the variables
        have already been defined in ExperimentVariables.
    :param exclude_list list(str): for variables that store very long lists of data, we don't we want to try to add
        these. they can lead to timeout errors are also generally not needed as variables within an experiment.
    :param exclude_keywords list(str): ignore variable names that contain keywords in this list. this is useful for
        auto-generated datasets such as the RF history datasets created in aom_feedback.py, so we don't have to
        remember to update the exclude_list for those when we create or remove feedback channels.
        For variables that store very long lists of data, we don't we want to try to add
        these. they can lead to timeout errors are also generally not needed as variables within an experiment.
    :return:
    """
    for var in experiment.variables:
        if not var in exclude_list and not np.product([x in var for x in exclude_keywords]):
            try:
                value = experiment.get_dataset(var)
                setattr(experiment, var, value)
            except Exception as e:
                if type(e) == KeyError:  # if the variable does not exist
                    print(f"Couldn't find variable {e}! Did you define it in vars_list in ExperimentVariables?")
                else:
                    print(f"Exception {e}") # todo: replace with raise statement

class ExperimentVariables(EnvExperiment):

    def build(self):

        Variable = namedtuple("Variable", "name value value_type kwargs group")
        """
        An object for defining an experiment variable. This simplifies creating the datasets for each var.
        
        name: the variable name. must be a valid python name
        value: the value of the variable. this is only used the very first time the code is run after you add
            a new variable. Every subsequent run, the value will be pulled from the dataset that is created.
        value_type: the artiq value function, e.g. NumberValue or StringValue
        kwargs: a dictionary of keyword arguments used by the value function.
        group: set to None or a string. If a string, all Variables with this same string will be grouped 
            in the GUI.
            
        Example:
        Variable("f_FORT", 210.0 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'}, "FORT switching AOM")
        """

        # the list of variables which we will use to make datasets of the same name stored in an hdf file.
        # these datasets can be referenced by other experiments instead of declaring the variable locally in each
        # experiment.
        self.vars_list = [

            Variable("n_measurements", 50, NumberValue, {'type': 'int', 'ndecimals':0, 'step':1, 'scale':1}, "general"),
            Variable("require_atom_loading_to_advance", False, BooleanValue, {}, "general"),
            Variable("require_atom_loading_to_advance_in_single_photon_exp", True, BooleanValue, {}, "single photon experiment"),

            Variable("n_excitation_attempts", 10, NumberValue, {'type': 'int', 'ndecimals': 0, 'step': 1, 'scale': 1},
                     "single photon experiment"),
            Variable("max_excitation_cycles", 200, NumberValue, {'type': 'int', 'ndecimals': 0, 'step': 1, 'scale': 1},
                     "single photon experiment"),
            Variable("atom_check_every_n", 5, NumberValue, {'type': 'int', 'ndecimals': 0, 'step': 1, 'scale': 1},
                     "single photon experiment"),

            Variable("record_every_shot", True, BooleanValue, {}, "single photon experiment"),

            # debugging
            Variable("dummy_variable", 0.0, NumberValue, {'type': 'float'}, "debugging"),

            # FORT AOM
            Variable("f_FORT", 80.0 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'}, "FORT AOM"),
            Variable("p_FORT_loading", -4, NumberValue,
                     {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1}, "FORT AOM"),
            Variable("p_FORT_RO", 1, NumberValue, {'type': 'float', 'unit': "(fractional)", 'scale': 1, 'ndecimals': 1},
                     "FORT AOM"),
            Variable("p_FORT_holding", 1, NumberValue,
                     {'type': 'float', 'unit': "(fractional)", 'scale': 1, 'ndecimals': 2},
                     "FORT AOM"),
            Variable("p_FORT_PGC", 1, NumberValue, {'type': 'float', 'unit': "(fractional)", 'scale': 1, 'ndecimals': 2},
                     "FORT AOM"),
            Variable("p_FORT_blowaway", 0.5, NumberValue,
                     {'type': 'float', 'unit': "(fractional)", 'scale': 1, 'ndecimals': 2},
                     "FORT AOM"),
            Variable("p_FORT_OP", 0.5, NumberValue,
                     {'type': 'float', 'unit': "(fractional)", 'scale': 1, 'ndecimals': 2},
                     "FORT AOM"),

            # Cooling double pass AOM
            Variable("f_cooling_DP_MOT", 111.0 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'},
                     "Cooling double pass AOM"),
            Variable("f_cooling_DP_MOT_phase2", 111.0 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'},
                     "Cooling double pass AOM"),
            Variable("f_cooling_DP_RO", 111.0 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'},
                     "Cooling double pass AOM"),
            Variable("f_cooling_DP_PGC", 111.0 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'},
                     "Cooling double pass AOM"),
            Variable("f_cooling_DP_resonant_2_to_3", 104.25 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'},
                     "Cooling double pass AOM"),
            Variable("f_cooling_DP_blowaway", 112.0 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'},
                     "Cooling double pass AOM"),
            Variable("p_cooling_DP_MOT", -4, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "Cooling double pass AOM"),
            Variable("p_cooling_DP_RO", 0.9, NumberValue, {'type': 'float', 'unit': "(fractional 0.0 to 1.0)",
                                                           'scale': 1, 'ndecimals': 1},
                     "Cooling double pass AOM"),
            Variable("p_cooling_DP_PGC", 0.9, NumberValue, {'type': 'float', 'unit': "(fractional 0.0 to 1.0)",
                                                            'scale': 1, 'ndecimals': 1},
                     "Cooling double pass AOM"),

            # D1 optical pumping, pumping repumper, and excitation
            Variable("f_D1_pumping_DP", 368 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz', 'ndecimals': 3},
                     "OP and excitation AOMs"),
            Variable("p_D1_pumping_DP", -9.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "OP and excitation AOMs"),
            Variable("f_pumping_repump", 345.10 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'},
                     "OP and excitation AOMs"),
            Variable("p_pumping_repump", -8.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "OP and excitation AOMs"),
            Variable("f_excitation", 150.50 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'},
                     "OP and excitation AOMs"),
            Variable("p_excitation", -8.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "OP and excitation AOMs"),

            # Fiber AOMs
            Variable("AOM_A1_freq", 78.5 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz', 'ndecimals': 3},
                     "Fiber AOMs"),
            Variable("p_AOM_A1", -13.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "Fiber AOMs"),
            Variable("AOM_A2_freq", 78.0 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz', 'ndecimals': 3},
                     "Fiber AOMs"),
            Variable("p_AOM_A2", -13.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "Fiber AOMs"),
            Variable("AOM_A3_freq", 78.503 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz', 'ndecimals': 3},
                     "Fiber AOMs"),
            Variable("p_AOM_A3", -13.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "Fiber AOMs"),
            Variable("AOM_A4_freq", 78.504 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz', 'ndecimals': 3},
                     "Fiber AOMs"),
            Variable("p_AOM_A4", -13.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "Fiber AOMs"),
            Variable("AOM_A5_freq", 78.505 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz', 'ndecimals': 3},
                     "Fiber AOMs"),
            Variable("p_AOM_A5", -13.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "Fiber AOMs"),
            Variable("AOM_A6_freq", 78.506 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz', 'ndecimals': 3},
                     "Fiber AOMs"),
            Variable("p_AOM_A6", -13.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "Fiber AOMs"),

            # # Pumping Repump
            # Variable("p_pumping_repump_A5", -20.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
            #          "Fiber AOMs"),
            # Variable("p_pumping_repump_A6", -8.0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
            #          "Fiber AOMs"),


            # Microwaves
            # assumes the microwave source is set at 6.5 GHz
            Variable("f_microwaves_dds", 334.682 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz', 'ndecimals': 7},
                     "Microwaves"),
            Variable("f_microwaves_detuning", 2.261 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz', 'ndecimals': 7},
                     "Microwaves"),
            Variable("p_microwaves", 0, NumberValue, {'type': 'float', 'unit': "dBm", 'scale': 1, 'ndecimals': 1},
                     "Microwaves"),

            # Coils - MOT
            Variable("AZ_bottom_volts_MOT", 4.896, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),
            Variable("AZ_top_volts_MOT", 5.405, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),
            Variable("AX_volts_MOT", 0.419, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),
            Variable("AY_volts_MOT", 0.255, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),

            Variable("AZ_bottom_volts_MOT_phase2", 4.896, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "MOT coil settings"),
            Variable("AZ_top_volts_MOT_phase2", 5.405, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "MOT coil settings"),
            Variable("AX_volts_MOT_phase2", 0.419, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "MOT coil settings"),
            Variable("AY_volts_MOT_phase2", 0.255, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "MOT coil settings"),

            # Coils - PGC
            Variable("AZ_bottom_volts_PGC", 0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),
            Variable("AZ_top_volts_PGC", 0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),
            Variable("AX_volts_PGC", -0.11, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),
            Variable("AY_volts_PGC", 0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),

            # Coils - Readout
            Variable("AZ_bottom_volts_RO", 0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),
            Variable("AZ_top_volts_RO", 0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),
            Variable("AX_volts_RO", 0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),
            Variable("AY_volts_RO", 0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals':3},
                     "MOT coil settings"),

            # Coils - state prep and science
            Variable("AZ_bottom_volts_OP", 0.0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            Variable("AZ_top_volts_OP", 0.0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            Variable("AX_volts_OP", 0.0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            # about 1 G along the parabolic mirror axis.
            Variable("AY_volts_OP", 2.105, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            Variable("AZ_bottom_volts_blowaway", 0.0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            Variable("AZ_top_volts_blowaway", 0.0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            Variable("AX_volts_blowaway", 0.0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            # about 1 G along the parabolic mirror axis.
            Variable("AY_volts_blowaway", 2.105, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            Variable("AZ_bottom_volts_microwave", 0.0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            Variable("AZ_top_volts_microwave", 0.0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            Variable("AX_volts_microwave", 0.0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),
            # about 1 G along the parabolic mirror axis.
            Variable("AY_volts_microwave", 0, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Science coil settings"),

            # Timing
            Variable("t_MOT_loading", 500*ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_MOT_phase2", 500 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_FORT_loading", 50*ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_FORT_drop", 10 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("t_FORT_modulation", 10 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_SPCM_exposure", 50 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_SPCM_first_shot", 10 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_SPCM_second_shot", 10 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_SPCM_recool_and_shot", 10 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_delay_between_shots", 20 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_PGC_in_MOT", 50 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_MOT_dissipation", 3 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_blowaway", 50 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("t_pumping", 50 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("t_depumping", 50 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("t_microwave_pulse", 100 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("t_pi_microwave_pulse", 1.494 * us, NumberValue, {'type': 'float', 'unit': 'us', 'ndecimals': 3}, "Timing"),
            Variable("t_photon_collection_time", 1 * ns, NumberValue, {'type': 'float', 'unit': 'ns'}, "Timing"),
            Variable("t_photon_collection_delay", 1 * ns, NumberValue, {'type': 'float', 'unit': 'ns'}, "Timing"),
            Variable("t_exp_trigger", 1*ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_UV_pulse", 10 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),
            Variable("t_OP_chop_period", 1 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("t_OP_chop_offset", 0.5 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("t_BA_chop_period", 1 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("duty_cycle_OP", 0.35, NumberValue, {'type': 'float', 'unit': "(fractional 0.0 to 1.0)",
                                                             'scale': 1, 'ndecimals': 2}, "Timing"),
            Variable("duty_cycle_FORT", 0.65, NumberValue, {'type': 'float', 'unit': "(fractional 0.0 to 1.0)",
                                                          'scale': 1, 'ndecimals': 2}, "Timing"),
            Variable("t_RO_chop_period", 1.5 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("t_RO_chop_offset", -0.2 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("t_RO_gate_offset", -0.2 * us, NumberValue, {'type': 'float', 'unit': 'us'}, "Timing"),
            Variable("duty_cycle_RO", 0.2, NumberValue, {'type': 'float', 'unit': "(fractional 0.0 to 1.0)",
                                                          'scale': 1, 'ndecimals': 2}, "Timing"),
            Variable("t_excitation_pulse", 100 * ns, NumberValue, {'type': 'float', 'unit': 'ns'}, "Timing"),
            Variable("gate_start_offset_mu", 0, NumberValue, {'type': 'int', 'ndecimals': 0, 'scale': 1, 'step':1},
                     "Timing"),
            Variable("gate_switch_offset", 0*ns, NumberValue, {'type': 'int', 'unit':'ns',
                                                            'ndecimals': 0, 'scale': 1, 'step': 1},
                     "Timing"),
            Variable("t_recooling", 1 * ms, NumberValue, {'type': 'float', 'unit': 'ms'}, "Timing"),

            # Booleans
            Variable("no_first_shot", False, BooleanValue, {}, "Booleans"),
            Variable("do_PGC_in_MOT", False, BooleanValue, {}, "Booleans"),
            Variable("pumping_light_off", False, BooleanValue, {}, "Booleans"),
            Variable("blowaway_light_off", False, BooleanValue, {}, "Booleans"),
            Variable("FORT_on_at_MOT_start", False, BooleanValue, {}, "Booleans"),
            Variable("D1_off_in_OP_phase", False, BooleanValue, {}, "Booleans"),
            Variable("D1_off_in_depump_phase", False, BooleanValue, {}, "Booleans"),
            Variable("require_D1_lock_to_advance", False, BooleanValue, {}, "Booleans"),
            Variable("use_chopped_readout", True, BooleanValue, {}, "Booleans"),
            Variable("chopped_RO_light_off", False, BooleanValue, {}, "Booleans"),
            Variable("verify_OP_in_photon_experiment", False, BooleanValue, {}, "Booleans"),

            # Thresholds and cut-offs
            Variable("single_atom_threshold", 10000.0, NumberValue, {'type': 'float'}, "Thresholds and cut-offs"),

            # Set points
            Variable("set_point_PD1_AOM_A1", 0.427, NumberValue, {'type':'float','ndecimals':4}, "Set points"),
            Variable("set_point_PD2_AOM_A2", 0.148, NumberValue, {'type': 'float','ndecimals':4}, "Set points"),
            Variable("set_point_PD3_AOM_A3", 0.293, NumberValue, {'type': 'float','ndecimals':4}, "Set points"),
            Variable("set_point_PD4_AOM_A4", 0.286, NumberValue, {'type': 'float','ndecimals':4}, "Set points"),
            Variable("set_point_PD5_AOM_A5", 0.214, NumberValue, {'type': 'float','ndecimals':4}, "Set points"),
            Variable("set_point_PD6_AOM_A6", 0.296, NumberValue, {'type': 'float','ndecimals':4}, "Set points"),
            Variable('set_point_FORT_APD_loading', 0.051, NumberValue, {'type': 'float', 'ndecimals': 3}, "Set points"),
            Variable('set_point_FORT_APD_science', 0.04, NumberValue, {'type': 'float', 'ndecimals': 3}, "Set points"),
            Variable('set_point_FORT_MM_loading', 0.272, NumberValue, {'type': 'float','ndecimals':3}, "Set points"),
            Variable('set_point_FORT_MM_science', 0.2, NumberValue, {'type': 'float', 'ndecimals': 3}, "Set points"),
            Variable('set_point_D1_SP', 0.175, NumberValue, {'type': 'float', 'ndecimals': 3}, "Set points"),
            Variable('set_point_excitation', 1.0, NumberValue, {'type': 'float', 'ndecimals': 3}, "Set points"),

            # Plotting
            Variable("MOT_beam_monitor_points", 100, NumberValue, {'type': 'int', 'ndecimals': 0, 'scale': 1, 'step':1},
                     "Plotting"),
            Variable("ignore_first_n_histogram_points", 10, NumberValue, {'type': 'int', 'ndecimals': 0, 'scale': 1,
                     'step': 1}, "Plotting"),
            Variable("FORT_MM_monitor_sci_ref_volts", 0.25, NumberValue, {'type': 'float', 'ndecimals': 3, 'scale': 1,
                                                                          'step': 0.05}, "Plotting"),
            Variable("second_shot_hist_color", 'r', StringValue, {}, "Plotting"),

            # Laser feedback
            Variable("Luca_trigger_for_feedback_verification", False, BooleanValue, {}, "Laser feedback"),
            Variable("enable_laser_feedback", False, BooleanValue, {}, "Laser feedback"),
            Variable("slow_feedback_dds_list", "['dds_AOM_A1', 'dds_AOM_A2', 'dds_AOM_A3', "
                                          "'dds_AOM_A4','dds_AOM_A4','dds_AOM_A5',"
                                      "'dds_AOM_A6','dds_cooling_DP']", StringValue, {}, "Laser feedback"),
            Variable("fast_feedback_dds_list", "['dds_FORT']", StringValue, {}, "Laser feedback"),
            Variable("aom_feedback_averages", 4, NumberValue, {'type': 'int', 'ndecimals': 0, 'scale': 1, 'step':1},
                     "Laser feedback"),
            Variable("aom_feedback_iterations", 4, NumberValue, {'type': 'int', 'ndecimals': 0, 'scale': 1, 'step':1},
                     "Laser feedback"),
            Variable("Rigol_carrier_frequency", 40 * kHz, NumberValue, {'type': 'float', 'unit': 'kHz', 'ndecimals': 3},
                     "Rigol DG1022Z settings"),
            Variable("Rigol_FM_deviation", 10 * kHz, NumberValue, {'type': 'float', 'unit': 'kHz', 'ndecimals': 3},
                     "Rigol DG1022Z settings"),
            Variable("Rigol_V_DC", 0.595, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Rigol DG1022Z settings"),
            Variable("Rigol_V_pp", 0.03, NumberValue, {'type': 'float', 'unit': 'V', 'ndecimals': 3},
                     "Rigol DG1022Z settings"),
            Variable("f_Rigol_modulation", 10 * kHz, NumberValue, {'type': 'float', 'unit': 'kHz', 'ndecimals': 3},
                     "Rigol DG1022Z settings"),
            # Thorlabs Devices
            #todo: delete these in ExperimentVariables. These are defined in device_db.py
            Variable("K10CR1_FORT_HWP_SN", 55000759, NumberValue, {'type': 'int', 'ndecimals': 0, 'scale': 1, 'step':1},
                     "Thorlabs Devices"),
            Variable("K10CR1_FORT_QWP_SN", 55000740, NumberValue,
                     {'type': 'int', 'ndecimals': 0, 'scale': 1, 'step': 1},
                     "Thorlabs Devices"),
            Variable("K10CR1_780_HWP_SN", 55422044, NumberValue,
                     {'type': 'int', 'ndecimals': 0, 'scale': 1, 'step': 1},
                     "Thorlabs Devices"),
            Variable("K10CR1_780_QWP_SN", 55420984, NumberValue,
                     {'type': 'int', 'ndecimals': 0, 'scale': 1, 'step': 1},
                     "Thorlabs Devices"),
            Variable("hwp_move_by_deg", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 2, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("hwp_move_to_deg", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 2, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("qwp_move_by_deg", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 2, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("qwp_move_to_deg", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 2, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("target_qwp_deg", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 2, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("target_hwp_deg", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 2, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("deg_to_pos", 136533, NumberValue,
                     {'type': 'int', 'ndecimals': 0, 'scale': 1, 'step': 1},
                     "Thorlabs Devices"),
            Variable("best_HWP_to_H", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 2, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("best_QWP_to_H", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 2, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("best_852HWP_to_max", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 3, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("best_852QWP_to_max", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 3, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices"),
            Variable("best_852_power", 0, NumberValue,
                     {'type': 'float', 'ndecimals': 3, 'scale': 1, 'step': 0.5},
                     "Thorlabs Devices")
        ]

        # can only call get_dataset in build, but can only call set_dataset in run. so
        # we need to see which datasets exist here, create new ones if there are new vars,
        # and then in run we will set all the datasets with the values from the GUI.
        for var in self.vars_list:
            try:  # get the value of the variable from the dataset, if it exists
                value = self.get_dataset(var.name)
                self.setattr_argument(var.name, var.value_type(value, **var.kwargs), var.group)
            except Exception as e:
                if type(e) == KeyError:  # if the variable does not exist, create the dataset
                    print(f"Found new variable {e}! Adding argument.")
                    self.setattr_argument(var.name, var.value_type(var.value, **var.kwargs), var.group)
                else:
                    print(f"Exception {e}")

        self.setattr_argument('which_node', EnumerationValue(['alice','bob','two_nodes']), "general")

    def run(self):

        for var in self.vars_list:
            self.set_dataset(var.name, getattr(self, var.name), broadcast=True, persist=True)

        self.set_dataset("which_node", self.which_node, broadcast=True, persist=True)

        # dependent quantities

        # note that these detunings are with respect to the bare atomic states, i.e. no FORT AC Stark shift
        detuning_MOT_units_Gamma = (-2*self.get_dataset("f_cooling_DP_MOT")+130e6+78.5e6)/6.065e6
        self.set_dataset("detuning_MOT_units_Gamma", detuning_MOT_units_Gamma, broadcast=True, persist=True)

        detuning_RO_units_Gamma = (-2 * self.get_dataset("f_cooling_DP_RO") + 130e6 + 78.5e6) / 6.065e6
        self.set_dataset("detuning_RO_units_Gamma", detuning_RO_units_Gamma, broadcast=True, persist=True)

        detuning_PGC_units_Gamma = (-2 * self.get_dataset("f_cooling_DP_PGC") + 130e6 + 78.5e6) / 6.065e6
        self.set_dataset("detuning_PGC_units_Gamma", detuning_PGC_units_Gamma, broadcast=True, persist=True)

        microwave_frequency_GHz = self.get_dataset("f_microwaves_dds")*1e-9 + 6.5
        self.set_dataset("microwave_frequency_GHz", microwave_frequency_GHz, broadcast=True, persist=True)

        # assumes 6.5 GHz
        self.set_dataset("f_dds_clock_resonance", 334.682*MHz, broadcast=True, persist=True)
