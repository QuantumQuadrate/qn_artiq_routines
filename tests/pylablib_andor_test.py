from time import sleep
import numpy as np
import matplotlib.pyplot as plt

import pylablib as pll
pll.par["devices/dlls/andor_sdk2"] = "C:\\Program Files\\Andor SOLIS"
from pylablib.devices import Andor

print("found",Andor.get_cameras_number_SDK2(),"camera(s)")

cam = Andor.AndorSDK2Camera(temperature=-20, fan_mode="full")

try:
    # whatever we want to do
    # print(dir(cam))

    # print out camera params
    print(cam.is_opened())
    print(cam.get_device_info())
    print(cam.get_status())
    print(cam.get_temperature_status())
    print(cam.get_temperature())
    print(cam.get_EMCCD_gain())
    print(cam.get_exposure())

    # warm up
    while not cam.get_temperature_status() == "stabilized":
      sleep(0.1)
    print("temperature stabilized")
    cam.set_trigger_mode("int") # software trigger not available for Luca
    cam.set_acquisition_mode("single")
    cam.set_exposure(0.1)
    cam.set_EMCCD_gain(100)
    cam.start_acquisition()
    acquired = 0
    (acquired, unread, skipped, size) = cam.get_frames_status()

    # this always says acquired is 0. weird
    for stat, lbl in zip((acquired, unread, skipped, size),
                         ("acquired", "unread", "skipped", "size")):
        print(lbl,": ",stat)
    cam.stop_acquisition()
    frame = np.array(cam.grab())
    print(frame.shape)
    plt.imshow(frame[0])
    plt.show()

    # show some ROIs
    fig, axes = plt.subplots(nrows=2,ncols=2)

    ROI1 = ((280,335), (500,800))
    ROI2 = ((540,595), (500,800))
    ROI3 = ((280,335), (500,800))
    ROI4 = ((280,330), (300,600))

    ROIs = [ROI1, ROI2, ROI3, ROI4]

    for ax,ROI in zip(axes.flat,ROIs):
        ROI_rows,ROI_cols = ROI
        ax.imshow(frame[0,slice(*ROI_rows),slice(*ROI_cols)])
        ax.set_aspect('auto')
    plt.show()

except KeyboardInterrupt:
    print("keyboard interrupt. quitting nicely.")

    """
    It is important to close all camera connections before finishing 
    your script. Otherwise, DLL resources might become permanently 
    blocked, and the only way to solve it would be to restart the PC.
    """
    cam.close()

"""
attributes of the AndorSDK2Camera instance:
['Error', 'FrameTransferError', 'NoParameterCaller', 'TimeoutError', '_TAcqTimings', 
'_TFrameInfo', '_TFramesStatus',  
'_acq_setup_requested', '_acqmode_caps', '_acqmode_desc', '_add_device_variable', 
'_add_info_variable', '_add_parameter_class', '_add_settings_variable', 
'_add_status_variable', '_adjustable_frameinfo_period', '_as_parameter_class', 
'_buffer_size', '_call_without_parameters', '_check_option', '_clear_pausing_acquisition',
'_close_on_error', '_convert_frame_format', '_convert_frame_info', '_convert_indexing',
'_cpar', '_default_acq_params', '_default_frame_format', '_default_frame_info', 
'_default_frameinfo_format', '_default_image_dtype', '_default_image_indexing', 
'_device_var_ignore_error', '_device_vars', '_device_vars_order', '_empty_frame_info',
'_ensure_acquisition_parameters', '_find_min_roi_end', '_frame_counter', 
'_frame_format', '_frame_info_to_namedtuple', '_frameinfo_fields', 
'_frameinfo_fields_mask', '_frameinfo_format', '_frameinfo_include_fields',
'_frameinfo_period', '_get_acquired_frames', '_get_capabilities_n', 
'_get_connection_parameters', '_get_data_dimensions_rc', '_get_device_variables',
'_get_grab_acquisition_parameters', '_get_updated_acquisition_parameters', 
'_has_option', '_image_indexing', '_initial_setup_ext_trigger',
'_initial_setup_temperature', '_int_to_enumlst', '_minh', '_minv', '_multiplex_func', 
'_opid', '_p_acq_mode', '_p_fan_mode', '_p_frame_format', '_p_frame_wait_mode', 
'_p_frameinfo_format', '_p_indexing', '_p_missing_frame', '_p_read_mode', 
'_p_shutter_mode', '_p_status', '_p_temp_status', '_p_trigger_mode', '_parameters', 
'_parse_capabilities', '_read_frames', '_readmode_caps', '_readmode_desc',
'_remove_device_variable', '_replace_parameter_class', '_restart_acq_after_pause',
'_select_camera', '_setup_default_settings', '_setup_parameter_classes', 
'_start_amp_mode', '_start_fan_mode', '_start_temperature', '_strict_option_check',
'_support_chunks', '_trigger_caps', '_trigger_desc', '_trim_images_range', 
'_truncate_amp_mode', '_truncate_roi_axis', '_truncate_roi_binning', 
'_update_device_variable_order', '_wait_for_next_frame', '_wait_sleep_period', 
'_wap', '_wip', '_wop', '_zero_frame', 'acquisition_in_progress', 'apply_settings',
'capabilities', 'clear_acquisition', 'close', 'device_info', 'dv',
'enable_frame_transfer_mode', 'get_EMCCD_gain', 'get_accum_mode_parameters', 
'get_acquisition_mode', 'get_acquisition_parameters', 'get_acquisition_progress',
'get_all_amp_modes', 'get_all_vsspeeds', 'get_amp_mode', 'get_buffer_size', 
'get_capabilities', 'get_channel', 'get_channel_bitdepth', 
'get_cont_mode_parameters', 'get_cycle_timings', 'get_data_dimensions',
'get_detector_size', 'get_device_info', 'get_device_variable', 'get_exposure', 
'get_ext_trigger_parameters', 'get_fan_mode', 'get_fast_kinetic_mode_parameters', 
'get_frame_format', 'get_frame_info_fields', 'get_frame_info_format', 
'get_frame_info_period', 'get_frame_period', 'get_frame_timings', 'get_frames_status', 
'get_full_info', 'get_full_status', 'get_hsspeed', 'get_hsspeed_frequency', 
'get_image_indexing', 'get_image_mode_parameters', 'get_keepclean_time', 
'get_kinetic_mode_parameters', 'get_max_vsspeed', 'get_min_shutter_times', 
'get_multi_track_mode_parameters', 'get_new_images_range', 'get_oamp', 'get_oamp_desc', 
'get_pixel_size', 'get_preamp', 'get_preamp_gain', 'get_random_track_mode_parameters', 
'get_read_mode', 'get_readout_time', 'get_roi', 'get_roi_limits', 'get_settings', 
'get_shutter', 'get_shutter_parameters', 'get_single_track_mode_parameters', 'get_status',
'get_temperature', 'get_temperature_range', 'get_temperature_setpoint', 
'get_temperature_status', 'get_trigger_level_limits', 'get_trigger_mode', 
'get_vsspeed', 'get_vsspeed_period', 'grab', 'handle', 'idx', 'ini_path', 
'init_amp_mode', 'is_acquisition_setup', 'is_cooler_on', 'is_frame_transfer_enabled',
'is_opened', 'open', 'pausing_acquisition', 'read_in_aux_port', 'read_multiple_images', 
'read_newest_image', 'read_oldest_image', 'send_software_trigger', 'set_EMCCD_gain', 
'set_acquisition_mode', 'set_amp_mode', 'set_cooler', 'set_device_variable', 
'set_exposure', 'set_fan_mode', 'set_frame_format', 'set_frame_info_format', 
'set_frame_info_period', 'set_frame_period', 'set_image_indexing', 'set_out_aux_port', 
'set_read_mode', 'set_roi', 'set_temperature', 'set_trigger_mode', 'set_vsspeed', 
'setup_accum_mode', 'setup_acquisition', 'setup_cont_mode', 'setup_ext_trigger', 
'setup_fast_kinetic_mode', 'setup_image_mode', 'setup_kinetic_mode', 
'setup_multi_track_mode', 'setup_random_track_mode', 'setup_shutter', 
'setup_single_track_mode', 'snap', 'start_acquisition', 'stop_acquisition', 
'wait_for_frame']
"""

