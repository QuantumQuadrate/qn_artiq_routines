"""
adapted from Thorlabs' Camera_Examples\Python\grab_single_frame.py
"""

import numpy as np
import os
import cv2
from thorlabs_tsi_sdk.tl_camera import TLCameraSDK, OPERATION_MODE
import matplotlib.pyplot as plt
from datetime import datetime as dt
from PIL import Image


NUM_FRAMES = 1  # adjust to the desired number of frames
dll_parent_dir = 'C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\thorlabs'

os.add_dll_directory(dll_parent_dir + "\\dlls")
os.environ['PATH'] = dll_parent_dir + "\\dlls\\" + os.pathsep + os.environ['PATH']

t_ThorCam_exposure_ms = 100

t_experiment_run = dt.now().strftime("%Y%m%d_%H%M%S")
data_dir = 'C:\\Networking Experiment\\artiq codes\\artiq-master\\repository\\qn_artiq_routines\\thorlabs\\images'
file_suffix = 'MonitorMOTandExternalBeamPositions'
file_prefix = os.path.join(data_dir, t_experiment_run + '_' + file_suffix)

with TLCameraSDK() as sdk:
    available_cameras = sdk.discover_available_cameras()
    if len(available_cameras) < 1:
        print("no cameras detected")

    with sdk.open_camera(available_cameras[0]) as camera:
        camera.exposure_time_us = int(1e3*t_ThorCam_exposure_ms)
        camera.frames_per_trigger_zero_for_unlimited = 1  # start camera in continuous mode
        camera.operation_mode = 1  # hardware triggering
        camera.image_poll_timeout_ms = 1000  # polling timeout
        camera.gain = 300

        camera.arm(2)

        i = 0
        while True:
            try:
                frame = camera.get_pending_frame_or_null()
                if frame is not None:
                    print("frame #{} received!".format(frame.frame_count))
                    frame.image_buffer
                    image_buffer_copy = np.copy(frame.image_buffer)
                    numpy_shaped_image = image_buffer_copy.reshape(camera.image_height_pixels, camera.image_width_pixels)
                    nd_image_array = np.full((camera.image_height_pixels, camera.image_width_pixels, 3), 0, dtype=np.uint8)
                    nd_image_array[:, :, 0] = numpy_shaped_image
                    nd_image_array[:, :, 1] = numpy_shaped_image
                    nd_image_array[:, :, 2] = numpy_shaped_image

                    im = Image.fromarray(nd_image_array)
                    fname = file_prefix+str(i)+'.bmp'
                    im.save(fname)
                    print("saved fname")
                    i += 1

                else:
                    print(f"Unable to acquire image. Is USB traffic too high (e.g. live video from Luca)?")

            except KeyboardInterrupt:
                print("keyboard interrupt. quitting nicely.")
                cv2.waitKey(0)
                camera.disarm()
                break

        #  Because we are using the 'with' statement context-manager, disposal has been taken care of.