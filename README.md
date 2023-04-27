# qn-artiq-routines

ARTIQ routines for the quantum network experiment in the Saffman group at UW-Madison.

## Usage notes

The scripts in the main directory can not be run as standalone experiments, but are
intended to be used within a parent artiq experiment, where the experiment instance
is passed by reference to each subroutine. These scripts are to serve the purpose of 
defining frequently used parts of an experiment sequence (e.g. optical pumping,
dipole trap loading, polarization gradient cooling) so they are easily deployable 
for reuse in a number of experiments.

The names of hardware channels, e.g. Urukul DDS channels, are defined in the device
database device_db.py. The subroutines here will only function when called by a parent
experiment which has declared all of the hardware of the same name which is used
in the subroutine (i.e. if the subroutine uses a DDS channel cooling_dds, then the
parent experiment should have initialized cooling_dds).

The examples folder has some standalone example code (in the sense that the script
contains the whole experiment workflow, including the parent class which calls the
subroutine). These can be run provided the device database uses the same names as 
are used in the example, as mentioned above.


