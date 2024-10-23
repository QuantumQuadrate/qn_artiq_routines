# qn-artiq-routines

ARTIQ routines for the quantum network experiment in the Saffman group at UW-Madison.

# I. Getting started
    if "Just give me a practical example already!": 
      Scroll to the BaseExperiment example subsection
    else:
      keep reading

When running experiments, it is good practice to checkout a branch other than main. Then, as you change and add new code, you are not messing with main, which should only contain tested code known to work. A good workflow is as follows:

1. Open git bash in C:\Networking Experiment\artiq codes\artiq-master\repository\qn_artiq_routines
2. If you are on branch main, do "checkout -b new_branchname". For example if I'm writing code to do an atom state tomography sequence, I might let new_branchname be "atom_tomography". This creates and switches to the new branch.
3. As you make changes day to day on the code, commit the changes by doing "git add my_modified_file.py" then 'git commit -m "my concise but clear message about the changes I made"'.
4. Occasionally push your changes to the remote hosted on GitHub by doing "git push" (if you created a new branch, you will be prompted to git push --set-upstream...).
5. When you are done with, e.g. your atom tomography code, go to the repo on GitHub and open a pull request. In most cases you should be able to automatically merge the pull request. Delete the old branch when given the option.
6. In git bash, checkout the main branch with "git checkout main" then update with "git pull".
7. Go back to step 2. Rinse and repeat.

# II. Experiment setup
In all but the most simple experiments, experiment code can get too expansive and too complicated if we initialize every hardware reference and experiment variable in one piece of code. Moreover, as there are often several experiments one will want to run, it is easy to make errors if there is not efficient code reuse and sharing of variables across experiments. The section below describes utility code aimed at streamlining both experiment management and design. 
## 1. ExperimentVariables 
For variables that you want to access across multiple experiments, they should be defined in ExperimentVariables.py. The ExperimentVariables experiment has a GUI for setting each of the variables you define, so you can update them there (and run the experiment with Submit) to update the variables for all of the experiments that use them. 
### Adding a new experiment variable
To add a new variable, open ExperimentVariables.py: 
1. In the experiment variables class, add an element to the list self.vars_list, following the form
  Variable(name, value, value_type, kwargs, group)
Variable is just a container for grouping the things we need to know about your variable. The arguments are
  * name: a valid python variable name
  * value: the value of the variable. 
  * value_type: the artiq value function, e.g. NumberValue or StringValue
  * kwargs: a dictionary of keyword arguments used by the value function.
  * group: set to None or a string. If a string, all Variables with this same string will be grouped in the GUI.
As an example, we could add a new variable for the frequency we want to set to a DDS channel for driving the cooling laser double-pass AOM:

Variable("f_cooling_DP_MOT", 111.0 * MHz, NumberValue, {'type': 'float', 'unit': 'MHz'},
                     "Cooling double pass AOM"),
                     
Save the file. Now, in the ARTIQ dashboard, open ExperimentVariables then click "Recompute arguments". The variable you added should now show up in the GUI. The GUI will save the value of the variable to a dataset. To change the value of the variable, change the value in the GUI and hit Submit. The value of 111.0 * MHz which you typed in the code is only used the very first time the code is run. After that, the value used will be whatever the last value was when ExperimentVariables was run (even if you close ARTIQ master and restart it).

### Using an experiment variable in your experiments
Use the BaseExperiment class to add all of the datasets created in ExperimentVariables as attributes of the same name to your experiment instance (i.e. if you define f_cooling_DP_MOT in ExperimentVariables, then you will be able to reference that as self.f_cooling_DP_MOT in your experiment). See the BaseExperiment section for details. 

There is also a low-level option if you insist on only importing a few of the variables. Do this in your build method:

    self.experiment.variables = [ # this goes in your build method at the beginning
              "f_cooling_DP_MOT", "p_cooling_DP_MOT",
              "AZ_bottom_volts_MOT", "AZ_top_volts_MOT", "AX_volts_MOT", "AY_volts_MOT",
              "t_MOT_loading"
          ]

    setattr_variables(self.experiment) # import this from ExperimentVariables

### Viewing the variables used in an experiment
Anytime an ARTIQ experiment is run using the variables that you defined in ExperimentVariables, the names and values of those variables will be stored in an h5 file generated when an experiment finishes. This file will be in results/todays_date/the_hour/ in your artiq directory. The file can be opened with HDFView. 
### Run an experiment with previous variable values
You can do this by clicking Load HDF in your experiment's GUI and navigating to the hdf file corresponding to the experiment whose values you want.
## 2. DeviceAliases
The names of hardware channels, e.g. the Urukul DDS channels, are defined in the device database device_db.py. The names of the hardware channels by default are, e.g. urukul0_ch0, but if ~ 10 different urukul channels are in use it is useful to reference them in experiments by self-documenting names such as dds_cooling_AOM. On the other hand, it is occasionally desired to run low-level tests that unambiguously reference a particular channel rather than its function; urukul0_ch0 is more useful now. DeviceAliases give the best of both worlds.
### Defining device aliases
A dictionary mapping aliases of the urukul channels and their defaults lives in utilities/config/your_node/. Entries can be modified, added, or removed, and the new alias mapping will be used everywhere the aliases are used in your experiments. The alias map looks like this:

    ALIAS_MAP = {
          "dds_FORT": "urukul0_ch0", # the key is the alias, and the value is the name given in device_db.py
          "dds_cooling_DP": "urukul0_ch1",
          "dds_D1_pumping_DP": "urukul0_ch2",
          ...
          "dds_AOM_A4": "urukul2_ch0",
          "dds_AOM_A5": "urukul2_ch1"
    }

and each channel that you define here must have a set of defined defaults which reference datasets (likely defined in ExperimentVariables), for example:

    "DDS_DEFAULTS": {
        "dds_FORT": {
          "frequency": "f_FORT",
          "power": "p_FORT_loading"
        },
        "dds_cooling_DP": {
          "frequency": "f_cooling_DP_MOT",
          "power": "p_cooling_DP_MOT"
        }

where p and f are shorthands for frequency (in Hz) and power (in dBm).
    
Note: you can technically do this with any of the hardware, but for now I have only used it for urukul channels, since other devices like the Sampler and Zotino are initialized by card but have many different functions for the channels on a given card. TTL channels are like urukul channels in this respect, and I plan to use aliases for them in the future: ttl_camera_trigger, ttl_debug, ttl_laser_unlocked, etc.

### Using the device aliases in experiments
Again, using the BaseExperiment class will add all of the hardware aliases as attributes to your experiment, as described below. 


If you insist on importing the aliases "manually", use the block below to set the device attributes (i.e. this block is used instead of statements like self.setattr_device("urukul0_ch1")):

        self.experiment.named_devices = DeviceAliases( # this goes in your build function. you need to import DeviceAliases from DeviceAliases.
            experiment=self.experiment,
            device_aliases=[ # the channels you want to use, which are already defined in ALIAS_MAP
                'dds_FORT',
                'dds_D1_pumping_DP',
                'dds_cooling_DP',
                'dds_pumping_repump',
                *[f'dds_AOM_A{i + 1}' for i in range(6)]  # the fiber AOMs
            ]
        ) # after this block carry on as usual. e.g. in run, do things like self.dds_FORT.init(), self.dds_FORT.set(210*MHz, amplitude=...)

Note that even after you have run the above block, you can still directly reference the channel by the name defined in the database device_db.py. This is because DeviceAliases is just creating attributes of your experiment that reference the device attribute with the database name. 

## III. BaseExperiment

The easiest way to set up an experiment which has access to all* (see build section below) of the experiment variables and hardare is to use the BaseExperiment class. By creating an instance of BaseExperiment inside your ARTIQ experiment, the experiment variables and devices will be added as attributes to your experiment. 

Note in particular the method calls for base.build, base.prepare, base.set_datasets_from_gui_args, base.initialize_hardware, and write_results.

  1. Example with BaseExperiment
  
    from artiq.experiment import *
    from BaseExperiment import BaseExperiment

    class MyExp(EnvExperiment):

      def build(self):
          # gets variables and devices
          self.base = BaseExperiment(experiment=self)
          self.base.build()

          # add some GUI arguments particular to this experiment
          self.setattr_argument("n_measurements", NumberValue(100, ndecimals=0, step=1))
          self.setattr_argument('scan_variable1_name', StringValue('t_blowaway'))
          self.setattr_argument("control_experiment", BooleanValue(False), "Control experiment")

          # finally, call this base method to add all of the GUI arguments as datasets so that they are saved to your h5 file at the end
          self.base.set_datasets_from_gui_args()

      def prepare(self):
          self.base.prepare()

      @kernel
      def run(self):
          self.base.initialize_hardware()
          # now we can do physics
          self.dds_cooling_DP.sw.on()
          
          ...
          self.dds_FORT.sw.on()
          delay(50*ms)
          self.dds_FORT.sw.off()

          # finally, call write_results()
          self.write_results()

  2. build:
  The variables and devices that your BaseExperiment includes are determined by what is defined in BaseExperiment.py. In the BaseExperiment class, in build, 
    a. the variables that get imported are added by name in the list self.experiment.variables. These names must exist as datasets created by ExperimentVariables.
    b. devices that you do not want to use with aliases are given by name in the list devices_no_alias.
    c. devices that you want to use by alias are given by name in the DeviceAliases call.
  3. set_datasets_from_gui_args:
  This does what it says. This sets your GUI arguments as datasets so they will be stored in the h5 save file, which would not otherwise be the case.
  4. prepare: 
  The prepare method takes care of any calculations or math that you want to do for every experiment. This includes things such as converting SI units to machine units (mu) and calculating amplitudes in volts for the urukul set method given a power in dBm.
  5. initialize_hardware: takes care of the remaining setup of hardware, i.e. things that happen on the kernel. Urukul default settings are defined in a dictionary in DeviceAliases, and are grabbed the base experiment when the initialize_hardware method is run. Any other devices that need initialization such as Samplers or Zotinos, or setting TTL channels as input or output, must happen here.
  6. write_results: this is a method of the experiment itself, but is added in BaseExperiment for reasons. This is a redundant filesave, because sometimes ARTIQ corrupts the h5 file in the process of saving, and you lose all of your data. This is not anything to do with the data itself, but a timeout when ARTIQ is doing cleanup. To avoid running into this, you just call this, which will save the file before ARTIQ does cleanup.
  
## IV. Node specific behavior

Node specific settings should be configured with the json files in utilities/config, e.g. for feedback settings and dds channel defaults. There are also several sections in BaseExperiment which are broken up into if-elif-else blocks to accomplish node-dependent setup where necessary. You should always strive to write your experiments so that node-dependent behavior is only apparent in the BaseExperiment class, rather than in the experiment itself. This helps ensure that behavior is the same between nodes and that we maximize the amount of code in common.

## V. Subroutines
The scripts in the main directory cannot be run as standalone experiments, but are
intended to be used within a parent artiq experiment, where the experiment instance
is passed by reference to each subroutine. These scripts are to serve the purpose of 
defining frequently used parts of an experiment sequence (e.g. optical pumping,
dipole trap loading, polarization gradient cooling) so they are easily deployable 
for reuse in a number of experiments.

Note that the subroutines in this repository will only function when called by a parent experiment which has declared all of the hardware of the same name which is used in the subroutine (i.e. if the subroutine uses a DDS channel cooling_dds, then the parent experiment should have initialized cooling_dds).

There is also a list of defaults for each channel which point to variables defined in ExperimentVariables.

## V. Old example experiments
The examples folder has some standalone example code (in the sense that the script
contains the whole experiment workflow, including the parent class which calls the
subroutine). These can be run provided the device database uses the same names as 
are used in the example, as mentioned above.


