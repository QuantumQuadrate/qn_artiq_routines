"""
Define dma waveform functions which can be imported in experiments

Notes:
* think about whether it makes sense to import these into experiments
directly or to add them to the base experiment stuff. an easier option
is to just write the functions assuming that the variables used are
defined by a parent experiment we pass by reference. but maybe I should
stop and figure out how to import experiment variables into the namespace
when the base experiment set up is run.
"""

# in some other file
def readout():
    dds_FORT.sw.on()

    print(...)




# in my exp file

import from function readout

class Example(EnvExperiment):

    def build(self):

        self.base.build()

        dds_FORT.sw.on()
    readout(experiment=self)
