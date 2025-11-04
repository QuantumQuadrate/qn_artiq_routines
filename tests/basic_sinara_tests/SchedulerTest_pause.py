from artiq.experiment import *


class SoftTerminateExperiment(EnvExperiment):

    def build(self):
        self.setattr_device('core')
        self.setattr_device('scheduler')

    def run(self):
        # This example stores progress to continue after a pause
        num_iterations = 20

        try:
            while num_iterations > 0:
                while self.scheduler.check_pause():
                    # Pause the scan
                    print('pausing')
                    self.core.comm.close()  # Close communications before pausing
                    self.scheduler.pause()  # Can raise a TerminationRequested exception
                    print('resuming')

                # Some long kernel
                num_iterations = self._run(num_iterations)

        except TerminationRequested:
            print('gracefully terminating such that hdf5 file will be written')

    @kernel
    def _run(self, num_iterations: TInt32):
        while num_iterations > 0:
            # Some non-interruptable work
            self.msg(num_iterations)
            num_iterations -= 1
            self.core.reset()
            delay(1 * s)
            self.core.wait_until_mu(now_mu())

            # Check if there is a pause condition once in a while
            if self.scheduler.check_pause():
                # Interrupt current work
                break

        # Return progress counter
        return num_iterations

    @rpc(flags={'async'})
    def msg(self, i):
        print('doing stuff... ({:d})'.format(i))
