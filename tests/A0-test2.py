"""
Dummy experiment to be deleted
"""
for variable2_value in self.scan_sequence2:

    self.set_dataset("iteration", iteration, broadcast=True)

    if self.scan_variable2 != None:
        setattr(self, self.scan_variable2, variable2_value)
        logging.info(f"current iteration: {self.scan_variable2_name} ={variable2_value}")

    self.experiment_function()
    self.write_results({'name': "parent_rid_" + f"{self.parent_rid}" + "_" + self.experiment_name[:-11] + "_scan_over_" + self.scan_var_filesuffix})

    iteration += 1

    ## Run Health Check
    if self.run_health_check_and_schedule:
        if self.health_check_every_n_ite and (iteration % self.every_n_ite) == 0:
            print(f"running health check after iteration #{iteration - 1}, scanning over {self.scan_variable1}: {variable1_value}")

            failed_scans = self.health_check_microwave_freqs()

            # write and overwrite the health check results
            self.write_results({'name': "parent_rid_" + f"{self.parent_rid}" + "_" + self.experiment_name[
                                                                                     :-11] + "_scan_over_" + self.scan_var_filesuffix})

            if failed_scans:
                print("These scans need re-optimisation:", failed_scans)
                ###todo: if health check fail, i) schedule optimization
                for scan_name in failed_scans:
                    print("Scheduling experiment to optimize: ", scan_name)
                    override_args = copy.deepcopy(self.override_arguments_for_scheduling_optimization_dict)
                    override_args[scan_name] = True
                    self.submit_optimization_scans(override_arguments=override_args)
                    self.override_arguments_for_scheduling_optimization_dict[scan_name] = False

                ###todo: if health check fail, ii) schedule resuming experiment
                self.submit_resume_scan_after_optimization(current_iteration=iteration - 1)

                ### after scheduling scans above, it terminates the current experiment
                self.scheduler.request_termination(self.scheduler.rid)

            else:
                print("All microwave health checks passed. Resuming next scan")
                ### update everything back to previous setting - done inside healthcheck if passed

                ### Override variables again to avoid conflict
                # override specific variables. this will apply to the entire scan, so it is outside the loops
                for variable, value in self.override_ExperimentVariables_dict.items():
                    setattr(self, variable, value)

            ### terminating the experiment if health_check failed.
            self.core.comm.close()  # placing the hardware in a safe state and disconnecting it from the core device
            self.scheduler.pause()

print("****************    General Variable Scan DONE   *****************")