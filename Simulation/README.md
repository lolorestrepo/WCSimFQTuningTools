### **Tunning simulations**

This directory contains scripts to generate and run the simulations.
The simulations are configured by the `config.jl` file, which defines the simulation parameters for each type of simulation: **CProfiles**, **STable**, ...

The simulations use **WCSimFQTuner**, which requires a macro of configuration parameters. This utility creates the configuration parameters, the tasks for each desired configuration and the jobs to be run.

The names of the configuration parameters are autoexplanatory. The logic is as follows: each simulation type requires its own set of configuration variables `config_variables`. For example the **CProfiles** simulations are performed as a function of the initial particle type and energy, and this utility will generate macro files and jobs for each particle-energy combination.

The procedure is as follows:
1) First modify the `config.jl` file parameters as required.
2) Modify the macro, task and job templates as desired. Notice that we define two types of macros, the *base* and the *config* macros. The *base* macro defines the set of fixed parameters for the particular simulation type, and the *config* the *seteable* parameters, which are defined with a preceding `$` symbol.
3) Run the scripts in the following order to produce the jobs
    1) `julia generate_macro_files.jl`
    2) `julia generate_task_files.jl`
    3) `julia generate_job_files.jl`
4) The jobs can be launched through the `launch_jobs.jl` script.

To run the simulations interactivelly in you local machine, see **Simulations/local**.