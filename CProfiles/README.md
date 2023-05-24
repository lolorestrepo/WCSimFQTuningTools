Tools for the production of **fiTQun** tuning files

## CProfiles (Cherenkov Profiles)

### Simulation

Contains scripts to generate de simulations.
The simulations use  **WCSimFQTuner** with **/fqTune/mode cherenkovProfile** to save the angle of the cherenkov photons and their emission point from the start of the track for a fixed initial energy. Two angles are saved, true angle and weighted angle. (The weighted angle is used for te electron profile?)

To run the simulations locally see **Simulation/local**.

The **Simulation** directory contains scripts to run the simulations in a HPC cluster. The procedure is as follows:
1) Modify the **config.jl** file parameters as required
2) Modify the macro, job and mjob templates as required
3) Run the scripts in the following order to produce the jobs

    1) `julia generate_macro_files.jl`
    2) `julia generate_job_files.jl`
    3) `julia generate_mjob_files.jl`
4) Finally modify and run the job launcher script `julia launch_jobs.jl`


### TuningFile
Produce a single file containing the (angle, s) distributions for each energy.
1) Setup Root 
2) Run the **genhist.cc** **ROOT** macro (hardcoding needed for file names and looping): it merges the cherenkov profiles for all the simulated momenta. To run it open a **ROOT** shell and run 
    ```
    -L genhist.cc
    genhist(PID)
    ```
    where **PID** is the particle id number from PDG. 

3) Compile and run **integcprofile.cc** (g++ -o Macro Macro.C `root-config --cflags --libs`): integrates the cherenkov profiles for each momenta to produce the $I_n = \int g(s, cos\theta) s^n \quad (n=0, 1, 2)$. Hardcoding needed for momenta looping, file name if differs from the **genhist.cc** output and the paramaters
    - **R0max**: maximum $r_0$ distance, ie distance between PMT and track vertex.
    - **nR0bin**: number of $r_0$ bins.
    - **nth0bin**: number of $\theta_0$ bins, where $\theta_0 \in [,]$.
    
    Run it through the compiled executable with `integcprofile PID 0`. The second argument **must** be 0, otherwise the default hardcoded HK parameter values will be used instead.

4) Compile and run **fitcprofile.cc**
5) Run the **writecprofile.cc** **ROOT** macro