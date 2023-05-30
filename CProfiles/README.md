Tools for the production of **fiTQun** tuning files

## CProfiles (Cherenkov Profiles)

### **Simulation**

Contains scripts to generate simulations.
The simulations use  **WCSimFQTuner** with **/fqTune/mode cherenkovProfile** to save the angle of the cherenkov photons and their emission point from the start of the track for a fixed initial energy. Two angles are saved, true angle and weighted angle. (The weighted angle is used for the electron profile?)

To run the simulations locally see **Simulation/local**.

The **Simulation** directory contains scripts to run the simulations in a HPC cluster. The procedure is as follows:
1) Modify the **config.jl** file parameters as required
2) Modify the macro, job and mjob templates as required
3) Run the scripts in the following order to produce the jobs

    1) `julia generate_macro_files.jl`
    2) `julia generate_job_files.jl`
    3) `julia generate_mjob_files.jl`
4) Finally modify and run the job launcher script in the parent directory `julia launch_jobs.jl`


### **TuningFile**
Produce a single file containing the $(\cos \theta, s)$ distributions for each energy.
1) ROOT setup (in order to be able to `import ROOT`) 
2) Run `merge_cprofiles.py` script: `python merge_cprofiles.py /path/to/input/files/ tr [-v]`, where the input files must be named as **cprofile_{energy}MeV_{tag}_{number}.root**. It outputs **cprofiles_merged.root** which merges the cherenkov profiles for all the simulated momenta. The second argument stands for **true** or **weigthed** profiles.

3) Run `integrate_cprofiles.py`: integrates the cherenkov profiles for each momenta to produce the $I_n = \int g(s, \cos\theta) s^n ~~ (n=0, 1, 2)$. Reads the **cprofiles_merged.root** (default) created at 2. Returns a **cprofiles_integrals.root** file.

    Run it using `python integrate_cprofiles.py r0max, nr0bins, nth0bins [-i merged_file.root]` where:
    - **r0max**: maximum $r_0$ distance in cm, ie distance between PMT and track vertex.
    - **nr0bins**: number of $r_0$ bins.
    - **nth0bins**: number of $\cos \theta_0$ bins, with $\cos \theta_0 \in [-1, +1]$.

4) Run `fit_cprofile_integrals.py`. For each pair $(r_0, \cos \theta_0)$, it performs a polynomial fit of the integral $I_n$ vs the logarithm of the energy $\log E[\text{MeV}]$.

    Run it using `python fit_cprofile_integrals.py npars [-i integrals_file.root]` where `npars` is the number of parameters of the polynomial ($I_n = \sum^{n}_{j=1} j_iE^i$). Reads the **cprofiles_integrals.root** (default) created at 3. The output is **cprofiles_fits.root** containing the fitted parameters $j_i$ for each pair of $(r_0, \cos \theta_0)$.

    The output file format is the same as the required by fiTQun, but must be renamed.

#### **To-Do**
- [ ] Add errors to integral values and parameter fits
- [ ] Implement sectioned fit, ie partitioned energy range (really needed?)


### **Analysis**

Contains jupyter-notebooks to check distributions for each production step and compare with previous implementation.