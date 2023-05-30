1) ROOT setup (in order to be able to `import ROOT`) 
2) Run the **genhist.cc** **ROOT** macro (hardcoding needed for file names and looping): it merges the cherenkov profiles for all the simulated momenta. To run it open a **ROOT** shell and run 
    ```
    .L genhist.cc
    genhist(PID)
    ```
    where **PID** is the particle id number from PDG.

3) Compile and run **integcprofile.cc** (`g++ -o integcprofile integcprofile.cc 'root-config --cflags --libs\'`): integrates the cherenkov profiles for each momenta to produce the $I_n = \int g(s, cos\theta) s^n ~~ (n=0, 1, 2)$. Hardcoding needed for: 1) momenta values in the loop; 2) input file name if it is different from the **genhist.cc** output; 3) the parameters
    - **R0max**: maximum $r_0$ distance, ie distance between PMT and track vertex.
    - **nR0bin**: number of $r_0$ bins.
    - **nth0bin**: number of $\cos \theta_0$ bins, where $\cos \theta_0 \in [-1, +1]$.
    
    Run it through the compiled executable with `integcprofile PID 0`. The second argument **must** be 0, otherwise the default hardcoded HK parameter values will be used instead.

4) Compile and run **fitcprofile.cc**. Script difficult to understand. Fits each Cherenkov profile, the number of partitions of each fit is hardcoded to 1. Run in through the executable `fitcprofile PID`
5) Run the **writecprofile.cc** **ROOT** macro. Also difficult to understand. Fits the rest of the variables: the isotropic integrals, the gNphot and gsthr.

### **Comments/Bugs/Unknon behaviours**

- `hprofinf`
    - the entry number 3 is hardcoded for each particle type
    - the entry number 4 is hardcoded to 0.2 if not passed as argument in **fitcprofile.cc**
    - very difficult to understand what entries 4, 5 are from **fitcprofile.cc**

- the upper bound of the last section in the fits for the 1D graphs `gNphot`, `gsthr`, `hI_iso_1` and `hI_iso_2` is wrongly taken from the 5 entry in `hprofinf`.

- only lower and upper section bounds are set in **writecprof.cc**, but middle bound 5.5 and 7.5 (line 24) are fixed. Not sure if it is the intended behaviour, but it might happen that the upper section bound is lower than 7.5 which might be probematic in fiTQun. Anyway the these values are hardcoded.