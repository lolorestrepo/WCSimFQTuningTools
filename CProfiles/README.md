Tools for the production of **fiTQun** tuning files

## CProfiles (Cherenkov Profiles)

### Simulation

Contains scripts to generate de simulations
The simulations use  **WCSimFQTuner** with **/fqTune/mode cherenkovProfile** to save the angle of the cherenkov photons and their emission point from the start of the track for a fixed initial energy. Two angles are saved, true angle and weighted angle. (The weighted angle is used for te electron profile?)

local: run the tuner interactively

### TuningFile
1) Produce a single file containing the (angle, s) distributions for each energy. This is done through the genhist.cc root macro.
2) 