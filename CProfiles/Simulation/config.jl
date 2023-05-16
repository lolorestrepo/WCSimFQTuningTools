prod_basedir   = "$(ENV["HOME"])/WCTE/Production/Tuning/CProfiles/Simulation/prod/"
tag            = "ccin2p3"
particle       = "mu-"
energies       = range(1000, 2000, step=10)
nevents        = 1
ntasks_per_job = 2
base_mac       = "templates/cprofile_base.mac"
job_template   = "templates/job_template.sh"
mjob_template  = "templates/mjob_template.sh"

fqtunerdir = "$(ENV["HOME"])/WCTE/Software/WCSimFQTuner/build"
wcsimdir   = "$(ENV["HOME"])/WCTE/Software/WCSim/install"
g4dir      = "$(ENV["HOME"])/Software/Geant4/install"
rootdir    = "$(ENV["HOME"])/Software/ROOT/install"

get_fnumber(fname::String) = parse(Int, match(r"[0-9]+", basename(fname)).match)