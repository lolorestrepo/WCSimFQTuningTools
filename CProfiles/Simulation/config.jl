module MyConfig

       prod_basedir      = "$(ENV["LUSTRE"])/CProfiles/"
       tag               = "ccin2p3"
       particle          = "e-"
       energies          = range(100, 1000, step=10)
       nevents_per_task  = 200
       ntasks_per_energy = 50
       ntasks_per_job    = 10
       base_mac          = "templates/cprofile_base.mac"
       job_template      = "templates/job_template.sh"
       mjob_template     = "templates/mjob_template.sh"

       fqtunerdir = "$(ENV["HOME"])/Software/WCSimFQTuner/build"
       wcsimdir   = "$(ENV["HOME"])/Software/WCSim/install"
       g4dir      = "$(ENV["HOME"])/Software/Geant4/install"
       condadir   = "$(ENV["HOME"])/Software/miniconda3"
       rootdir    = "$(ENV["HOME"])/Software/Root/install"

       get_energy_and_index(fname::String) = 
       map(x->parse(Float64, x.match), getindex(collect(eachmatch(r"\d+\.*\d*+", basename(fname))), [1, 2]))

       export prod_basedir, tag, particle, energies, nevents_per_task, ntasks_per_energy, ntasks_per_job, base_mac, 
       job_template, mjob_template, fqtunerdir, wcsimdir, g4dir, condadir, rootdir, get_energy_and_index
end
