module MyConfig

    export prod_basedir, tag, particle, energies, nevents_per_task, ntasks_per_energy, ntasks_per_job,
           base_mac, job_template, mjob_template, fqtunerdir, wcsimdir, g4dir, condadir, rootdir,
           get_energy_and_index

    prod_basedir      = "$(ENV["LUSTRE"])/CProfiles/"
    tag               = "ccin2p3"
    particle          = "mu-"
    energies          = range(1000, 3000, step=50)
    nevents_per_task  = 100
    ntasks_per_energy = 20
    ntasks_per_job    = 5
    base_mac          = "templates/cprofile_base.mac"
    job_template      = "templates/job_template.sh"
    mjob_template     = "templates/mjob_template.sh"

    fqtunerdir = "$(ENV["HOME"])/Software/WCSimFQTuner/build"
    wcsimdir   = "$(ENV["HOME"])/Software/WCSim/install"
    g4dir      = "$(ENV["HOME"])/Software/Geant4/install"
    condadir   = "$(ENV["HOME"])/Software/miniconda3"
    rootdir    = "$(ENV["HOME"])/Software/Root/install"

    get_energy_and_index(fname::String) = map( x->parse(Int, x.match)
                                            , getindex(collect(eachmatch(r"[0-9]+", basename(fname))), [1, 2]))

end
