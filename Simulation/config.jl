module Parameters

       export fqtunerdir, wcsimdir, g4dir, condadir, rootdir, prod_basedir
       export verbose, nevents_per_task, nsubtasks, ntasks_per_job
       export base_mac, config_mac, task_template, job_template
       export config_variables, queue_command, max_jobs_queue

       fqtunerdir   = "$(ENV["HOME"])/Software/WCSimFQTuner/build"
       wcsimdir     = "$(ENV["HOME"])/Software/WCSim/install"
       g4dir        = "$(ENV["HOME"])/Software/Geant4/install"
       condadir     = "$(ENV["HOME"])/Software/miniconda3"
       rootdir      = "$(ENV["HOME"])/Software/ROOT/install"
       prod_basedir = "$(ENV["LUSTRE"])/CProfiles/"

       verbose           = true
       nevents_per_task  = 200
       nsubtasks         = 100
       ntasks_per_job    = 10
       base_mac          = abspath("templates/cprofile_base.mac")
       config_mac        = abspath("templates/cprofile_config.mac")
       task_template     = abspath("templates/task_template.sh")
       job_template      = abspath("templates/job_template.sh")

       config_variables  = Dict( "energy"   => range(100, 1000, step=10)
                               , "particle" => ["e-", "mu-"])

       queue_command  = pipeline(`squeue -ah`, `wc -l`)
       max_jobs_queue = 100
end

# dont need to modify, just needed by the scripts
sorter(fname) = [parse(Float64, m.match) for m in eachmatch(r"\d+(\.\d+)?", basename(fname))]
