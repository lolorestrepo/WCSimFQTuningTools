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
       nevents_per_task  = 50
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

# dont need to modify, it is needed by the generate_job_files.jl and launch_jobs.jl scripts
function sorter(fname)
       # returns the values of the config_variables for a given fname
       m    = match(r"_(.+)\.", basename(fname))
       @assert length(m.captures) == 1
       vars = split(m.captures[1], "_")
       return [(tryparse(Float64, v) === nothing) ? v : parse(Float64, v) for v in vars]
end
