include("config.jl")
using .MyConfig

mjobs_fname = joinpath("$prod_basedir/$particle/mjobs", "mjob_jobid.sh")
job_files   =  readdir("$prod_basedir/$particle/jobs", join=true)
sort!(job_files, by=get_energy_and_index)

mkpath(dirname(mjobs_fname))
mkpath("$prod_basedir/$particle/logs/")
mkpath("$prod_basedir/$particle/errs/")

jobtemplate = 
let fin = open(abspath(mjob_template), "r")
    read(fin, String)
end

jobtemplate =
replace( jobtemplate
       , "NTASKS" => ntasks_per_job)

for mjobfiles in Iterators.partition(job_files, ntasks_per_job)

    energy, idx = get_energy_and_index(mjobfiles[1])
    jobid = join([string(energy), string(Int(idx))], "_")

    JOBS = ""
    for f in mjobfiles
        JOBS *= "srun -n 1 $f &\n"
    end
    JOBS = chopsuffix(JOBS, "&\n")
    JOBS *= "\nwait"

    job = 
    replace( jobtemplate
           , "JOBNAME"       => "$jobid"
           , "LOGFILENAME"   => "$prod_basedir/$particle/logs/$jobid.log"
           , "ERRORFILENAME" => "$prod_basedir/$particle/errs/$jobid.err"
           , "jobid"         => jobid
           , "JOBS"          => JOBS)

    write(replace(mjobs_fname, "jobid" => jobid), job)
end

