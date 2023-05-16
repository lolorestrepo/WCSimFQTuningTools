include("config.jl")

mjobs_fname = joinpath("$prod_basedir/mjobs/$particle/", "mjob_jobid.sh")
job_files   =  readdir("$prod_basedir/jobs/$particle/", join=true)
sort!(job_files, by=get_fnumber)

mkpath(dirname(mjobs_fname))
mkpath("$prod_basedir/logs/")
mkpath("$prod_basedir/errs/")

jobtemplate = 
let fin = open(abspath(mjob_template), "r")
    read(fin, String)
end

jobtemplate =
replace( jobtemplate
       , "NTASKS" => ntasks_per_job)

for mjobfiles in Iterators.partition(job_files, ntasks_per_job)

    jobid = get_fnumber(mjobfiles[1])

    JOBS = ""
    for f in mjobfiles
        JOBS *= "srun -n 1 $f &\n"
    end
    JOBS = chopsuffix(JOBS, "&\n")
    JOBS *= "\nwait"

    job = 
    replace( jobtemplate
           , "JOBNAME"       => "$jobid"
           , "LOGFILENAME"   => "$prod_basedir/logs/$jobid.log"
           , "ERRORFILENAME" => "$prod_basedir/errs/$jobid.err"
           , "jobid"         => jobid
           , "JOBS"          => JOBS)

    write(replace(mjobs_fname, "jobid" => jobid), job)
end

