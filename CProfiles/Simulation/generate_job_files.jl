include("config.jl")
using .MyConfig

jobs_fname  = joinpath("$prod_basedir/jobs/$particle/", "job_jobid.sh")
macro_files =  readdir("$prod_basedir/mac/$particle/", join=true)
sort!(macro_files, by=get_energy_and_index)

mkpath(dirname(jobs_fname))

jobtemplate = 
let fin = open(abspath(job_template), "r")
    read(fin, String)
end

jobtemplate =
replace( jobtemplate
       , "PROD_BASEDIR"       => abspath(prod_basedir)
       , "FQTUNER_INSTALLDIR" => abspath(fqtunerdir)
       , "WCSIM_INSTALLDIR"   => abspath(wcsimdir)
       , "G4_INSTALLDIR"      => abspath(g4dir)
       , "CONDA_INSTALLDIR"   => abspath(condadir)
       , "ROOT_INSTALLDIR"    => abspath(rootdir))

for macrofile in macro_files

    jobid = join(get_energy_and_index(macrofile), "_")

    job = 
    replace(jobtemplate
           , "jobid"     => jobid
           , "macrofile" => macrofile)

    job_fname_ = replace(jobs_fname, "jobid" => jobid)
    write(job_fname_, job)
    chmod(job_fname_, 0o744)
end

