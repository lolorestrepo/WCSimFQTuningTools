jobsdir        = "$(ENV["HOME"])/WCTE/Production/Tuning/CProfiles/Simulation/prod/mjobs/mu-"
queue_command  = `squeue -ah`
max_jobs_queue = 10

get_fnumber(fname::String) = parse(Int, match(r"[0-9]+", basename(fname)).match)
get_njobs_in_queue()       = parse(Int, replace(readchomp(pipeline(queue_command, `wc -l`)), " " => ""))

job_files = readdir(jobsdir, join=true)
sort!(job_files, by=get_fnumber)

for job in job_files

    free_queue = @task while 
        get_njobs_in_queue() > max_jobs_queue
        sleep(10)
    end

    schedule(free_queue)
    wait(free_queue)

    run(`sbatch $job`)
end