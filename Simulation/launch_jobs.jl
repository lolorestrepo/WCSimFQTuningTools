include("config.jl")
using .Parameters

get_njobs_in_queue() = parse(Int, replace(readchomp(queue_command), " " => ""))

jobs_dir  = abspath(prod_basedir * "/jobs")
job_files = readdir(jobs_dir, join=true)
sort!(job_files, by=sorter)

for job in job_files

    free_queue = @task while 
        get_njobs_in_queue() > max_jobs_queue
        sleep(10)
    end

    schedule(free_queue)
    wait(free_queue)

    run(`sbatch $job`)

    println("$job launched")
end
