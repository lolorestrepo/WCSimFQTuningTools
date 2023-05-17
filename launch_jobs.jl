jobsdir        = "$(ENV["LUSTRE"])/CProfiles/mjobs/mu-"
queue_command  = `squeue -ah`
max_jobs_queue = 100
include("./CProfiles/Simulation/config.jl")
using .MyConfig: get_energy_and_index as sorter

get_njobs_in_queue() = parse(Int, replace(readchomp(pipeline(queue_command, `wc -l`)), " " => ""))

job_files = readdir(jobsdir, join=true)
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
