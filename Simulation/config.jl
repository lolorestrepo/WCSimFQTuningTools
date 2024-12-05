module Parameters

       export fqtunerdir, wcsimdir, g4dir, rootdir, prod_basedir
       export verbose, nevents_per_task, nsubtasks, ntasks_per_job
       export base_mac, config_mac, task_template, job_template
       export config_variables, queue_command, max_jobs_queue

       fqtunerdir   = "$(ENV["HOME"])/Software/WCSimFQTuner/install/"
       wcsimdir     = "$(ENV["HOME"])/Software/WCSim/install-Linux_x86_64-gcc_13.2.0-python_3.10.13/"
       g4dir        = "$(ENV["THRONG_DIR"])/Software/Geant4/geant4_v10.5.1/"
       rootdir      = "$(ENV["THRONG_DIR"])/Software/ROOT/root_v6-28-00/"
       prod_basedir = "$(ENV["LUSTRE"])/Pyfitqun/center/e-/"

       verbose           = true
       nevents_per_task  = 100
       nsubtasks         = 20
       ntasks_per_job    = 10
       base_mac          = abspath("templates/gps_base.mac")
       config_mac        = abspath("templates/gps_config.mac")
       task_template     = abspath("templates/task_template.sh")
       job_template      = abspath("templates/job_template.sh")

       # CProfiles/GPS // use momentum for CProfiles and energy for GPS
       config_variables  = Dict( "energy"   => range(150, 600, step=50)
                               , "particle" => ["e-"])

       # # Charge PDF
       # config_variables = Dict("true_charge" => vcat( range(0.1, 1.9, step=0.1)
       #                                              , range(2.0, 4.5, step=0.5)
       #                                              , range(5.0, 9.0, step=1)
       #                                              , range(  10, 19, step=2)
       #                                              , range(  20, 50, step=5)))

       # # STable and Angular
       # config_variables = Dict("particle" => ["e-"])

       # # Time PDF
       # config_variables = Dict("particle" => ["mu-"]
       #                        ,"energy"   => vcat( range(40, 90, step=10)
       #                                           , range(100, 575, step=25) 
       #                                           , range(600, 900, step=50)))

       # # Time PDF
       # config_variables = Dict( "particle" => ["mu-"]
       #                        , "energy"   => range(120, 360, step=50))

       # # Prefit
       # # y vertical distance
       # ymax = 120.
       # rmax = 120.
       # ys     = collect(range(0, ymax, 3))
       # rs     = collect(range(0, rmax, 3))
       # thetas = collect(range(0, pi/2., 3))
       
       # positions = []
       # for y in ys
       #        for r in rs
       #             x = round.(r * cos.(thetas))
       #             z = round.(r * sin.(thetas))
       #             append!(positions, collect(zip(x, repeat([y], length(x)), z)))
       #        end
       # end
       # unique!(positions)
       # positions = [join(position, " ") for position in positions]

       # config_variables = Dict( "particle"  => ["e-"]
       #                        , "energy"    => range(100, 250, step=50)
       #                        , "position"  => positions
       #                        , "direction" => ["1 0 0"])
       
       # # Prefit
       # # y vertical distance
       # ymax = 120.
       # rmax = 120.
       # ys     = collect(range(0, ymax, 3))
       # rs     = collect(range(0, rmax, 3))
       # thetas = collect(range(0, pi/2., 3))

       # ys = [0.0, 50.0]
       # rs = [0.0, 50.0]
       # thetas = [0., pi/2.]
       
       # positions = []
       # for y in ys
       #        for r in rs
       #             x = round.(r * cos.(thetas))
       #             z = round.(r * sin.(thetas))
       #             append!(positions, collect(zip(x, repeat([y], length(x)), z)))
       #        end
       # end
       # unique!(positions)
       # positions = [join(position, " ") for position in positions]

       # config_variables = Dict( "particle"  => ["e-"]
       #                        , "energy"    => [150]
       #                        , "position"  => positions
       #                        , "direction" => ["1 0 0", "-1 0 0"])

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
