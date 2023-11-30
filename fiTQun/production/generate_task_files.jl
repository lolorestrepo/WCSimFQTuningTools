using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "indir"
            help = "Input directory"
            required = true
        "outdir"
            help = "Output directory"
            required = true
        "parameters"
            help = "Parameters override filename"
            required = true
        "--task_template", "-t"
            help = ""
            arg_type = String
            default  = "task_template.sh"
        "--verbose", "-v"
            help = "verbose"
            action = :store_true
    end
    return parse_args(s)
end

function sorter(fname)
    m    = match(r"_(.+)\.", basename(fname))
    @assert length(m.captures) == 1
    vars = split(m.captures[1], "_")
    return [(tryparse(Int    , v) === nothing) ?
           ((tryparse(Float64, v) === nothing) ? v : parse(Float64, v)) : parse(Int, v) for v in vars]
end

function main()
    
    parsed_args = parse_commandline()

    indir         = abspath(parsed_args["indir"])
    outdir        = abspath(parsed_args["outdir"])
    parfile       = abspath(parsed_args["parameters"])
    task_template = abspath(parsed_args["task_template"])
    verbose       =         parsed_args["verbose"]

    if !isdir(indir)    error("Input directory does not exists") end
    if !isfile(parfile) error("Override parameters file does not exist") end

    if verbose
        println("Input directory: ", indir)
        println("Output directory: ", outdir)
        println("Override parameters filename: ", parfile)
    end 

    tasktemplate = 
    let fin = open(task_template, "r")
        read(fin, String)
    end

    tasktemplate =
    replace( tasktemplate
           , "PARAMETER_FILE"    => parfile)

    inputfiles = readdir(indir, join=true)
    inputfiles = [file for file in inputfiles if (match(r"out_(.+)_(\d*\.?\d+)_(\d+).root", basename(file)) !== nothing)]
    sort!(inputfiles, by=sorter)

    # # group by energies
    # energies = []
    # for infile in inputfiles push!(energies, sorter(infile)[2]) end
    # unique!(energies)

    mkpath(joinpath(outdir, "out"))
    mkpath(joinpath(outdir, "tasks"))

    for infile in inputfiles

        taskid = join(sorter(infile), "_")

        outfile = replace(basename(infile), "out" => "fq")
        outfile = joinpath(outdir, "out", outfile)

        task = replace( tasktemplate
                      , "input_file"  => infile
                      , "output_file" => outfile)

        task_fname = joinpath(outdir, "tasks", "task_$taskid.sh")
        if verbose println("Writting $task_fname") end
        write(task_fname, task)
        chmod(task_fname, 0o744)
    end

end


main()