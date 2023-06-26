path      = "$(ENV["LUSTRE"])/CProfiles/out/"
filenames = readdir(path, join=true)


for filename in filenames
    
    m = match(r"out_e-_(\d+\.?\d+)_(\d+).root", basename(filename))

    if m === nothing continue end

    energy  = parse(Float64, m.captures[1])
    subtask = parse(    Int, m.captures[2])

    println(basename(filename), "  --->  ", "cprofile_$(energy)MeV_$subtask.root")

    mv(filename, joinpath(dirname(filename), "cprofile_$(energy)MeV_$subtask.root"))
end