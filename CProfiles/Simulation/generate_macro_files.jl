include("config.jl")
using .MyConfig

out_mac_fname  = joinpath("$prod_basedir/mac/$particle", "cprofile_energyMeV_idx.mac")
out_root_fname = joinpath("$prod_basedir/out/$particle", "cprofile_energyMeV_idx_tag.root")

mkpath(dirname(out_mac_fname))
mkpath(dirname(out_root_fname))

global rseed = 1
for energy in energies

    for idx in range(1, ntasks_per_energy)

        out_root_fname_ = replace( out_root_fname
                                 , "energy"   => energy
                                 , "idx"      => idx
                                 , "tag"      => tag)

        mac = "/control/execute   $(abspath(base_mac)) \n\
               /WCSim/random/seed $rseed               \n\
               /gun/particle      $particle            \n\
               /gun/energy        $energy MeV          \n\
               /gun/direction     0 0 1                \n\
               /WCSimIO/RootFile  $out_root_fname_     \n\
               /run/beamOn        $nevents_per_task"
        global rseed += 1

        write(replace(out_mac_fname, "energy" => energy, "idx" => idx), mac)
    end
end

