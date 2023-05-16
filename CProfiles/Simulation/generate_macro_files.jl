include("config.jl")

out_mac_fname  = joinpath("$prod_basedir/mac/$particle", "cprofile_energyMeV.mac")
out_root_fname = joinpath("$prod_basedir/out/$particle", "cprofile_particle_energyMeV_tag.root")

mkpath(dirname(out_mac_fname))
mkpath(dirname(out_root_fname))

for energy in energies

    out_root_fname_ = replace( out_root_fname
                             , "particle" => particle
                             , "energy"   => energy
                             , "tag"      => tag)

    mac = "/control/execute   $(abspath(base_mac)) \n\
           /WCSim/random/seed $energy              \n\
           /gun/particle      $particle            \n\
           /gun/energy        $energy MeV          \n\
           /WCSimIO/RootFile  $out_root_fname_     \n\
           /run/beamOn        $nevents"

    write(replace(out_mac_fname, "energy" => energy), mac)
end

