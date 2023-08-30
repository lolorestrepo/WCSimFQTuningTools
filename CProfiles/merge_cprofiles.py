import re
import glob
import argparse
import itertools
import uproot
import ROOT
import numpy  as np

from itertools import groupby
from os.path   import expandvars, realpath, join, basename

def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "merge_cprofiles"
                                    , description = "description"
                                    , epilog      = "Text at the bottom of help")
    
    parser.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument( "indir", type=str, nargs=1, help = "directory containing produced files")
    parser.add_argument(  "type", type=str, nargs=1, help = "type of cherenkov profile, true (tr) or weighted (wt)") # modify this, it can be tr/wt
    
    args = parser.parse_args()
    ##########################################

    # check if type is tr or wt
    if args.type[0] not in ["tr", "wt"]: raise Exception("Type must be one of tr or wt")

    # sorter function
    get_energy_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))

    # get all files in input dir, filter out the '_flat.root' and sort them based on (energy, index)
    infiles = glob.glob(join(args.indir[0], "*"))
    infiles = [f for f in infiles if re.match("cprofile_\d+(?:\.\d+)?MeV_\d+.root", basename(f))]
    infiles = sorted(infiles, key=get_energy_and_index)

    # split input files in energy groups
    energies = np.unique([get_energy_and_index(f)[0] for f in infiles])
    groups = [list(group) for key, group in groupby(infiles, key=lambda x: get_energy_and_index(x)[0])]

    # open output file to write in it
    fout = ROOT.TFile("cprofiles_merged.root", "RECREATE")

    # loop through each energy group
    nphotons = [] # mean nphotons per energy
    sthrs    = [] # s value at which 90% of photons are emitted per energy
    for energy, files in zip(energies, groups):

        if args.verbose: print(f"Processing {energy} MeV files...".ljust(50))

        nevents = 0
        # start with array to accumulate histogram entries
        with uproot.open(files[0]) as fin:
            nevents += fin["nevents"].counts()[0]
            h, thbins, sbins = fin[f"{args.type[0]}g"].to_numpy()
        # loop through the rest of files
        for file in files[1:]:
            with uproot.open(file) as fin:
                nevents += fin["nevents"].counts()[0]
                h_, _, _ = fin[f"{args.type[0]}g"].to_numpy()
                h += h_ # accumulate entries
                
        # save mean nphotons for this energy
        nphotons.append(h.sum()/nevents)

        # create, fill and save 2D histogram for this energy
        thbinw = thbins[1] - thbins[0]
        sbinw =   sbins[1] -  sbins[0]
        h *= (1/h.sum())*(1/(thbinw*sbinw)) # normalize histogram

        th2d = ROOT.TH2D( f"g_{energy}", f"g_{energy}"
                        , len(thbins)-1, thbins
                        , len(sbins) -1, sbins)
        for ix, iy in itertools.product(range(1, len(thbins)), range(1, len(sbins))): th2d.SetBinContent(ix, iy, h[ix-1, iy-1])
        fout.WriteObject(th2d, f"g_{energy}")

        # compute s distance at which 90% of photons were emitted
        # (needed in the parabolic approximation for the angular response)
        sproj = h.sum(axis=0)*thbinw # projected pdf in s
        
        assert len(sproj) == len(sbins) - 1

        scum = np.cumsum(sproj)*sbinw
        sindex = np.argwhere(scum>0.9).flatten()[0]
        sthrs.append((sbins[sindex] + sbins[sindex+1])/2.)


    # save mean number of photons per energy
    g = ROOT.TGraph(len(energies), energies, np.array(nphotons))
    g.SetTitle("Mean number of photons per event")
    fout.WriteObject(g, "gNphot")

    # save s thresholds
    g = ROOT.TGraph(len(energies), energies, np.array(sthrs))
    g.SetTitle("Length of track")
    fout.WriteObject(g, "gsthr")

    # close output file
    fout.Close()
    return


if __name__ == "__main__":
    main()