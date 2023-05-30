import sys
import re
import glob
import argparse
import itertools
import uproot
import ROOT
import numpy  as np

from itertools import groupby
from os.path   import expandvars, realpath, join, basename


def histogram2d_to_func(hist, xbins, ybins):
    '''
    Converts a numpy 2D histogram into a function
    Example: to-do
    '''
    def getter(x, y):
        xbin = np.digitize(x, xbins) - 1
        ybin = np.digitize(y, ybins) - 1
        return hist[xbin, ybin]
    return getter

def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "integrate_cprofiles"
                                    , description = "description"
                                    , epilog      = "Text at the bottom of help")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument( "-i", "--infile", type=str, nargs="?", help = "", default="cprofiles_merged.root")

    parser.add_argument(    "r0max", type=float, nargs=1, help = "")
    parser.add_argument(  "nr0bins", type=int  , nargs=1, help = "")
    parser.add_argument( "nth0bins", type=int  , nargs=1, help = "")
    
    args = parser.parse_args()
    ##########################################

    # collect the arguments
    r0max    = args.r0max   [0]
    nr0bins  = args.nr0bins [0]
    nth0bins = args.nth0bins[0]

    # define r0, th0 bining
    r0bins  = np.linspace(0, r0max*nr0bins/(nr0bins-1), nr0bins+1)
    th0bins = np.linspace(-1, 1, nth0bins+1)

    # open input file and get energies
    fin = uproot.open(args.infile)
    energies = np.sort([float(re.findall(r'\d+(?:\.\d+)?', key)[0]) for key in fin.keys() if re.match("g_\d+", key)])
    # energy bins such that energies are the bin centers
    ebins = (energies[1:] + energies[:-1])/2.
    ebins = np.insert(ebins,          0, energies [0] - (ebins [0] - energies [0]))
    ebins = np.insert(ebins, len(ebins), energies[-1] - (ebins[-1] - energies[-1]))

    # define 3D histograms to save the integrals
    for n in range(3):
        locals()[f"th3d_{n}"] = ROOT.TH3D( f"I_{n}", f"I_{n}"
                                         ,       nr0bins,  r0bins
                                         ,      nth0bins, th0bins
                                         , len(energies), ebins)
        if n == 0: continue
        locals()[f"th1d_{n}"] = ROOT.TH1D( f"I_iso_{n}", f"I_iso_{n}"
                                         , len(energies), ebins)

    # compute integrals as a function of (energy, r0, th0)
    if args.verbose: print(f"Performing Cherenkov integrals...")
    for ebin, energy in enumerate(energies, 1):
        if args.verbose: print(f"Integrating {energy} MeV...")

        # get cherenkov profile 2D histogram 
        th2d = fin[f"g_{energy}"]
        g, thbins, sbins = th2d.to_numpy()
        g = histogram2d_to_func(g, thbins, sbins)

        # values of s for evaluation
        s  = (sbins[1:] + sbins[:-1])/2.
        ds = sbins[1] - sbins[0]

        # compute bin centers (avoids errors in computing th)
        th0s = (th0bins[1:] + th0bins[:-1])/2.
        r0s  = ( r0bins[1:] +  r0bins[:-1])/2.

        # loop over th0 and r0 values
        for th0bin, th0 in enumerate(th0s, 1):
            for r0bin, r0 in enumerate(r0s, 1):
                # from simple trigonometry
                r2 = r0**2 + s**2 - 2*r0*s*th0
                th = (r0*th0 - s) / np.sqrt(r2)

                # compute integrals
                i = g(th, s)
                I = [sum(i*s**n)*ds for n in range(3)]

                # save computed values
                for n in range(3): locals()[f"th3d_{n}"].SetBinContent(r0bin, th0bin, ebin, I[n])


    # compute isotropic integrals as a function of energy
    gsthr = []
    if args.verbose: print(f"Performing Isotropic integrals...")
    for ebin, energy in enumerate(energies, 1):
        if args.verbose: print(f"Integrating {energy} MeV...")

        # get cherenkov profile 2D histogram 
        th2d = fin[f"g_{energy}"]
        g, thbins, sbins = th2d.to_numpy()

        # values of s for evaluation
        s  = (sbins[1:] + sbins[:-1])/2.
        ds = sbins[1] - sbins[0]

        # get normalized 1D isotropic projection (s)
        iso = np.sum(g, axis=0) * (thbins[1] - thbins[0])

        # compute integrals
        I = [sum(iso*s**n)*ds for n in range(1, 3)]

        # save computed values
        for n in range(1, 3): locals()[f"th1d_{n}"].SetBinContent(ebin, I[n-1])

        # compute gsthr
        gsthr.append(s[np.cumsum(iso)*ds>0.9][0])

    # open output file to save data
    fout = ROOT.TFile("cprofiles_integrals.root", "RECREATE")

    # write 3D integral histograms to output file
    for n in range(3): fout.WriteObject(locals()[f"th3d_{n}"], f"I_{n}")

    # write 1D integral histograms to output file
    for n in range(1, 3): fout.WriteObject(locals()[f"th1d_{n}"], f"I_iso_{n}")

    # write gsthr
    fout.WriteObject(ROOT.TGraph(len(energies), energies, np.array(gsthr)), "gsthr")

    # copy mean number of photons per energy (gNphot)
    x, y = fin["gNphot"].values()
    g = ROOT.TGraph(len(x), x, y)
    g.SetTitle("Mean number of photons per event")
    fout.WriteObject(g, "gNphot")

    fout.Close()
    fin .close()
    return


if __name__ == "__main__":
    main()