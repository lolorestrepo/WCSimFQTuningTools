import sys
import re
import glob
import argparse
import uproot
import ROOT
import numpy  as np
import matplotlib.pyplot as plt

from itertools import groupby
from os.path   import expandvars, join, basename


def get_DigiHitQs_and_nHits(filename):
    """Get DigiHit charges and number of hits for each PMT and event in the root file"""
    rootf = ROOT.TFile(filename, "read")
    tree  = rootf.GetKey("wcsimT").ReadObj()
    nevents = tree.GetEntries()

    Qs   = []
    nHit = []
    for event in range(nevents):
        tree.GetEvent(event)
        trigger  = tree.wcsimrootevent.GetTrigger(0)
        # fill DigiHit charges
        DigiHits = trigger.GetCherenkovDigiHits()
        for dhit in DigiHits: Qs.append(dhit.GetQ())
        # fill nHit
        nHit.append(trigger.GetNumTubesHit()) # same as trigger.GetNcherenkovhits()

    return np.array(Qs), np.array(nHit)


def get_nPMTs(filename):
    """Get number of PMTs in the detector"""
    return 0


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "create 2D histos and unhitP"
                                    , description = "Compute 2D histogram (predicted vs measured charge) and predicted Q vs unhitP"
                                    , epilog      = """""")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument( "indir", type=str, nargs=1, help = "directory containing input files")
    parser.add_argument( "--wcsimlib",   type=str, nargs="?", help = "WCSim lib path", default="$HOME/Software/WCSim/install/lib")
    
    args = parser.parse_args()
    ##########################################

    # Load WCSimRoot.so
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    # in what follows, "mu" stands for predicted and "q" for measured charge
    # sorter function
    get_mu_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))[:2]))

    # get all files in input dir, filter out the '_flat.root' and sort them based on (mu, index)
    infiles = glob.glob(join(args.indir[0], "*"))
    infiles = [f for f in infiles if re.match("out_\d+(?:\.\d+)?_\d+.root", basename(f))]
    infiles = sorted(infiles, key=get_mu_and_index)

    # split input files in mu groups
    mus = np.unique([get_mu_and_index(f)[0] for f in infiles])
    groups = [list(group) for key, group in groupby(infiles, key=lambda x: get_mu_and_index(x)[0])]

    # loop over mu and fill histograms for each one
    Hq = [] # 2D histogram q vs mu
    qbins = np.loadtxt("qbins_wcte.txt")
    for mu, files in zip(mus, groups):
        Qs    = []
        nHits = []
        for filename in files:
            Qs_, nHits_ = get_DigiHitQs_and_nHits(filename)
            Qs   .extend(Qs_   )
            nHits.extend(nHits_)
        Qs    = np.array(Qs)
        nHits = np.array(nHits)

        hq, _ = np.histogram(Qs, bins=qbins)
        Hq.append(hq)

        plt.figure()
        plt.title(rf"Predicted charge $\mu$ = {mu}")
        plt.hist(Qs, bins=qbins, histtype="step")
        plt.show(block=True)

        break


    # normalize nHits/nPMTs

    # compute mu bins such that mus array lie at the center
    mubins = (mus[:-1] + mus[1:])/2.
    mubins = np.insert(mubins,           0, mus[0]  - (mus[0]  - mus[0]))
    mubins = np.insert(mubins, len(mubins), mus[-1] - (mus[-1] - mus[-1]))

    return 


if __name__ == "__main__":
    main()