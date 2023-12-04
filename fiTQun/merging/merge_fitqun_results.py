import sys
import os
import re
import glob
import argparse
import uproot
import ROOT
import numpy  as np
import pandas as pd

from itertools   import groupby
from os.path     import expandvars, join, basename
from wcsimreader import utils as wcr


def read_time_window_data(filename):

    columns = ["Event", "ntws", "cluster", "t0", "t1", "pt", "px", "py", "pz", "nsub", "t", "G"]
    df = pd.DataFrame(columns=columns)
    
    f = uproot.open(filename)
    t = f["fiTQun"]
    fqntwnd          = t["fqntwnd"]         .array().to_numpy()-1
    fqtwnd_iclstr    = t["fqtwnd_iclstr"]   .array().to_numpy()
    fqtwnd_npeak     = t["fqtwnd_npeak"]    .array().to_numpy()
    fqtwnd_prftt0    = t["fqtwnd_prftt0"]   .array().to_numpy()
    fqtwnd_prftpos   = t["fqtwnd_prftpos"]  .array().to_numpy()
    fqtwnd           = t["fqtwnd"]          .array().to_numpy()
    fqtwnd_peakt0    = t["fqtwnd_peakt0"]   .array().to_numpy()
    fqtwnd_peakiness = t["fqtwnd_peakiness"].array().to_numpy()

    nevents = len(fqntwnd)

    for ev, event in enumerate(range(nevents)):
        nc = fqntwnd[ev]
        for cluster in fqtwnd_iclstr[ev]:
            ts   = fqtwnd        [ev, cluster]
            pt   = fqtwnd_prftt0 [ev, cluster]
            ppos = fqtwnd_prftpos[ev, cluster]
            nsub = fqtwnd_npeak  [ev, cluster]
            for subp in range(nsub):
                t = fqtwnd_peakt0   [ev, cluster, subp]
                G = fqtwnd_peakiness[ev, cluster, subp]
                df.loc[len(df)] = (event, nc, cluster, *ts, pt, *ppos, nsub, t, G)
    f.close()
    return df


def read_subevent_data(filename):

    columns = ["Event", "nsub", "tw", "peak", "N", "Q", "Q0R", "nll0R", "n50", "q50"]
    df = pd.DataFrame(columns=columns)
    
    f = uproot.open(filename)
    t = f["fiTQun"]

    fqnse     = t["fqnse"]    .array().to_numpy()
    fqitwnd   = t["fqitwnd"]  .array().to_numpy()
    fqipeak   = t["fqipeak"]  .array().to_numpy()
    fqnhitpmt = t["fqnhitpmt"].array().to_numpy()
    fqtotq    = t["fqtotq"]   .array().to_numpy()
    fq0rtotmu = t["fq0rtotmu"].array().to_numpy()
    fq0rnll   = t["fq0rnll"]  .array().to_numpy()
    fqn50     = t["fqn50"]    .array().to_numpy()
    fqq50     = t["fqq50"]    .array().to_numpy()

    nevents = len(fqnse)

    for ev, event in enumerate(range(nevents)):
        nsub = fqnse[ev]
        for subp in range(nsub):
            tw    = fqitwnd  [ev, subp]
            peak  = fqipeak  [ev, subp]
            N     = fqnhitpmt[ev, subp]
            Q     = fqtotq   [ev, subp]
            Q0R   = fq0rtotmu[ev, subp]
            nll0R = fq0rnll  [ev, subp]
            n50   = fqn50    [ev, subp]
            q50   = fqq50    [ev, subp]
            df.loc[len(df)] = (event, nsub, tw, peak, N, Q, Q0R, nll0R, n50, q50)
    f.close()
    return df


def read_1Ring_data(filename, pids, event_counter=0):

    columns = ["Event", "peak", "pid", "pc", "p", "t", "x", "y", "z", "theta", "phi", "Q1R", "nll1R", "L", "Eloss"]
    df = pd.DataFrame(columns=columns)
    
    f = uproot.open(filename)
    t = f["fiTQun"]

    fqnse   = t["fqnse"]  .array()
    fqipeak = t["fqipeak"].array()
    nevents = len(fqnse)

    fq1rpcflg = t["fq1rpcflg"].array()
    fq1rtotmu = t["fq1rtotmu"].array()
    fq1rnll   = t["fq1rnll"]  .array()
    fq1rmom   = t["fq1rmom"]  .array()
    fq1rt0    = t["fq1rt0"]   .array()
    fq1rpos   = t["fq1rpos"]  .array()
    fq1rdir   = t["fq1rdir"]  .array()
    fq1rdconv = t["fq1rdconv"].array()
    fq1reloss = t["fq1reloss"].array()

    for ev in range(nevents):
        event = event_counter + ev
        nsub = fqnse[ev]
        for subp in range(nsub):
            peak = fqipeak [ev, subp]
            for pid in pids:
                pc    = fq1rpcflg[ev, subp, pid]
                mom   = fq1rmom  [ev, subp, pid]
                t     = fq1rt0   [ev, subp, pid]
                pos   = fq1rpos  [ev, subp, pid]
                dir   = fq1rdir  [ev, subp, pid]
                Q1R   = fq1rtotmu[ev, subp, pid]
                nll1R = fq1rnll  [ev, subp, pid]
                L     = fq1rdconv[ev, subp, pid]
                Eloss = fq1reloss[ev, subp, pid]

                theta = np.arccos(dir[2])
                phi   = np.arctan2(dir[1], dir[0])
                if (np.sign(phi)<0): phi += 2.*np.pi

                df.loc[len(df)] = (event, peak, pid, pc, mom, t, *pos, theta, phi, Q1R, nll1R, L, Eloss)
    f.close()
    return df



def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "merge_cprofiles"
                                    , description = "description"
                                    , epilog      = "Text at the bottom of help")
    
    parser.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument( "indir", type=str, nargs=1, help = "input directory")
    parser.add_argument("outdir", type=str, nargs=1, help = "output directory")
    parser.add_argument("--wcsimlib", type=str, nargs="?", help = "WCSim lib path", default="$HOME/Software/WCSim/install/lib")
    
    args = parser.parse_args()
    ##########################################

    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    indir  = expandvars(args.indir [0])
    outdir = expandvars(args.outdir[0])
    os.makedirs(outdir, exist_ok=True)

    get_momentum_and_index = lambda fname: list(map(float, re.findall(r'\d+(?:\.\d+)?', basename(fname))))[:2]
    
    # group files by momentum
    filenames = glob.glob(join(indir, "*"))
    momenta   = np.unique([get_momentum_and_index(fname)[0] for fname in filenames])
    filenames = sorted(filenames, key=get_momentum_and_index)
    grouped_filenames = [list(group) for key, group in groupby(filenames, key=lambda x: get_momentum_and_index(x)[0])]

    # momentum fit
    for momentum, filenames in zip(momenta, grouped_filenames):
        if args.verbose: print("Processing momentum ", momentum)
        fits = pd.DataFrame()
        event_counter = 0
        for filename in filenames:
            # tw = read_time_window_data(filename)
            # se = read_subevent_data   (filename)
            R1 = read_1Ring_data      (filename, pids=[1, 2, 3], event_counter=event_counter)
            event_counter += len(R1.Event.unique())
            fits = pd.concat((fits, R1))
        fits.to_csv(join(outdir, f"fq_{momentum}.csv"), index=False)

    return


if __name__ == "__main__":
    main()
    