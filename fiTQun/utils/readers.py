import uproot
import numpy  as np
import pandas as pd


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
