import sys
import re
import glob
import argparse
import itertools
import uproot
import ROOT
import numpy  as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from os.path   import expandvars, join, basename


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "create 2D histos and unhitP"
                                    , description = "Compute 2D histogram (predicted vs measured charge) and predicted Q vs unhitP"
                                    , epilog      = """""")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument(   "--infile",   type=str, nargs="?", help = "charge 2D file", default="charge2D_and_unhit.root")
    parser.add_argument( "--wcsimlib",   type=str, nargs="?", help = "WCSim lib path", default="$HOME/Software/WCSim/install/lib")
    
    args = parser.parse_args()
    ##########################################

    # Load WCSimRoot.so
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    # in what follows, "mu" stands for predicted and "q" for measured charge

    infilename = args.infile
    npars      = [5] # args.npars 
    nqranges   = len(npars)
    qranges    = 0   # args.qranges

    # read 2D histogram
    with uproot.open(infilename) as file: 
        Hq, mubins, qbins = file["charge2D"].to_numpy()

    # hardcoded
    qranges = np.linspace(qbins[0], qbins[-1], nqranges)

    # compute bin centers
    mus = (mubins[1:] + mubins[:-1])/2. 
    qs  = ( qbins[1:] +  qbins[:-1])/2.

    # output file
    outfilename = "fitted_cpdf.root"
    fout = ROOT.TFile(outfilename, "RECREATE")

    # define q fit ranges and save to output
    qranges  = np.linspace(qbins[0], qbins[-1], nqranges+1)
    hCPDFrange = ROOT.TH1D("hCPDFrange_type0", "", nqranges, qranges)
    for i, n in enumerate(npars): hCPDFrange.SetBinContent(i, n)
    fout.WriteObject(hCPDFrange, "hCPDFrange_type0")

    # output file
    fout = ROOT.TFile("fitted_cpdf.root", "RECREATE")

    # define output graphs:
    # - gParams: for each range and parameter in this range: q vs parameter
    # - gmuthrs : for each range and parameter in this range: q vs low and q vs high thresholds
    for rang in range(nqranges):
        # select qs in range
        qmin, qmax = qranges[rang], qranges[rang+1]
        selqbins = np.argwhere((qmin <= qbins) & (qbins <= qmax)).flatten()[:-1]

        # get polynomia degree
        deg = npars[rang]-1

        # graphs to fill for this range
        gParams = [ROOT.TGraph() for par in range(npars[rang])] # one graph per parameter
        gmuthrs = [ROOT.TGraph() for i in range(2)] # lower and upper bounds

        # fit for each q in range
        print("Number qbins:", len(qbins)-1)
        for pi, qindex in enumerate(selqbins, 0): # pi stands for point-index
            # get the qvalue, define and fill thresholds
            q = qs[qindex]
            # get pdf values for this q
            pdf    = Hq[:, qindex]
            # remove zeros (to avoid inf when doing the log)
            sel = ~(pdf == 0)
            pdf  = pdf[sel]
            mus_ = mus[sel]
            # fit log(pdf) instead of pdf
            pdflog = np.log(pdf)

            print(f"Processing bin {pi+1}, q = {q}")

            # define mu thresholds (what are them??)
            mulow = q - 4.*np.sqrt(q)
            muup  = q + 4.*np.sqrt(q)
            if mulow<mus_[0] : mulow = mus_[0]
            if muup >mus_[-1]: muup  = mus_[-1]
            gmuthrs[0].SetPoint(pi, q, mulow)
            gmuthrs[1].SetPoint(pi, q, muup)

            # print(mulow, "<", q, "<", muup)

            # fit
            # sel = (mulow<=mus_) & (mus_<=muup)
            # if np.any(sel): pars = np.polyfit(mus_[sel], pdflog[sel], deg)
            # else          : pars = np.polyfit(mus_, pdflog, deg)
            pars = np.polyfit(mus_, pdflog, deg)

            # save fitted parameters
            for par, gPar in zip(pars, gParams): gPar.SetPoint(pi, q, par)

            # save to output file
            for  pari, gPar   in enumerate(gParams): fout.WriteObject(  gPar, f"gParam_type0_Rang{rang}_{pari}")
            for bound, gmuthr in enumerate(gmuthrs): fout.WriteObject(gmuthr, f"gmuthr_type0_Rang{rang}_{bound}")
            
            if (pi+1) == 20:
                plt.ion()

                plt.figure()
                plt.title(f"q = {q}")
                
                plt.scatter(mus_, pdflog, color="k")
                x = np.linspace(0, np.max(mus_), 1000)
                # x = x[(mulow<=x) & (x<=muup)]
                poly = np.poly1d(pars)
                plt.plot(x, poly(x), color="r")

                plt.axvline(mulow, color="k")
                plt.axvline( muup, color="k")

                plt.xlabel(r"$\mu$")
                plt.ylabel("log(pdf)")
                plt.draw()

                input("Press Enter to continue...")
                plt.ioff()
                plt.close()
                break

    # copy PUnhit to output file (TGraph copy not implemented in uproot)
    # fout["PUnhit"]           = fin["PUnhit"]
    fout.Close()

    # copy 2D charge histogram and PUnhit parameters to output file
    with (uproot.open(infilename) as fin, uproot.update(outfilename) as fout):
        fout["charge2D"]         = fin["charge2D"]
        fout["hPunhitPar_type0"] = fin["hPunhitPar_type0"]

    return 


if __name__ == "__main__":
    main()