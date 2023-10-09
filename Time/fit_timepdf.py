import os
import re
import argparse
import warnings
import uproot
import ROOT
import numpy  as np
from scipy.optimize import curve_fit

from os.path import expandvars

# register RankWarning as exception
warnings.filterwarnings("error", category=np.RankWarning)

# define gaussian pdf
def gaussian(x, mu, sig):
    return (1/(sig*np.sqrt(2.*np.pi)))*np.exp(-(x-mu)**2/(2*sig**2))


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = ""
                                    , description = ""
                                    , epilog      = """""")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument(  "--infile",   type=str, nargs="?", help = "time 2D file", default="tres_trueq_2Dhistogram.root")

    parser.add_argument( "--npars_gauss", type=int  , nargs="+", help = "number of parameters", default=[7])
    parser.add_argument( "--npars_pars" , type=int  , nargs="+", help = "number of parameters", default=[7])
    
    args = parser.parse_args()
    ##########################################

    # Config variables
    infilename  = args.infile
    npars_gauss = args.npars_gauss[0]
    npars_pars  = args.npars_pars [0]

    # read file containing 2D histograms for each momentum
    f  = uproot.open(infilename)
    # get momentum values and sort them
    histonames = list(f.classnames().keys())
    momenta = []
    for name in histonames:
        m = re.match(r"htimepdf_(\d+\.\d+)", name)
        if m: momenta.append(float(m.group(1)))
    momenta = np.array(momenta)
    momenta.sort() # inplace sort

    # define momentum bins, such that each value in momenta is the bin center
    pbins = (momenta[:-1] + momenta[1:])/2.
    pbins = np.insert(pbins,            0, momenta[0]  - (pbins[0]  - momenta[0]))
    pbins = np.insert(pbins, len(momenta), momenta[-1] - (pbins[-1] - momenta[-1]))

    # define outfile
    fout = ROOT.TFile("fitted_timepdf.root", "RECREATE")

    # read one of the 2D histograms to get the binning and define bin centers
    _, tresbins, trueqbins = f[f"htimepdf_{momenta[0]}"].to_numpy()
    tress  = (tresbins [1:] + tresbins [:-1])/2.
    trueqs = (trueqbins[1:] + trueqbins[:-1])/2.

    # define TH2D of p and trueq to save fit mean and sigma
    th2ds = [ ROOT.TH2D( "mean",  "mean", len(pbins)-1, pbins, len(trueqbins)-1, trueqbins)
            , ROOT.TH2D("sigma", "sigma", len(pbins)-1, pbins, len(trueqbins)-1, trueqbins)]

    # First loop in momenta: Gaussian fits, for each (p, trueq)
    if args.verbose: print("1) Performing gaussian fits...")
    for pi, p in enumerate(momenta, 1):
        # read 2D histogram for this momentum
        H, _, _ = f[f"htimepdf_{p}"].to_numpy()

        # loop in true charge
        for tqi, tq in enumerate(trueqs, 1):
            # get residual time projection
            proj = H[:, tqi-1]
            norm = np.sum(proj * (tresbins[1:] - tresbins[:-1]))

            # do the fit (only if number of entries is enought)
            if (proj.sum()<45): pars = [np.nan, np.nan]
            else:
                try:
                    pars, cov = curve_fit(gaussian, tress, proj/norm, p0=[0, 0.1], method="trf")
                    if ((pars[0]<tresbins[0])&(tresbins[-1]<pars[0])): pars = [np.nan, np.nan]
                except: pars = [np.nan, np.nan] # if RuntimeError: minimization fails

            # save result
            th2ds[0].SetBinContent(pi, tqi, pars[0])
            th2ds[1].SetBinContent(pi, tqi, pars[1])

    # save gaussian fit results
    fout.WriteObject(th2ds[0], "mean")
    fout.WriteObject(th2ds[1], "sigma")
    if args.verbose: print("==> 1) gaussian fits finished")

    # Second loop in momenta: Polynomial fits, for each p
    if args.verbose: print("2) Performing polynomial fits for each momentum...")
    graphs = [[ROOT.TGraph(len(momenta)) for i in range(npars_gauss)] for i in range(2)] # graphs for mean and sigma
    for i, graph in enumerate(graphs):
        for pi, p in enumerate(momenta, 1):
            # get parameter (mean or sigma)
            gausspar = np.array([th2ds[i].GetBinContent(pi, j) for j in range(1, len(trueqs)+1)])

            # polynomial fit
            sel = ~np.isnan(gausspar)
            pars = np.flip(np.polyfit(trueqs[sel], gausspar[sel], npars_gauss-1))
            for j in range(npars_gauss): graph[j].SetPoint(pi-1, p, pars[j])

    # save parameters to output file
    for i in range(npars_gauss):
        fout.WriteObject(graphs[0][i], f"gtcmnpar_{i}")
        fout.WriteObject(graphs[1][i], f"gtcsgpar_{i}")
    if args.verbose: print("==> 2) polynomial fits for each momentum finished")


    # Second loop in momenta: Polynomial fits, for each p
    if args.verbose: print("3) Performing polynomial fits for each parameter...")
    tpdfpars = [ ROOT.TH2D("htpdfparmn", "Time pdf mean fit parameters" , npars_gauss, 0, npars_gauss, npars_pars, 0, npars_pars)
               , ROOT.TH2D("htpdfparsg", "Time pdf sigma fit parameters", npars_gauss, 0, npars_gauss, npars_pars, 0, npars_pars)]
    
    for i, tpdfpar in enumerate(tpdfpars):
        for j in range(npars_gauss):
            # get fitted parameters
            fitted_pars = np.array([graphs[i][j].GetPointY(k) for k in range(len(momenta))])

            # fit parameter vs momentum
            npars = npars_pars+1
            while (npars>2):
                npars-=1
                if (npars < 2): raise Exception("Cannot perform 3) polynomial fits for each momentum")
                try:
                    pars = np.flip(np.polyfit(momenta, fitted_pars, npars-1))
                    break
                except np.RankWarning: continue
            
            # fill histogram with parameter values
            for k in range(npars)            : tpdfpar.SetBinContent(j+1, k+1, pars[k])
            for k in range(npars, npars_pars): tpdfpar.SetBinContent(j+1, k+1, 0)

    # save parameters to output file
    fout.WriteObject(tpdfpars[0], "htpdfparmn")
    fout.WriteObject(tpdfpars[1], "htpdfparsg")
    if args.verbose: print("==> 3) polynomial fits for each parameter finished")
        
    # Fill info histogram and write it
    htpdfinfo = ROOT.TH1D("htpdfinfo","", 10, 0, 10)
    htpdfinfo.SetBinContent(1, npars_gauss)
    htpdfinfo.SetBinContent(2, trueqs[0])
    htpdfinfo.SetBinContent(3, trueqs[-1])
    htpdfinfo.SetBinContent(4, npars_pars)
    htpdfinfo.SetBinContent(5, momenta[0])
    htpdfinfo.SetBinContent(6, momenta[-1])
    htpdfinfo.SetBinContent(7, 0.) # momentum offset
    fout.WriteObject(htpdfinfo, "htpdfinfo")

    fout.Close()

if __name__ == "__main__":
    main()