import argparse
import uproot
import ROOT
import numpy  as np


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
    parser.add_argument( "-i", "--infile", type=str, nargs="?", help = "", default="cprofiles_integrals.root")

    parser.add_argument( "npars", type=int, nargs=1, help = "")
    
    args = parser.parse_args()
    ##########################################

    npars = args.npars[0]
    # open input and output files
    fin  = uproot.open(args.infile)
    fout = ROOT.TFile("cprofiles_fits.root", "RECREATE")

    # fit parameters are saved in 3D histo, thus define the bins
    # add +2 bins to save the fit bounds (to be fiTQun compatible)
    parsbins = np.arange(-0.5, npars+0.5+2, 1)
    
    # loop over s exponent in the integrals
    for n in range(3):

        if args.verbose: print(f"Fitting n={n}...")

        h, r0bins, th0bins, mbins = fin[f"I_{n}"].to_numpy()

        th3d = ROOT.TH3D( f"Js_{n}", f"Js_{n}"
                        ,  len(r0bins)-1,  r0bins
                        , len(th0bins)-1, th0bins
                        ,        npars+2, parsbins)
        
        th2d_chi2 = ROOT.TH2D( f"chi2_{n}", f"chi2_{n}"
                             ,  len(r0bins)-1,  r0bins
                             , len(th0bins)-1, th0bins)
        
        # This is a dummy ones 2D histogram for compatibility with fiTQun
        th2d_nsect = ROOT.TH2D( f"nsect_{n}", f"nsect_{n}"
                              ,  len(r0bins)-1,  r0bins
                              , len(th0bins)-1, th0bins)

        # fit I_n vs log(p)
        log_momenta = np.log((mbins[1:] + mbins[:-1])/2.)

        # loop over r0 and th0 bins
        for th0bin, _ in enumerate(th0bins[:-1], 0):
            for r0bin, _ in enumerate(r0bins[:-1], 0):
                # get integrals for (r0, th0) bins
                I = h[r0bin, th0bin]
                # polynomial fit
                js = np.polyfit(log_momenta, I, npars-1)

                # fill 3D histogram with results
                for i, ji in enumerate(np.flip(js), 1): th3d.SetBinContent(r0bin+1, th0bin+1, i, ji)
                # save bounds (to be fiTQun compatible)
                th3d.SetBinContent(r0bin+1, th0bin+1, npars+1, log_momenta[0])
                th3d.SetBinContent(r0bin+1, th0bin+1, npars+2, log_momenta[-1])

                # compute and save chi2
                p = np.poly1d(js)
                chi2 = sum((p(log_momenta) - I)**2)
                th2d_chi2.SetBinContent(r0bin+1, th0bin+1, chi2)

                # save 1 in nsect histo (to be fiTQun compatible)
                th2d_nsect.SetBinContent(r0bin+1, th0bin+1, 1)

        fout.WriteObject(      th3d, f"hI3d_par_{n}")
        fout.WriteObject( th2d_chi2, f"hI3d_par_{n}_chi2")
        fout.WriteObject(th2d_nsect, f"hI3d_nsect_{n}")
    
    # Add hprofinf
    hprofinf = ROOT.TH1F("hprofinf", "hprofinf", 6,  np.arange(-0.5, 6.5, 0.5))
    hprofinf.SetBinContent(1, npars)
    hprofinf.SetBinContent(2, 50)   # nsectmax, hardcoded (?), not really used in fiTQun
    hprofinf.SetBinContent(3, 0)    # momentum ofset
    hprofinf.SetBinContent(4, 0.2)  # momentum min step, hardcoded (?)
    hprofinf.SetBinContent(5, np.log(mbins[0]))  # min momentum
    hprofinf.SetBinContent(6, np.log(mbins[-1])) # max momentum
    fout.WriteObject(hprofinf, "hprofinf")
    

    # fit and save isotropic integrals, gNphot and gsthr
    # add +2 bins to save the fit bounds (to be fiTQun compatible)
    parsbins = np.arange(-0.5, npars+0.5+2, 1)
    # fit and save isotopric integrals
    for n in range(1, 3):
        h, mbins = fin[f"I_iso_{n}"].to_numpy()
        momenta = (mbins[1:] + mbins[:-1])/2.
        pars = np.polyfit(np.log(momenta), np.log(h), npars-1)
        th1d = ROOT.TH1D(f"gI_iso_{n}_pars", f"gI_iso_{n}_pars",  len(parsbins)-1, parsbins)
        for i, pi in enumerate(np.flip(pars), 1): th1d.SetBinContent(i, pi)
        th1d.SetBinContent(npars+1, np.log(momenta[0]))
        th1d.SetBinContent(npars+2, np.log(momenta[-1]))
        fout.WriteObject(th1d, f"gI_iso_{n}_pars")

    # mean number of photons per momentum (gNphot)
    # copy 
    x, y = fin["gNphot"].values()
    g = ROOT.TGraph(len(x), x, y)
    g.SetTitle("Mean number of photons per event")
    fout.WriteObject(g, "gNphot")
    # fit, fill histogram and save
    pars = np.polyfit(np.log(x), np.log(y), npars-1)
    th1d = ROOT.TH1D("gNphot_pars", "gNphot_pars",  len(parsbins)-1, parsbins)
    for i, pi in enumerate(np.flip(pars), 1): th1d.SetBinContent(i, pi)
    th1d.SetBinContent(npars+1, np.log(x[0]))
    th1d.SetBinContent(npars+2, np.log(x[-1]))
    fout.WriteObject(th1d, "gNphot_pars")

    # copy gsthr
    # copy
    x, y = fin["gsthr"].values()
    g = ROOT.TGraph(len(x), x, y)
    fout.WriteObject(g, "gsthr")
    # fit, fill histogram and save
    pars = np.polyfit(np.log(x), np.log(y), npars-1)
    th1d = ROOT.TH1D("gsthr_pars", "gsthr_pars",  len(parsbins)-1, parsbins)
    for i, pi in enumerate(np.flip(pars), 1): th1d.SetBinContent(i, pi)
    th1d.SetBinContent(npars+1, np.log(x[0]))
    th1d.SetBinContent(npars+2, np.log(x[-1]))
    fout.WriteObject(th1d, "gsthr_pars")
    fout.Close()

    # copy the original integral histograms to allow easier checking
    with uproot.update("cprofiles_fits.root") as fout:
        for n in range(3)   : fout[f"I_{n}"]      = fin[f"I_{n}"]
        for n in range(1, 3): fout[f"hI_iso_{n}"] = fin[f"I_iso_{n}"]
    fin .close()
    return

if __name__ == "__main__":
    main()