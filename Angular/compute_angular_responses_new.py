import sys
import glob
import itertools
import argparse
import uproot
import ROOT
import numpy  as np

from os.path   import expandvars, join, basename

sys.path.append("../STable/")
from STable_tools import read_wcsim_geometry, split_tubeids


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "compute Angular responses"
                                    , description = "description"
                                    , epilog      = """""")
    
    parser.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument(      "indir",   type=str,   nargs=1, help = "directory containing produced files")
    parser.add_argument(    "--nbins",   type=int,   nargs=2, help = "histogram bins for the r and eta variables", default=[10, 10])
    parser.add_argument(    "--vaxis",   type=int, nargs="?", help = "detector vertical axis (0)", default=2)
    parser.add_argument(    "--rmax", type=float, nargs="?", help = "histogram edge for radial (r) dimensions in cm"  , default=172.05-20.) # default wcte
    parser.add_argument( "--wcsimlib",   type=str, nargs="?", help = "WCSim lib path"         , default="$HOME/Software/WCSim/install/lib")
    
    args = parser.parse_args()
    ##########################################

    # Load WCSimRoot.so
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    # get simulation filenames, removing those ending in '_flat.root'
    infiles = glob.glob(join(args.indir[0], "*"))
    infiles = sorted([f for f in infiles if "_flat.root" not in basename(f)])

    # Select rotation matrix based on vertical axis
    vaxis = args.vaxis
    if   vaxis == 0:
        R = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    elif vaxis == 1:
        R = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    elif vaxis == 2: pass
    else: raise Exception("Invalid value for vertical axis (vaxis)")

    # read PMT information (ids, positions, orientations)
    _, pmts_df = read_wcsim_geometry(infiles[0])
    pmts_df = pmts_df.set_index("TubeNo")
    
    # get ids of the tubes normal to the vessel (assumed to be the central PMTs)
    inormal = pmts_df[pmts_df.mPMT_PMTNo == 1].index.values.astype(int)
    
    ibottom, itop, iside = split_tubeids(infiles[0], vaxis=vaxis)
    inormal_bottom = ibottom[np.isin(ibottom, inormal)]
    inormal_top    = itop   [np.isin(   itop, inormal)]
    inormal_side   = iside  [np.isin(  iside, inormal)]

    pmtpos = pmts_df.loc[inormal, ("Position_x0", "Position_x1", "Position_x2")]*0.1
    pmtpos.rename({"Position_x0": 0, "Position_x1": 1, "Position_x2": 2}, axis=1, inplace=True)
    if vaxis != 2: pmtpos = np.matmul(R, pmtpos.T).T

    # assume detector radius and lenght are given by extreme pmt position
    det_rad  = np.max(np.sum(pmtpos.loc[:, (0, 1)]**2, axis=1)**0.5)
    det_hlen = np.max(pmtpos.loc[:, 2])
    rmax     = args.rmax

    if args.verbose: 
        print("Estimated detector radius:", det_rad)
        print("Estimated detector half-length:", det_hlen)
    assert (det_rad  > rmax)
    assert (det_hlen > rmax)

    # select pmts with rmax shell
    # side
    sel = (rmax + abs(pmtpos.loc[inormal_side, 2]) <= det_hlen)
    if sel.sum() == 0: raise Exception("rmax is too large")
    inormal_side = inormal_side[sel]
    # top
    sel = (rmax + np.sum(pmtpos.loc[inormal_top, (0, 1)]**2, axis=1)**0.5 <= det_rad)
    if sel.sum() == 0: raise Exception("rmax is too large")
    inormal_top = inormal_top[sel]
    # bottom
    sel = (rmax + np.sum(pmtpos.loc[inormal_bottom, (0, 1)]**2, axis=1)**0.5 <= det_rad)
    if sel.sum() == 0: raise Exception("rmax is too large")
    inormal_bottom = inormal_bottom[sel]

    # redefine inormal
    inormal = inormal[np.isin(inormal, inormal_side)|
                      np.isin(inormal, inormal_top) |
                      np.isin(inormal, inormal_bottom)]
    
    # define histograms for 2D histogram
    rbins     = np.linspace(0, rmax, args.nbins[0]+1)
    etabins   = np.linspace(0,    1, args.nbins[1]+1)
    
    # empty histograms to be filled
    Hside   = np.zeros((len(inormal_side)  , len(rbins)-1, len(etabins)-1))
    Htop    = np.zeros((len(inormal_top)   , len(rbins)-1, len(etabins)-1))
    Hbottom = np.zeros((len(inormal_bottom), len(rbins)-1, len(etabins)-1))
    Hall    = np.zeros((len(inormal)       , len(rbins)-1, len(etabins)-1))

    # Loop over simulation files and fill the histograms
    for n, filename in enumerate(infiles, 1):
        if args.verbose: print(f"Processing file {n}/{len(infiles)}...", basename(filename))
        # read data
        with uproot.open(filename) as f:
            # contains the following data if the optical photon reaches a PMT
            # either directly or indirectly
            # photon emission position
            oppos   = f["sttree/oppos"]  .array().to_numpy()*0.1
            # PMT number/pos
            ihPMT   = f["sttree/ihPMT"]  .array().to_numpy()
            tubepos = f["sttree/tubepos"].array().to_numpy()*0.1
            # ISTORY=(# of refl)*1000 + (# of scat); 0 if direct hit
            isct = f["sttree/isct"].array().to_numpy()

        # select direct photons
        sel_direct = (isct == 0)
        oppos   = oppos  [sel_direct]
        ihPMT   = ihPMT  [sel_direct]
        tubepos = tubepos[sel_direct]

        # select only PMTs normal to the vessel (central PMTs) and with rmax shell
        sel_normal = np.isin(ihPMT, inormal)
        oppos   = oppos  [sel_normal]
        ihPMT   = ihPMT  [sel_normal]
        tubepos = tubepos[sel_normal]

        # get PMT orientation
        tubedir = pmts_df.loc[ihPMT, ["Orientation_x0", "Orientation_x1", "Orientation_x2"]].values.astype(float)

        # compute variables of interest:
        #    - distance between PMT and emission point (mm to cm)
        #    - cosine of the angle of incidence w.r.t. the PMT normal
        # don't need to rotate since module and angles are invariant
        r    = oppos - tubepos
        modr = np.sqrt(np.sum(r**2, axis=1))
        
        # compute angles
        eta  = np.array([np.dot(v1, v2) for (v1, v2) in zip(r, tubedir)])/modr

        # fill histograms 
        for pmtid in inormal:
            sel = (ihPMT == pmtid)
            h, _, _ = np.histogram2d(modr[sel], eta[sel], bins=[rbins, etabins])

            index = np.argwhere(inormal == pmtid)
            if len(index) == 1:    Hall[index[0, 0]] += h
            else: continue

            index = np.argwhere(inormal_side == pmtid)
            if len(index) == 1:   Hside[index[0, 0]] += h; continue

            index = np.argwhere(inormal_top == pmtid)
            if len(index) == 1:    Htop[index[0, 0]] += h; continue

            index = np.argwhere(inormal_bottom == pmtid)
            if len(index) == 1: Hbottom[index[0, 0]] += h; continue


    # compute efficiencies
    assert all(  Hside.sum(axis=(1, 2)) > 0)
    assert all(   Htop.sum(axis=(1, 2)) > 0)
    assert all(Hbottom.sum(axis=(1, 2)) > 0)

    Hside   =   Hside /   Hside.sum(axis=(1, 2))[:, np.newaxis, np.newaxis]
    Htop    =    Htop /    Htop.sum(axis=(1, 2))[:, np.newaxis, np.newaxis]
    Hbottom = Hbottom / Hbottom.sum(axis=(1, 2))[:, np.newaxis, np.newaxis]
    Hall    =    Hall /    Hall.sum(axis=(1, 2))[:, np.newaxis, np.newaxis]

    # compute mean values
    eff_side   =   Hside.mean(axis=0)
    eff_top    =    Htop.mean(axis=0)
    eff_bottom = Hbottom.mean(axis=0)
    eff_all    =    Hall.mean(axis=0)

    # save photon tables
    if args.verbose: print("Writting histograms...")
    
    # define root 2D histograms
    th2d_side   = ROOT.TH2D("AngularResponse_side", "AngularResponse_side", len(rbins)-1, rbins, len(etabins)-1, etabins)
    th2d_top    = ROOT.TH2D("AngularResponse_top", "AngularResponse_top", len(rbins)-1, rbins, len(etabins)-1, etabins)
    th2d_bottom = ROOT.TH2D("AngularResponse_bottom", "AngularResponse_bottom", len(rbins)-1, rbins, len(etabins)-1, etabins)
    th2d_all    = ROOT.TH2D("AngularResponse_all", "AngularResponse_all", len(rbins)-1, rbins, len(etabins)-1, etabins)

    # fill histograms
    for ix, iy in itertools.product(range(1, len(rbins)), range(1, len(etabins))):
        th2d_side  .SetBinContent(ix, iy,   eff_side[ix-1, iy-1])
        th2d_top   .SetBinContent(ix, iy,    eff_top[ix-1, iy-1])
        th2d_bottom.SetBinContent(ix, iy, eff_bottom[ix-1, iy-1])
        th2d_all   .SetBinContent(ix, iy,    eff_all[ix-1, iy-1])

    # write output
    fout = ROOT.TFile("angular.root", "RECREATE")
    fout.WriteObject(  th2d_side, "AngularResponse_side")
    fout.WriteObject(   th2d_top, "AngularResponse_top")
    fout.WriteObject(th2d_bottom, "AngularResponse_bottom")
    fout.WriteObject(   th2d_all, "AngularResponse_all")
    fout.Close()
    
    return


if __name__ == "__main__":
    main()