import sys
import glob
import argparse
import uproot
import ROOT
import numpy  as np

from os.path   import expandvars, join, basename

from STable_tools import read_wcsim_geometry, split_tubeids, azimuth_angle


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "compute STable"
                                    , description = "description"
                                    , epilog      = """ The table variables are the following:
                                                    zs: source z position
                                                    Rs: source R position
                                                    (bottom and top table) R_PMT: R position of the PMT
                                                    (side table) z_PMT: zposition of the PMT
                                                    phi: angle between source and PMT
                                                    zd : source z direction
                                                    theta: source direction angle w.r.t PMT""")
    
    parser.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument(    "indir", type=str, nargs=1, help = "directory containing produced files")
    parser.add_argument(    "nbins", type=int, nargs=7, help = "histogram bins in the following order for variables: zs,Rs,R_PMT,z_PMT,phi,zd,theta"
                                                      , default=[35, 16, 8, 16, 16, 16, 16])
    parser.add_argument(    "-vaxis",    "--vaxis", type=int, nargs="?", help = "detector vertical axis (0)", default=2)
    parser.add_argument( "-wcsimlib", "--wcsimlib", type=str, nargs="?", help = "WCSim lib path"            , default="$HOME/Software/WCSim/install/lib")
    parser.add_argument(   "-fitqun",   "--fitqun", type=str, nargs="?", help = "fiTQun source directory"   , default="$HOME/Software/fiTQun/fiTQun/")

    parser.add_argument( "-o", "--outpath", type=str, nargs="?", help = ".hdf5 file path", default=".")
    
    args = parser.parse_args()
    ##########################################

    # Load TScatTable class from fiTQun
    if args.verbose: print("Compiling fiTQun TScatTableF.cc...")
    ROOT.gROOT.SetMacroPath(expandvars(args.fitqun))
    ROOT.gROOT.LoadMacro("TScatTableF.cc++")

    # Load WCSimRoot.so
    ROOT.gSystem.AddDynamicPath(expandvars(args.wcsimlib))
    ROOT.gSystem.Load          ("libWCSimRoot.dylib" if sys.platform == "darwin" else "libWCSimRoot.so")

    # get simulation filenames, removing those ending in '_flat.root'
    infiles = glob.glob(join(args.indir[0], "*"))
    infiles = [f for f in infiles if "_flat.root" not in basename(f)]

    # Get horizontal indexes (vaxis defines the index of the vertical direction)
    vaxis = args.vaxis
    axes  = np.array([0, 1, 2])
    haxes = axes[axes != vaxis]

    # get bottom, top and side tube-ids from the first file
    tubeid_bottom, tubeid_top, tubeid_side = split_tubeids(infiles[0], vaxis=vaxis)

    # Define TScatTable 6D histogram
    df, _     = read_wcsim_geometry(infiles[0])
    tube_ztop = df.loc["WCCylLength", "WC"]/2.  # detector cylinder half-length
    tube_rad  = df.loc["WCPMTRadius", "WC"]     # PMT module radius
    cyl_rad   = df.loc["WCCylRadius", "WC"]     # detector cylinder radius

    # Define scatering table binning
    # 6 variables
    # zs: source z position
    # Rs: source R position
    # (bottom and top table) R_PMT
    # (side table) z_PMT: zposition of the PMT
    # phi: angle between source and PMT
    # zd : source z direction
    # theta: source direction angle w.r.t PMT
    nzs    = args.nbins[0]
    nRs    = args.nbins[1]
    nR_PMT = args.nbins[2]
    nz_PMT = args.nbins[3]
    nphi   = args.nbins[4]
    nzd    = args.nbins[5]
    ntheta = args.nbins[6]
    
    zs    = [-(tube_ztop - tube_rad), tube_ztop - tube_rad]
    Rs    = [0., cyl_rad - tube_rad]
    R_PMT = [0, cyl_rad]
    z_PMT = [-(tube_ztop - tube_rad), tube_ztop - tube_rad]
    phi   = [-np.pi, +np.pi]
    zd    = [-1., 1.]
    theta = [-np.pi, +np.pi]

    # Define scattering tables, notice the only change beside name and title is the use of R_PMT (bottom and top) and z_PMT (side)
    # bottom
    Tbottom_indirect = ROOT.TScatTableF( "bottomnphot", "Number of Photons for Tubes on the Bottom of the Tank"
                                       , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nR_PMT, R_PMT[0], R_PMT[1]
                                       , nphi, phi[0], phi[1], nzd, zd[0], zd[1], ntheta, theta[0], theta[1])
    
    Tbottom_direct   = ROOT.TScatTableF( "bottomnphot", "Number of Photons for Tubes on the Bottom of the Tank"
                                       , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nR_PMT, R_PMT[0], R_PMT[1], nphi, phi[0], phi[1])
    
    # top
    Ttop_indirect    = ROOT.TScatTableF( "topnphot", "Number of Photons for Tubes on the Top of the Tank"
                                       , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nR_PMT, R_PMT[0], R_PMT[1]
                                       , nphi, phi[0], phi[1], nzd, zd[0], zd[1], ntheta, theta[0], theta[1])
    
    Ttop_direct      = ROOT.TScatTableF( "topnphot", "Number of Photons for Tubes on the Top of the Tank"
                                       , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nR_PMT, R_PMT[0], R_PMT[1], nphi, phi[0], phi[1])

    # side
    Tside_indirect = ROOT.TScatTableF( "sidenphot", "Number of Photons for Tubes on the Side of the Tank"
                                     , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nz_PMT, z_PMT[0], z_PMT[1]
                                     , nphi, phi[0], phi[1], nzd, zd[0], zd[1], ntheta, theta[0], theta[1])
    
    Tside_direct   = ROOT.TScatTableF( "sidenphot", "Number of Photons for Tubes on the Side of the Tank"
                                     , nzs, zs[0], zs[1], nRs, Rs[0], Rs[1], nz_PMT, z_PMT[0], z_PMT[1], nphi, phi[0], phi[1])

    # Loop over simulation files in order to fill the table histograms
    for filename in infiles:
        if args.verbose: print("Processing", basename(filename))
        # read data
        with uproot.open(filename) as f:
            # contains the following data if the optical photon reaches a PMT
            # either directly or indirectly
            # parent pos and dir
            srcpos = f["sttree/srcpos"].array().to_numpy()
            srcdir = f["sttree/srcdir"].array().to_numpy()
            # photon pos and dir
            oppos  = f["sttree/oppos"] .array().to_numpy()
            opdir  = f["sttree/opdir"] .array().to_numpy()
            # PMT number/pos
            ihPMT   = f["sttree/ihPMT"]  .array().to_numpy()
            tubepos = f["sttree/tubepos"].array().to_numpy()
            # distance covered? (not needed)
            deltaL = f["sttree/deltaL"].array().to_numpy()
            # ISTORY=(# of refl)*1000 + (# of scat); 0 if direct hit
            isct = f["sttree/isct"].array().to_numpy()

        # Compute scattering table variables
        zs    =  srcpos[:, vaxis]
        Rs    = np.sum(srcpos[:, haxes]**2, axis=1)**0.5
        # z_PMT used for side table, R_PMT used for bottom and top
        z_PMT = tubepos[:, vaxis]
        R_PMT = np.sum(tubepos[:, haxes]**2, axis=1)**0.5
        phi   = np.abs(azimuth_angle(srcpos) - azimuth_angle(tubepos))
        zd    = opdir[:, vaxis]
        theta = np.abs(azimuth_angle(srcpos) - azimuth_angle(tubepos-srcpos))

        # transform angles from [0, 2pi] to [-pi, pi]
        phi  [  phi>np.pi] -= 2.*np.pi
        theta[theta>np.pi] -= 2.*np.pi

        # Fill indirect tables
        sel_indirect = (isct != 0) # select indirect photons
        # bottom
        sel = sel_indirect & np.isin(ihPMT, tubeid_bottom)
        for vals in zip(zs[ sel], Rs[ sel], R_PMT[ sel], phi[ sel], zd[ sel], theta[ sel]):Tbottom_indirect.Fill(*vals)
        # top
        sel = sel_indirect & np.isin(ihPMT, tubeid_top)
        for vals in zip(zs[ sel], Rs[ sel], R_PMT[ sel], phi[ sel], zd[ sel], theta[ sel]):Ttop_indirect   .Fill(*vals)
        # side
        sel = sel_indirect & np.isin(ihPMT, tubeid_side)
        for vals in zip(zs[ sel], Rs[ sel], z_PMT[ sel], phi[ sel], zd[ sel], theta[ sel]):Tside_indirect  .Fill(*vals)

        # Fill direct tables
        # bottom
        sel = ~sel_indirect & np.isin(ihPMT, tubeid_bottom)
        for vals in zip(zs[ sel], Rs[ sel], R_PMT[ sel], phi[ sel]):Tbottom_direct.Fill(*vals)
        # top
        sel = ~sel_indirect & np.isin(ihPMT, tubeid_top)
        for vals in zip(zs[ sel], Rs[ sel], R_PMT[ sel], phi[ sel]):Ttop_direct   .Fill(*vals)
        # side
        sel = ~sel_indirect & np.isin(ihPMT, tubeid_side)
        for vals in zip(zs[ sel], Rs[ sel], z_PMT[ sel], phi[ sel]):Tside_direct  .Fill(*vals)

    # save photon tables
    if args.verbose: print("Writting tables...")

    fout = ROOT.TFile("scattables.root", "RECREATE")
    fout.WriteObject(Tbottom_indirect, "bottom_indirect")
    fout.WriteObject(Tbottom_direct  , "bottom_direct")
    fout.WriteObject(Ttop_indirect   , "top_indirect")
    fout.WriteObject(Ttop_direct     , "top_direct")
    fout.WriteObject(Tside_indirect  , "side_indirect")
    fout.WriteObject(Tside_direct    , "side_direct")

    Tbottom_indirect.DivideUnnormalized4D(Tbottom_direct)
    Ttop_indirect   .DivideUnnormalized4D(Ttop_direct)
    Tside_indirect  .DivideUnnormalized4D(Tside_direct)

    fout.WriteObject(Tbottom_indirect, "botscattable")
    fout.WriteObject(Ttop_indirect   , "topscattable")
    fout.WriteObject(Tside_indirect  , "sidescattable")
    
    return


if __name__ == "__main__":
    main()