import sys
import argparse
import uproot
import ROOT
import numpy  as np
import matplotlib.pyplot as plt


def main():

    ############ Program arguments ############
    parser = argparse.ArgumentParser( prog        = "fit angular response"
                                    , description = "description"
                                    , epilog      = """""")
    
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-i", "--infile", type=str, nargs="?", help = "", default="angular.root")

    parser.add_argument("polydeg", type=int, nargs=1, help = "")

    args = parser.parse_args()
    ##########################################

    polydeg = args.polydeg[0]

    # read angular responses file
    with uproot.open(args.infile) as infile:
        # read angular response histogram computed from all the PMTs (side and caps)
        h, rbins, etabins = infile["AngularResponse_all"].to_numpy()
        # select only histogram for the first bin in radious (closest to the PMT)
        h = h[0]
    
    # compute bin centers and fit angular response
    etas = (etabins[1:] + etabins[:-1])/2.
    pars = np.polyfit(etas, h, polydeg)
    
    # create polynomial TF1 and set the computed parameters
    tf1 = ROOT.TF1("angResp", f"pol{polydeg}", etabins[0], etabins[-1])
    tf1.SetParameters(*np.flip(pars))

    if args.verbose:
        # plot results
        plt.figure(figsize=[5, 3])
        plt.scatter(etas, h, color="k")
        x = np.linspace(etabins[0], etabins[-1], 100)
        fitted = [tf1.Eval(eta) for eta in x]
        plt.plot(x, fitted, color="k")
        plt.xlabel(r"$\cos \eta$")
        plt.ylabel(r"$\epsilon(\cos \eta)$")
        plt.tight_layout()
        plt.show()

    # write output
    fout = ROOT.TFile("fitted_angular.root", "RECREATE")
    fout.WriteObject(tf1, "angResp")
    fout.Close()
    
    return


if __name__ == "__main__":
    main()