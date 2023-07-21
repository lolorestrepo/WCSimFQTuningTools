# Introduction

FiTQun requires an accurate description of the photosensor charge response to accurately reconstruct events. This involves obtaining a probability density function for the "observed" (i.e., digitized) charge, given a prediction of μ photoelectrons having been produced at the photocathode. The procedure is repeated for a range of μ, and the resulting distributions are fitted with polynomial functions (with Gaussian tails). In addition, the probability of registering a hit as a function of μ is plotted and fitted with a modified Poisson expectation. The parameters resulting from the fit constitute the tune and are given to fiTQun at runtime.

# Tuning procedure
These are the instructions to produce a charge PDF tune for a WCSim photosensor.

The following packages are required:
 * [WCSim](https://github.com/nuPRISM/WCSim) with [WCSimFQTuner](https://github.com/fiTQun/WCSimFQTuner):
 * [fiTQun](https://github.com/fiTQun/fiTQun)

## Setting up
The shell environment should be set up such that the following variables `ROOTSYS`, `WCSIMDIR` and `FITQUN_ROOT` are appropriately set.

### Compiling
If you want to use the diagnostics tool (highly recommended, more on that later), `chargePDFdiagnostics` must be compiled. Runing `make` in the `Utilities/chrgpdf` should do the trick.

## Producing the Monte Carlo
The first step in the procedure is to produce a dedicated MC sample that contains no physical particles, only photoelectrons added directly to the raw hit collection, according to Poisson distributions with means **μ**. A list containing the **μ** that should be generated is given in `Utilities/chrgpdf/mutbl.txt`.

To produce the MC, add the lines: `/mygen/pmtPoisson true` and `/mygen/poissonMean 3.0`, for μ = 3.0. In addition it is recommended<sup>1</sup> to turn off dark noise and primary particle generation (or set to something below Cherenkov threshold). The number of events necessary depends heavily on the value of μ and on the number of photosensors in the chosen geometry. For a Super-K like tank (~11k tubes), 80 events are enough for μ up to 30, and that number can be reduced to 20 for μ above 30.

A sample WCSim macro file to generate the point μ = 3.0 for the `20inchBandL` in the `Cylinder_60x74_20inchBandL_14perCent` detector configuration can be found in `Utilities/chrgpdf/chargePDF_SAMPLE.mac`

A bash script that generates the full set of samples is given, for reference, in:
`Utilities/chrgpdf/SukapSampleScripts/runWCSimChargePDFs.sh`

## Making the observed charge PDF plots
Once the MC has been generated, the standard WCSim output is used to generate the charge PDF plots. The ROOT macro `Utilities/chrgpdf/makeChargePDFplot.C` will produce the charge PDF plot for a given WCSim file. It can be used, for example, by running: `root -b -q 
Utilities/chrgpdf/makeChargePDFplot.C'(inFile.root, outFile.root)'`. Where `inFile.root` is the WCSim output and `outFile.root` is a root file where the histogram will be stored. To make things easier from here onwards it is recommended that `outFile.root` is named `3.0_pdf.root` for μ = 3.0, or `100_pdf.root` for μ = 100, etc and that they are placed in directories `Mu_<μ>`.

Note that these histograms have a variable bin width, and also that the x-axis is the logarithm of the observed charge in p.e. An example of the raw output of `makeChargePDFplot.C` for μ = 5.0 is shown below.
![](http://nngroup.physics.sunysb.edu/~cvilela/FilesForGitHubDoc/sampleRawCPDFhist.png)

An example of a bash script that runs `makeChargePDFplot.C` over a complete MC set is given in:
`Utilities/chrgpdf/SukapSampleScripts/makeChargePDFPlots.sh`

## Running the fits
Once the PDF plots have been created for all values of μ, a series of scripts need to be run which will combine the PDF distributions into a 2D histogram, fit the curves with the prescribed functional forms, and store the fitted parameters in a cPDFpar.root file that can be interpreted by fiTQun.

All these scripts assume that the files containing the PDFs are in the `Utilities/chrgpdf` directory and that they are named `<μ>_pdf.root`. The first step is to create symbolic links in the `Utilities/chrgpdf` directory pointing to the relevant files. For example:
```bash
while read mu
do
    ln -s /dir/where/I/put/my/mc/Mu_${mu}/${mu}_pdf.root Utilities/chrgpdf/
done < mutbl.txt
```

With the links in `Utilities/chrgpdf`, the following root scripts should be run, in order:
```
root -b -q gen2d.cc
root -b -q plotChrgPDF.cc
root -b -q fitpdf.cc'(0,4)'
root -b -q fitpdf.cc'(1,4)'
root -b -q MakecPDFparFile.cc
```
The first script produces a 2D histogram containing all the charge PDFs generated. Here the hit probabilities are also fitted with a modified Poisson curve. One error that has occurred in the past is that the modified Poisson fit gives unphysical results. In that case `gen2d.cc` will quit with the message `ERROR: Bad Phit fit.` In that case, get in touch with @cvilelasbu to discuss how you can get around this.

The second script will identify the boundaries of the charge PDF distributions so that the points where the fit goes from polynomial to Gaussian can be defined. At this stage, a `plt_pdf2d_hst2d_type0_logz.pdf` file will be produced which will display the 2D plot of the charge PDFs, along with the boundary. The boundary sometimes comes out a bit buggy, but we don't think this is a problem. An example is shown below.
![](http://nngroup.physics.sunysb.edu/~cvilela/FilesForGitHubDoc/plt_pdf2d_hst2d_type0_logz.png)

The `fitpdf.cc` script will perform the actual fits to the charge PDFs. The first argument is the PMT "type". This is an SKDETSIM relic, where two types of PMT exists: old and new. In the future we will be generalizing this to work in possible WCSim multi-PMT configurations. For now we run on the same files for "old" and "new" PMTs so that all histograms get filled. The second argument is the SK era. We use `4` as the charge threshold in the early eras is low compared to anything we might expect for a future detector.

Finally, `MakecPDFparFile.cc` will compile the results of the fit into the file `cPDFpar.root`. This is the file that fiTQun needs.

An example of a script which performs all the steps above (and the check below) is given in:
`Utilities/chrgpdf/SukapSampleScripts/makeChargePDFFile.sh`

## Checking the results
In order to verify that the procedure above is giving the intended results, a test was set up where fiTQun is used to fit the MC events generated in the first step of this procedure with a uniform charge distribution over all photosensors. The fit should return a value close to the true μ.

The `chargePDFdiagnostics` interfaces with fiTQun libraries to produce plots of μ<sub>reco</sub>/μ<sub>true</sub>. For the tool to work out of the box, the WCSim output should be in the following directory structure: `SOME_BASE_DIR/Mu_<μ>/<μ>_pdf.root`

To produce the verification plot, the default `cPDFpar.root` in fiTQun should be replaced by the new one. Since we are not using a "parameter override" file, the fiTQun default is `cPDFpar_sk1.root`.

Replace it:
```
mv ${FITQUN_ROOT}/const/cPDFpar_sk1.root ${FITQUN_ROOT}/const/cPDFpar_sk1.root.backup
cp cPDFpar.root ${FITQUN_ROOT}/const/cPDFpar_sk1.root
```

and run `chargePDFdiagnostics`:

`Utilities/chrgpdf/chargePDFdiagnostics Utilities/chrgpdf/SukapSampleScripts/mutbl_reduced.txt SOME_BASE_DIR FitMuOut.root`

Where `mutbl_reduced.txt` is an abbreviated version of `mutbl.txt`. It takes a long time to run `chargePDFdiagnostics` on the entire set.

The output plot should look something like this:
![](http://nngroup.physics.sunysb.edu/~cvilela/FilesForGitHubDoc/FitMuOut.png)

The systematic pulls of < ± 5% at low charge are known, and consistent with what is obtained for Super-K PMT tunes using SKDETSIM. The last two points in the plot are above the threshold used to produce this tune (SK4 threshold is ~1000 p.e.)

Don't forget to move `cPDFpar_sk1.root` back once you're done.
```
mv ${FITQUN_ROOT}/const/cPDFpar_sk1.root.backup ${FITQUN_ROOT}/const/cPDFpar_sk1.root
```

## Adding tune to fiTQun
To use the new charge PDF tune in fiTQun, rename it `cPDFpar_<WCSimPMTType>.root`, where `<WCSimPMTType>` should be the PMT name in WCSim, for example, `20inchBandL`. Then specify the PMTType in the parameter override file for the WCSim configuration you want to run on. See `${FITQUN_ROOT}/ParameterOverrideFiles` for examples. To include your tune in the official fiTQun releases, get in touch with @cvilelasbu.

***
<sup>1</sup> For the moment, the `/mygen/pmtPoisson` flag should take care of all this, but that might not be the case in the future.
