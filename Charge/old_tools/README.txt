
makeChargePDFplot.C:

  gSystem->Load("/pbs/home/g/gdiazlop/Software/WCSim/install/lib/libWCSimRoot.so");
  .L makeChargePDFplot.C++
  makeChargePDFplot("../../Simulation/local/out/out.root", "ChargePDF.root")


gen2d.cc:
  .L gen2d.cc++
  gen2d()

fitpdf.cc:
  .L /pbs/home/g/gdiazlop/Software/fiTQun/fiTQun/fQChrgPDF.cc++
  .L /pbs/home/g/gdiazlop/Software/fiTQun/fiTQun/fQChrgPDF_cc.so
  .L fitpdf.cc++







SKDETSIM:
Charge pdf generation mode: setenv MEAN_PE [mu value], then chrgpdf_mode=1
Set darkrate to 0 in card file
absorb all optical photons in SGPMT
in gutrev, call genpoispe after gtreve is called!
in dsres_multi: pltchdst fills histograms for unhit prob calculation!



How to:

Just run:
./RunAll_detsimjobs.csh [SK Ver]
  Reads in list of mu values, at each mu, creates [mu].dat(zbs)

Go to the output directory.
./Run_filljob.csh [SK Ver]
  Reads in list of mu values, and genchrgpdf at each mu, creates [mu]_pdf.root(contains hchpdf? which is the 1d Q_obs distribution of old/new PMT at given mu value)
  gen2d.cc reads in list of mu values, read 1D charge distributions from [mu]_pdf.root and make it into 2d histo hst2d_type?, obtain Punhit coefficients, write them all in pdf2d.root
  plotChrgPDF.cc produces graphs of mu lower/upper threshold from 2d pdf histogram and outputs in muthresh.root
  fitpdf.cc fits the pdf at different q range for each PMT type
  MakecPDFparFile.cc combines the parameters for all PMT types, and produces the final root file, cPDFpar.root.


================================================================================
  MAKING WCSIM CHARGE PDFs
================================================================================

 - Environment
   - Get a version of WCSim that has the pmtPoisson option. Currently (30 Apr 2015) you can get it by running:
     - git clone -b feature/fiTQunTuning git@github.com:cvilelasbu/WCSim.git
   - Set FITQUN_ROOT to point at your installation of fiTQun
   - Setup ROOT, Geant4, WCSim (e.g., by sourcing the nuPRISM software setup script)

 - To generate the MC:
   - Look at the script SukapSampleScripts/runWCSimChargePDFs.sh and "base" WCSim macro file SukapSampleScripts/chargePDF_base.mac

 - To make the charge PDF histograms run the makeChargePDFplot.C macro on all files generated
   - Script: SukapSampleScripts/makeChargePDFPlots.sh

 - To fit the histograms and generate a cPDFpar file that can be used with fiTQun, link all the histogram root files resulting from the previous step to the chrgpdf directory and then run the following sequence of macros: gen2d.cc, plotChrgPDF.cc, fitpdf.cc'(0,4)', fitpdf.cc'(1,4)' and MakecPDFparFile.cc
   - Note: the '4' in the fitpdf.cc is for SK4. This determines the threshold to be used and SK4 is the highest option available. Edit fitpdf.cc and change this number to something high.... the threshold is not simulated by WCSim
   - Sample script: SukapSampleScripts/makeChargePDFFile.sh
   - This script will also run the FitMu diagnostics tool if it has been previously compiled (make)
