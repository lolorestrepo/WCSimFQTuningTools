#include "TScatTable.h"
#include "TROOT.h"
#include "TFile.h"
#include <iostream>

using namespace std;

void mergetables() {

  gROOT->SetMacroPath("/pbs/home/g/gdiazlop/Software/WCSimFQTuningTools/STable/old_tools");  
  gROOT->ProcessLine(".L TScatTable.cc++");

  const int nfiles = 100;
  const int ntypes = 6;
  TScatTable *stblmrgd[ntypes],*nphot[ntypes];
  TString typenames[ntypes] = {"topnphot","botnphot","sidenphot","topisodir","botisodir","sideisodir"};

  //  TString workdir="/scratch/k/konaka/wilking/production/scat_Cylinder_60x74_20inchBandL_40perCent_e-_ene3_100000evts/scatfiles";
  TString workdir="/pbs/home/g/gdiazlop/Software/WCSimFQTuningTools/STable/old_tools/";

  TFile* outfile = new TFile(workdir+"/scattbl_merged.root","RECREATE");

  for (int ifile=0; ifile<nfiles; ifile++) {
    //TFile* tf = new TFile(workdir+Form("/scat_%d/scattbl_out.root",ifile));
    TFile* tf = new TFile(workdir+"/scattbl_out.root");
    if (tf->IsZombie()) {
      cout << "Input file does not exist!";
      exit(-1);
    }
    cout << "File #" << ifile << endl;
    for (int itype=0; itype<ntypes; itype++) {
      tf->GetObject(typenames[itype],nphot[itype]);

      if (ifile==0) {
        outfile->cd();
        cout << "File" << endl;

        stblmrgd[itype] = new TScatTable(nphot[itype]);
        
        cout << "File" << ifile << endl;

        stblmrgd[itype]->Reset();
        tf->cd();
      }

      cout << "Adding " << typenames[itype] << endl;
      stblmrgd[itype]->Add(nphot[itype]);
      delete nphot[itype];
    }
    delete tf;
  }

  outfile->cd();
  cout << "Writing..." << endl;
  for (int itype=0; itype<ntypes; itype++) {
    stblmrgd[itype]->SetName(typenames[itype]);
    stblmrgd[itype]->Write();
  }
  outfile->Close();

}
