#include "TScatTable.h"
#include <iostream>
#include "TROOT.h"
#include "TFile.h"

void makeScatTable() {

  gROOT->SetMacroPath("/pbs/home/g/gdiazlop/Software/WCSimFQTuningTools/STable/old_tools/");
  gROOT->ProcessLine(".L TScatTable_cc.so");

  bool printStuff = true;
  TString name,title;

  TString workdir="/pbs/home/g/gdiazlop/Software/WCSimFQTuningTools/STable/old_tools/";
  TString scatfile = "scattbl_out.root";

  const int nfiles = 2;
  const int ntypes = 3;
  TScatTable *nphot[nfiles][ntypes],*oscttbl[ntypes];
  TString typenames[ntypes] = {"top","bot","side"};
  TString typetitles[ntypes] = {"Top","Bottom","Side"};

  if (printStuff) {
    cout << "opening files" << std::endl;
  }
  TFile* tf = new TFile(workdir+"/"+scatfile);
  if (tf->IsZombie()) {
    cout << "Input file does not exist!";
    exit(-1);
  }
  if (printStuff) cout << "opening output files" << std::endl;
  TString outfilename = workdir + "/fiTQun_scattables.root";
  TFile* outfile = new TFile(outfilename,"RECREATE");
  
    
  for (int itype=0; itype<ntypes; itype++) {
    name = typenames[itype];
    name += "nphot";
    tf->GetObject(name,nphot[0][itype]);
    name = typenames[itype];
    name += "isodir";
    tf->GetObject(name,nphot[1][itype]);
    
    oscttbl[itype] = new TScatTable(nphot[0][itype]);
    name = typenames[itype];
    name += "scattable";
    oscttbl[itype]->SetName(name);
    title = typetitles[itype];
    title += " PMTs Scattering Table";
    oscttbl[itype]->SetTitle(title);
    TScatTable* nphot4d = nphot[1][itype];
    if (printStuff) cout << "Dividing 4D in calculation of " << title << std::endl;
    oscttbl[itype]->DivideUnnormalized4D(nphot4d);
    oscttbl[itype]->PrintMinElement();
    if (printStuff) cout << "writing..." << std::endl;
    oscttbl[itype]->Write();
    
    delete nphot[0][itype];
    delete nphot[1][itype];
  }
//  if (printStuff) cout << "writing..." << std::endl;
//  outfile->Write();
  outfile->Close();

}
