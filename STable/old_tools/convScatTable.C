#include <iostream>
#include "TScatTable.h"
#include "TScatTableF.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"

void convScatTable() {
// convert TScatTable to TScatTableF
  
  gROOT->SetMacroPath("/pbs/home/g/gdiazlop/Software/WCSimFQTuningTools/STable/old_tools/");
  gROOT->ProcessLine(".L TScatTable_cc.so");
  gROOT->ProcessLine(".L TScatTableF_cc.so");

  TString workdir="/pbs/home/g/gdiazlop/Software/WCSimFQTuningTools/STable/old_tools/";
  bool printStuff = true;
  TString name,title;
  
  TFile* tf = new TFile(workdir+"/fiTQun_scattables.root");
  if (tf->IsZombie()) {
    cout << "Input file does not exist!";
    exit(-1);
  }
  TFile* outfile = new TFile(workdir+"/fiTQun_scattablesF.root","RECREATE");
  
  const int ntypes = 3;
  TScatTable *tmpscttbl[ntypes];
  TScatTableF *oscttbl[ntypes];
  TString typenames[ntypes] = {"top","bot","side"};
  TString typetitles[ntypes] = {"Top","Bottom","Side"};
    
  for (int itype=0; itype<ntypes; itype++) {
    name = typenames[itype];
    name += "scattable";
    tf->GetObject(name,tmpscttbl[itype]);
    oscttbl[itype] = new TScatTableF(tmpscttbl[itype]);
    oscttbl[itype]->SetName(name);
    if (printStuff) cout << "writing..." << name << std::endl;
    oscttbl[itype]->Write();
  }
//  if (printStuff) cout << "writing..." << std::endl;
//  outfile->Write();
  outfile->Close();

}
