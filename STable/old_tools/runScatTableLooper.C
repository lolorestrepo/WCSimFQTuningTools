#include "scatTableLooper.h"

void runScatTableLooper(TString strpath="./", bool isNuPRISM=1, bool ismPMT=1)
{
  cout << "Working in: " << strpath << endl;
  
  gSystem->Load("/pbs/home/g/gdiazlop/Software/WCSim/install/lib/libWCSimRoot.so");

  gROOT->SetMacroPath("/pbs/home/g/gdiazlop/Software/WCSimFQTuningTools/STable/old_tools");  
  gROOT->ProcessLine(".L TScatTable.cc+");
  gROOT->ProcessLine(".L TScatTableF.cc+");
  gROOT->ProcessLine(".L scatTableLooper.C+");

  scatTableLooper stl;
  stl.Run(strpath, isNuPRISM, ismPMT);//set the path where the *sttree.root files are

}


