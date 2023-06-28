//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 12 10:33:06 2011 by ROOT version 5.28/00
// from TChain sttree/
//////////////////////////////////////////////////////////

#ifndef scatTableLooper_h
#define scatTableLooper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TScatTable.h"
#include <iostream>
#include "TSystem.h"

#include "TVector3.h"
#include "TRotation.h"

class scatTableLooper {
public :
   TChain*         fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   Int_t nDimensions;
   TScatTable* top_table;
   TScatTable* bot_table;
   TScatTable* side_table;

   TScatTable* top_isodir;
   TScatTable* bot_isodir;
   TScatTable* side_isodir;

   TFile* outfile;

   // Declaration of leaf types
   Float_t         scatvars_zs;
   Float_t         scatvars_rs;
   Float_t         scatvars_ast;
   Float_t         scatvars_ct;
   Float_t         scatvars_phi;
   Float_t         scatvars_zt;
   Float_t         scatvars_rt;
   Float_t         scatvars_ctcor;
   Float_t         scatvars_phicor;

   Float_t         myvars_coseta;
   Float_t         myvars_stdist;

   Float_t srcpos[3];
   Float_t srcdir[3];
   Float_t         oppos[3];
   Float_t         tubepos[3];
   Int_t ihPMT;
   Int_t isct;
   Float_t         opdir[3];

   static const Int_t ntuberows = 1;
   //   double tubezpos[ntuberows];
   TH1D* hphotdist[ntuberows];

//   static const int ntubecolumns = 150;
//   double tubephipos[ntubecolumns];

   static const Int_t nsamplingarcs = 13;
   Double_t arcdist[nsamplingarcs];
   Double_t arcwidth;
   TH1D* hincidentangle[nsamplingarcs];

   Double_t cylRadius;
   Double_t tuberadius;
   Double_t tubezpos;

   // List of branches
   TBranch        *b_scatvars;   //!

   scatTableLooper();
   virtual ~scatTableLooper();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Run(TString outfilenameprefix, bool isNuPRISM, bool ismPMT);
   virtual void     Init(TString outfilenameprefix, int sta, bool isNuPRISM, bool ismPMT);
   virtual void     Loop(bool ismPMT);
   virtual void     Show(Long64_t entry = -1);
   virtual void     MakeTransHistos();
   virtual void     Write();

   // Coordinate transformations:
   TVector3 TfmDetToGlblCoord(TVector3 tmpVct, bool flgPos, bool flgInverse=false);
   TVector3 TransDetToGlbl;
   TRotation RotDetToGlbl;
};

#endif
