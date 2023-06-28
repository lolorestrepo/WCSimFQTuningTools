#define scatTableLooper_cxx
#include "scatTableLooper.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TSystem.h"
#include "TChain.h"
#include "TChainElement.h"
#include <cstdio>
#include "/pbs/home/g/gdiazlop/Software/WCSim/install/include/WCSimRootGeom.hh"

scatTableLooper::scatTableLooper()
{
  outfile = NULL;
  fChain = NULL;
  top_table = NULL;
  bot_table = NULL;
  side_table = NULL;
}

scatTableLooper::~scatTableLooper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void scatTableLooper::Run(TString outfilenamestring, bool isNuPRISM, bool ismPMT){

  nDimensions = 6;
  Init(outfilenamestring,0,isNuPRISM,ismPMT);
  Loop(ismPMT);

  Write();
  top_table->PrintMinElement();
  bot_table->PrintMinElement();
  side_table->PrintMinElement();

}

Int_t scatTableLooper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t scatTableLooper::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
     fCurrent = chain->GetTreeNumber();
     cout << "Processing Tree #" << fCurrent << ", ";
     cout << "# of entries : " << (fChain->GetTree())->GetEntries() << endl;
   }
   return centry;
}


void scatTableLooper::Init(TString outfilenamestring, int sta, bool isNuPRISM, bool ismPMT)

{

  // Set branch addresses and branch pointers
  if (fChain) {
    delete fChain;
    fChain = NULL;
  }

  fChain = new TChain("sttree","");
  
  fChain->Add(outfilenamestring);
  
  if (!(fChain->GetNtrees()>0)) {
    cout << "Trees not loaded!" << endl;
    exit(-1);
  }
  
  //  std::cout << "Added " << fChain->GetNtrees() << " files to TChain" << std::endl;

  cout << "deleting outfile" << endl;
  if (outfile) {
    delete outfile;
    outfile = NULL;
  }

  gStyle->SetOptFit(111111);

  //  TString outfilename = outfilenameprefix;
  //  outfilename += "/scattbl_out.root";
  TString outfilename = "scattbl_out.root";
  cout << "Making outfile " << outfilename << endl;
  outfile = new TFile(outfilename,"RECREATE");

  // cout << "deleting outfile" << endl;
  // if (outfile) {
  //   delete outfile;
  //   outfile = NULL;
  // }

  gStyle->SetOptFit(111111);

// //  TString outfilename = outfilenameprefix;
// //  outfilename += "/scattbl_out.root";
//   TString outfilename = "scattbl_out.root";
//   cout << "Making outfile " << outfilename << endl;
//   outfile = new TFile(outfilename,"RECREATE");

  if (top_table) {
    delete top_table;
    top_table = NULL;
  }
  if (bot_table) {
    delete bot_table;
    bot_table = NULL;
  }
  if (side_table) {
    delete side_table;
    side_table = NULL;
  }

  // Extract geometry information from WCSim file
  TChain * tempGeoChain = new TChain("wcsimGeoT","tempGeoChain");
  tempGeoChain->Add(outfilenamestring, 1);
  
  WCSimRootGeom*  geoWCSim = new WCSimRootGeom();
  TBranch * bWCSimRootGeom = tempGeoChain->GetBranch("wcsimrootgeom");
  bWCSimRootGeom->SetAddress(&geoWCSim);
  tempGeoChain->GetBranch("wcsimrootgeom")->SetAutoDelete(kTRUE);
  tempGeoChain->GetEntry(0);
  
  //  tuberadius = 10.16;
  tuberadius = geoWCSim->GetWCPMTRadius();
  if(ismPMT) tuberadius = 20.; // OK -- only used to define table dimensions. CV
  //  tubezpos = 500.; 
  tubezpos = geoWCSim->GetWCCylLength()/2.;
  // 5010.545
  //  cylRadius = 300.;
  cylRadius = geoWCSim->GetWCCylRadius();
  if(isNuPRISM){
    if(ismPMT){
      tubezpos  = 136.95;
      cylRadius = 172.05;
    }else{
      tubezpos = 500.;
      cylRadius = 400.;
    }
  }

  std::cout<<"PMT radius: "<<tuberadius<<"cm, Detector half height: "<<tubezpos
            << "cm, Detector radius: "<<cylRadius<<"cm"<<std::endl;

  double rmax = cylRadius - tuberadius;
  double zmin = - (tubezpos - tuberadius);
  double zmax = tubezpos - tuberadius;

  
  // Set up detector transformation
  TVector3 DetZinGlbl(0.,0.,1.);// Detector local Z-axis in global coordinate system
  TVector3 DetXinGlbl(1.,0.,0.);// Detector local X-axis in global coordinate system
  
  float WCXRotation[3] = {1.,0.,0.};
  float WCYRotation[3] = {0.,1.,0.};
  float WCZRotation[3] = {0.,0.,1.};
  
  TChain * settingsChain = new TChain("Settings","settings");
  settingsChain->Add(outfilenamestring, 1);

  if (!settingsChain){

    std::cout<<"There is no Settings TTree in "<<outfilenamestring<<std::endl;
  } else {
    TBranch *bWCXRotation  = settingsChain->GetBranch("WCXRotation");
    if (bWCXRotation) {
      bWCXRotation->SetAddress(WCXRotation);
    } else {
      WCXRotation[0] = 1.;
      WCXRotation[1] = 0;
      WCXRotation[2] = 0;
    }
    TBranch *bWCYRotation  = settingsChain->GetBranch("WCYRotation");
    if (bWCYRotation) {
      bWCYRotation->SetAddress(WCYRotation);
    } else {
      WCYRotation[0] = 0;
      WCYRotation[1] = 1.;
      WCYRotation[2] = 0;
    }
    TBranch *bWCZRotation  = settingsChain->GetBranch("WCZRotation");
    if (bWCZRotation) {
      bWCZRotation->SetAddress(WCZRotation);
    } else {
      WCZRotation[0] = 0;
      WCZRotation[1] = 0;
      WCZRotation[2] = 1.;
    }
    settingsChain->GetEntry(0);
  }

  TransDetToGlbl.SetXYZ(geoWCSim->GetWCOffset(0),geoWCSim->GetWCOffset(1),geoWCSim->GetWCOffset(2));
  
  DetZinGlbl.SetXYZ(WCZRotation[0],WCZRotation[1],WCZRotation[2]);
  DetXinGlbl.SetXYZ(WCXRotation[0],WCXRotation[1],WCXRotation[2]);
  RotDetToGlbl.SetZAxis(DetZinGlbl,DetXinGlbl);// Rotation matrix which transforms detector local frame to global frame
  
  
  int nrbinss = 16;
  int nrbinst = 16;
  if(ismPMT) nrbinst = 8;
  int nzbinss = 35;
  int nzbinst = 35;
  if(ismPMT) nzbinst = 16;

  int nangbins = 16;
  int nctbins = 16;

//  tubezpos = 5010.545*.1;
//  tuberadius = 10.16;
//  int nrbins = 16;
//  double rmax = 480.;
//  int nzbins = 18;
//  double zmin = -280.;
//  double zmax = 280.;
//  int nangbins = 16;
//  int nctbins = 16;

  // tubezpos = 1818.921111111;
  // tuberadius = 25.4;
  // int nrbins = 16;
  // double rmax = 1600.;
  // int nzbins = 35;
  // double zmin = -1750.;
  // double zmax = 1750.;
  // int nangbins = 16;
  //  int nctbins = 16;

  // Hyper=K Large Cylinders
  

  
  //  tuberadius = 25.4;


  double pinumber = TMath::Pi();
  double oneplusepsilon = 1.00001;

  arcwidth = 50.;
  for (int isa=0; isa<nsamplingarcs; isa++) {
    arcdist[isa] = 500. + isa*50.;
  }

//  double deltaphi = 2.*TMath::Pi()/ntubecolumns;
//  for (int iphi=0; iphi<ntubecolumns; iphi++) {
//    tubephipos[iphi] = (-TMath::Pi()+deltaphi/2.) + iphi*deltaphi;
//    //    std::cout << "phipos = " << tubephipos[iphi] << std::endl;
//  }

  for (int ipd=0; ipd<ntuberows; ipd++) {
    TString hname = "hphotdist";
    hname += ipd;
    TString htitle = "Distance between source position and tube position";
    hphotdist[ipd] = new TH1D(hname,htitle,1000,0.,3500.);
    hphotdist[ipd]->GetXaxis()->SetTitle("Distance to PMT (cm)");
  }

  for (int isa=0; isa<nsamplingarcs; isa++) {
    TString hname = "hincidentangle";
    hname += isa;
    TString htitle = Form("Cosine of photon incident angle for a source/PMT distance between %8.1f and %8.1f cm",arcdist[isa],arcdist[isa]+arcwidth);
    hincidentangle[isa] = new TH1D(hname,htitle,100,-1.,1.);
    hincidentangle[isa]->GetXaxis()->SetTitle("Cos(Incident angle)");
  }

  top_table = new TScatTable("topnphot","Number of Photons for Tubes on the Top of the Tank",
           nzbinss, zmin, zmax, nrbinss, 0., rmax, nrbinst, 0., rmax,
           nangbins, -pinumber*oneplusepsilon, pinumber*oneplusepsilon,
           nctbins, -oneplusepsilon, oneplusepsilon,
           nangbins, -pinumber*oneplusepsilon, pinumber*oneplusepsilon);
  bot_table = new TScatTable("botnphot","Number of Photons for Tubes on the Bottom of the Tank",
           nzbinss, zmin, zmax, nrbinss, 0., rmax, nrbinst, 0., rmax,
           nangbins, -pinumber*oneplusepsilon, pinumber*oneplusepsilon,
           nctbins, -oneplusepsilon, oneplusepsilon,
           nangbins, -pinumber*oneplusepsilon, pinumber*oneplusepsilon);
  side_table = new TScatTable("sidenphot","Number of Photons for Tubes on the Side of the Tank",
           nzbinss, zmin, zmax, nrbinss, 0., rmax, nzbinst, zmin, zmax,
           nangbins, -pinumber*oneplusepsilon, pinumber*oneplusepsilon,
           nctbins, -oneplusepsilon, oneplusepsilon,
           nangbins, -pinumber*oneplusepsilon, pinumber*oneplusepsilon);
  
  top_isodir = new TScatTable("topisodir","Number of Photons for Tubes on the Top of the Tank",
            nzbinss, zmin, zmax, nrbinss, 0., rmax, nrbinst, 0., rmax,
           nangbins, -pinumber*oneplusepsilon, pinumber*oneplusepsilon);
  bot_isodir = new TScatTable("botisodir","Number of Photons for Tubes on the Bottom of the Tank",
           nzbinss, zmin, zmax, nrbinss, 0., rmax, nrbinst, 0., rmax,
           nangbins, -pinumber*oneplusepsilon, pinumber*oneplusepsilon);
  side_isodir = new TScatTable("sideisodir","Number of Photons for Tubes on the Side of the Tank",
           nzbinss, zmin, zmax, nrbinss, 0., rmax, nzbinst, zmin, zmax,
           nangbins, -pinumber*oneplusepsilon, pinumber*oneplusepsilon);

  fCurrent = -1;
  cout << "Setting Make Class" << endl;
  fChain->SetMakeClass(1);
  
  //  fChain->SetBranchStatus("srcpos",0);
  //  fChain->SetBranchStatus("srcdir",0);
  //  fChain->SetBranchStatus("ihPMT",0);

  cout << "Setting Branch Address" << endl;
  //  fChain->SetBranchAddress("scatvars", &scatvars_zs, &b_scatvars);
  fChain->SetBranchAddress("srcpos",srcpos);
  fChain->SetBranchAddress("srcdir",srcdir);
  fChain->SetBranchAddress("oppos", oppos);
  fChain->SetBranchAddress("tubepos", tubepos);
  fChain->SetBranchAddress("ihPMT",&ihPMT);
  fChain->SetBranchAddress("isct",&isct);
  fChain->SetBranchAddress("opdir", opdir);

}

void scatTableLooper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

void scatTableLooper::Loop(bool ismPMT)
{
  
  if (fChain == 0) return;

//  Long64_t nentries = fChain->GetEntries();//this is slow!
//  cout << "chain has " << nentries << " entries" << endl;

  Long64_t jentry=0;
  while (1) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) {//no more entries!
      cout << "Total " << jentry << " entries" << endl;
      break;
    }
    
    fChain->GetEntry(jentry);

    Float_t capthr = 0.99;
    if(ismPMT) capthr = 0.98;
    

    TVector3 sptmp(srcpos);
    TVector3 sdtmp(srcdir);
    TVector3 tptmp(tubepos);

    TVector3 sp = TfmDetToGlblCoord(sptmp, true, true);
    TVector3 sd = TfmDetToGlblCoord(sdtmp, false, true);
    TVector3 tp = TfmDetToGlblCoord(tptmp, true, true);

    
    sp *= 0.1;
    tp *= 0.1;

    scatvars_zs = sp.Z();
    scatvars_rs = sp.Perp();
    scatvars_ast = tp.DeltaPhi(sp);
    scatvars_ct = sd.Z();
    scatvars_phi = sd.DeltaPhi(tp-sp);
    scatvars_zt = tp.Z();
    scatvars_rt = tp.Perp();


    TVector3 td;
    if (scatvars_zt > tubezpos*capthr) {
      td.SetXYZ(0.,0.,-1.);
    } else if (scatvars_zt < -tubezpos*capthr) {

      td.SetXYZ(0.,0.,1.);
    } else {
      td = -tp;
      td.SetZ(0.);
      td = td.Unit();
    }
    myvars_coseta = -td.Dot((tp-sp).Unit());
    myvars_stdist = (tp-sp).Mag();
//    if ((fabs(scatvars_zt)<600.)&&(myvars_coseta > 0.99)) {
//      TVector3 udiff = (tp-sp).Unit();
//      std::cout << "WTF? coseta = " << myvars_coseta << ", stdist = " << myvars_stdist << ", td, tp, sp, (tp-sp).Unit() = (" << td.X() << "," << td.Y() << "," << td.Z() << "), (" << tp.X() << "," << tp.Y() << "," << tp.Z() << "), (" << sp.X() << "," << sp.Y() << "," << sp.Z() << "), (" << udiff.X() << "," << udiff.Y() << "," << udiff.Z() << ")" << std::endl;
//    }

    //    std::cout << "outside, scatvars_zs = " << scatvars_zs << std::endl;


    if (isct != 0) {//indirect photon
      if (scatvars_zt > 110.) {
        top_table->Fill(scatvars_zs,scatvars_rs,scatvars_rt,scatvars_ast,scatvars_ct,scatvars_phi);
      } else if (scatvars_zt < -110.) {

        bot_table->Fill(scatvars_zs,scatvars_rs,scatvars_rt,scatvars_ast,scatvars_ct,scatvars_phi);
      } else {
        side_table->Fill(scatvars_zs,scatvars_rs,scatvars_zt,scatvars_ast,scatvars_ct,scatvars_phi);
      }
      //      cout << "Fills direct " << jentry << endl;
    } else {//direct photon
      if (scatvars_zt > 110.) {
        top_isodir->Fill(scatvars_zs,scatvars_rs,scatvars_rt,scatvars_ast);
      } else if (scatvars_zt < -110.) {

        bot_isodir->Fill(scatvars_zs,scatvars_rs,scatvars_rt,scatvars_ast);
      } else {
        side_isodir->Fill(scatvars_zs,scatvars_rs,scatvars_zt,scatvars_ast);
	//	std::cout << "ihPMT = " << ihPMT << std::endl;
	//	if ((ihPMT>=1888)&&(ihPMT<=1938)) {
	MakeTransHistos();

	  //}
      }
      //      cout << "Fills indirect "<< jentry << endl;
    }
    jentry++;
  }
}

void scatTableLooper::MakeTransHistos() {

//  tubezpos[1] = 1696.79993;
//  tubezpos[0] = 1763.79993;
//  std::cout << "inside, scatvars_zs = " << scatvars_zs << std::endl;

  if (fabs(sin(scatvars_ast)*scatvars_rs) < tuberadius) {  // Photons in front of the PMT (in xy plane)
    if (fabs(scatvars_zs-scatvars_zt) < tuberadius) {  // Photons in front of the PMT (in z axis)
      for (int izpos=0; izpos<ntuberows; izpos++) {
	//	if (scatvars_zt < tubezpos[izpos] + tuberadius) {
	//	  double tsq = scatvars_zt*scatvars_zt+1690.*1690.;
	  double tsq = scatvars_zt*scatvars_zt+cylRadius*cylRadius;
	  double ssq = scatvars_zs*scatvars_zs+scatvars_rs*scatvars_rs;
	  //      	  double tdots = scatvars_zt*scatvars_zs+1690.*scatvars_rs*cos(scatvars_ast);

	  double tdots = scatvars_zt*scatvars_zs+cylRadius*scatvars_rs*cos(scatvars_ast);
	  double rsq = tsq+ssq-2.*tdots;
	  //	  hphotdist[izpos]->Fill(sqrt(rsq),rsq);
	  hphotdist[0]->Fill(sqrt(rsq));

	  //	  hphotdist[izpos]->Fill(1690.-cos(scatvars_ast)*scatvars_rs);
	  //}
      }
    }
    //    if (fabs(scatvars_zt)<600.) {
    if (fabs(scatvars_zt) < tubezpos) {

      if ((myvars_stdist > arcdist[0]) && (myvars_stdist < arcdist[nsamplingarcs-1])) {
	for (int isa=0; isa<nsamplingarcs; isa++) {
	  if ((myvars_stdist > arcdist[isa]) && (myvars_stdist < arcdist[isa]+arcwidth)) {
	    hincidentangle[isa]->Fill(myvars_coseta);
	    break;
	  }
	}
      }
    }
  }

}

void scatTableLooper::Write() {

  cout << "writing outfile" << endl;
  if (outfile) {
    outfile->cd();
    top_table->Write();
    bot_table->Write();
    side_table->Write();
    top_isodir->Write();
    bot_isodir->Write();
    side_isodir->Write();
//    outfile->Write();
    for (int ipd=0; ipd<ntuberows; ipd++) {
      TString funcname = "myfunc";
      funcname += ipd;
      TF1* tf1 = new TF1(funcname,"[0]*exp(-x/[1])/(x*x+645.16)",100.,3500.);
      tf1->SetParameter(0,1.e8);
      tf1->SetParameter(1,7000.);
      hphotdist[ipd]->Fit(tf1,"","",500.,3000.);
      hphotdist[ipd]->Write();
    }
    for (int isa=0; isa<nsamplingarcs; isa++) {
      hincidentangle[isa]->Write();
    }
    outfile->Close();
    delete outfile;
    outfile = NULL;
  } else {
    cout << "Why is there no outfile?" << endl;
  }
  cout << "done writing outfile" << endl;

}

TVector3 scatTableLooper::TfmDetToGlblCoord(TVector3 tmpVct, bool flgPos, bool flgInverse){
  
  if (flgInverse) {
    if (flgPos) tmpVct-=TransDetToGlbl;
    tmpVct*=RotDetToGlbl.Inverse();
  }
  else {
    tmpVct*=RotDetToGlbl;
    if (flgPos) tmpVct+=TransDetToGlbl;
  }
  
  return tmpVct;
  
}
