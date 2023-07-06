#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TSystem.h"
#include "TVector3.h"
#include "/pbs/home/g/gdiazlop/Software/WCSim/install/include/WCSimRootGeom.hh"

TVector3 TransDetToGlbl;
TRotation RotDetToGlbl;
TVector3 TfmDetToGlblCoord(TVector3 tmpVct, bool flgPos, bool flgInverse);

void angularResponsePlotter_v1 ( std::string inList, std::string outFName, int isNuPRISM = 0, int ismPMT = 0){

  gSystem->Load("/pbs/home/g/gdiazlop/Software/WCSim/install/lib/libWCSimRoot.so");
  std::cout << "WCSim ROOT lib loaded" << std::endl;

  std::ifstream inListFS(inList.c_str());
  std::string line;
  std::vector<std::string> rootFileNames;
  while (std::getline(inListFS, line)){
    rootFileNames.push_back(line);
  }

  if (!rootFileNames.size()){
    std::cout << "No root file names provided. Quitting..." << std::endl;
    return;
  }

  gStyle->SetOptStat(0);
  
  TFile * fWCSim = new TFile (rootFileNames.at(0).c_str());
  TTree * tWCSimGeo = (TTree*) fWCSim->Get("wcsimGeoT");
  WCSimRootGeom*  geoWCSim = new WCSimRootGeom();
  

  TBranch * bWCSimRootGeom = tWCSimGeo->GetBranch("wcsimrootgeom");
  bWCSimRootGeom->SetAddress(&geoWCSim);
  tWCSimGeo->GetBranch("wcsimrootgeom")->SetAutoDelete(kTRUE);

  tWCSimGeo->GetEntry(0);

  double detRadius = geoWCSim->GetWCCylRadius();
  double detLength = geoWCSim->GetWCCylLength();

  if(isNuPRISM){
    if(ismPMT){
      detRadius = 371.;
      detLength = 1042.;
    }else{
      detRadius = 400.;
      detLength = 500.*2.;
    }
  }
  
  detRadius =  172.05-20.;
  detLength = (136.95-20.) * 2.;

  cout << "Detector radius " << detRadius << " Detector Length " << detLength << endl;

  std::map<int,TVector3> posMap;
  std::map<int,TVector3> dirMap;
  std::map<int,TVector3>::iterator iterPos;
  std::map<int,TVector3>::iterator iterDir;
  
  for (int i = 0; i < geoWCSim->GetWCNumPMT(); i++){
    posMap[geoWCSim->GetPMT(i).GetTubeNo()] = TVector3(geoWCSim->GetPMT(i).GetPosition(0),
						       geoWCSim->GetPMT(i).GetPosition(1),
						       geoWCSim->GetPMT(i).GetPosition(2));

    dirMap[geoWCSim->GetPMT(i).GetTubeNo()] = TVector3(geoWCSim->GetPMT(i).GetOrientation(0),
						       geoWCSim->GetPMT(i).GetOrientation(1),
						       geoWCSim->GetPMT(i).GetOrientation(2));
    
  }

  // Set up detector transformation
  TVector3 DetZinGlbl(0.,0.,1.);// Detector local Z-axis in global coordinate system
  TVector3 DetXinGlbl(1.,0.,0.);// Detector local X-axis in global coordinate system

  float WCXRotation[3] = {1.,0.,0.};
  float WCYRotation[3] = {0.,1.,0.};
  float WCZRotation[3] = {0.,0.,1.};
  
  TTree * fTSet = (TTree*) fWCSim->Get("Settings");
  if (!fTSet){
    std::cout<<"There is no Settings TTree in "<<rootFileNames.at(0).c_str()<<std::endl;
  } else {
    TBranch *bWCXRotation  = fTSet->GetBranch("WCXRotation");
    if (bWCXRotation) {
      bWCXRotation->SetAddress(WCXRotation);
    } else {
      WCXRotation[0] = 1.;
      WCXRotation[1] = 0;
      WCXRotation[2] = 0;
    }
    TBranch *bWCYRotation  = fTSet->GetBranch("WCYRotation");
    if (bWCYRotation) {
      bWCYRotation->SetAddress(WCYRotation);
    } else {
      WCYRotation[0] = 0;
      WCYRotation[1] = 1.;
      WCYRotation[2] = 0;
    }
    TBranch *bWCZRotation  = fTSet->GetBranch("WCZRotation");
    if (bWCZRotation) {
      bWCZRotation->SetAddress(WCZRotation);
    } else {
      WCZRotation[0] = 0;
      WCZRotation[1] = 0;
      WCZRotation[2] = 1.;
    }
    fTSet->GetEntry(0);
  }

  TransDetToGlbl.SetXYZ(geoWCSim->GetWCOffset(0),geoWCSim->GetWCOffset(1),geoWCSim->GetWCOffset(2));
  
  DetZinGlbl.SetXYZ(WCZRotation[0],WCZRotation[1],WCZRotation[2]);
  DetXinGlbl.SetXYZ(WCXRotation[0],WCXRotation[1],WCXRotation[2]);
  RotDetToGlbl.SetZAxis(DetZinGlbl,DetXinGlbl);// Rotation matrix which transforms detector local frame to global frame

  for (int i=0; i<3; i++){
      cout << "ROT matrix" << endl;
      cout << RotDetToGlbl(i, 0) << "  " << RotDetToGlbl(i, 1) << "  " << RotDetToGlbl(i, 2) << endl;
      cout << "----------" << endl;
  }

  delete tWCSimGeo;
  delete geoWCSim;
  delete fWCSim;
  
  cout << "Done reading geometry" << endl;
  
  TChain * WCSimChain = new TChain("sttree");

  for (int i = 0; i < rootFileNames.size(); i++) WCSimChain->Add(rootFileNames.at(i).c_str());

  vector<double> r;
  double dr;
  
  if (!isNuPRISM){
    cout << "Good bining" << endl;
    dr = 50;
    r.push_back(100);
    r.push_back(200);
    r.push_back(300);
    r.push_back(400);
    r.push_back(500);
    r.push_back(600);
    r.push_back(700);
    r.push_back(800);
    r.push_back(900);
    r.push_back(1000);
    r.push_back(1100);
    r.push_back(1200);
    r.push_back(1300);
    r.push_back(1400);
    r.push_back(1500);
  } else {
    cout << "Bad bining" << endl;
    dr = 50;
    r.push_back(50);
    r.push_back(100);
    r.push_back(150);
    r.push_back(200);
    r.push_back(250);
    r.push_back(300);
  }
  
  
  cout << "Sorted out radii" << endl;

  TFile * outRootFile = new TFile(outFName.c_str(), "RECREATE");
  
  std::vector<TH1F*> hAllPe;
  std::vector<TH1F*> hBarrelPe;
  std::vector<TH1F*> hEndCapPe;

  char nbuff[500];
  char buff[500];
  int nBins = 25;
  TH1F * temp;
  
  for (int i =0; i < r.size(); i++){

    sprintf(nbuff, "angRespAll_%i", (int) r.at(i));
    sprintf(buff, "Angular response function");
    temp = new TH1F(nbuff, buff, nBins, 0, 1);
    hAllPe.push_back(temp);
    hAllPe.back()->Sumw2();
    
    sprintf(nbuff, "angRespBarrel_%i", (int) r.at(i));
    //    sprintf(buff, "Angular response function at %i cm (only barrel PMTs)", (int) r.at(i).);
    sprintf(buff, "Angular response function(only barrel PMTs)");
    temp = new TH1F(nbuff, buff, nBins, 0, 1);
    hBarrelPe.push_back(temp);
    hBarrelPe.back()->Sumw2();
    
    sprintf(nbuff, "angRespEndCap_%i", (int) r.at(i));
    //    sprintf(buff, "Angular response function at %i cm(only end-cap PMTs)", (int) r.at(i).);
    sprintf(buff, "Angular response function(only endcap PMTs)");
    temp = new TH1F(nbuff, buff, nBins, 0, 1);
    hEndCapPe.push_back(temp);
    hEndCapPe.back()->Sumw2();
  }

  cout << "Newed histograms" << endl;
  
  float srcpos_[3];
  float oppos_[3];
  int ihPMT_;
  int isct_;

  WCSimChain->SetBranchAddress("srcpos", &srcpos_);
  WCSimChain->SetBranchAddress("oppos", &oppos_);
  WCSimChain->SetBranchAddress("ihPMT",  &ihPMT_);
  WCSimChain->SetBranchAddress("isct",   &isct_);
  cout << "Sets addresses" << endl;

  TVector3 thisTrackPos;
  TVector3 thisTrackPosTemp;
  TVector3 thisPMTdir;
  TVector3 thisPMTpos;

  TVector3 relPos;

  cout << " About to start loop " << endl;
  cout << "Number of entries in chain is " << WCSimChain->GetEntries() << endl;

  bool inFiducialZ;
  bool inFiducialR;
  
  for (int evt = 0; evt < WCSimChain->GetEntries(); evt++){
  //  for (int evt = 0; evt < 10 ; evt++){
    WCSimChain->GetEntry(evt);
    
    if(isct_ != 0) continue;

    if(!(evt%10000)) cout << rootFileNames.at(0) << "\t" << evt << std::endl;
    
    //    thisTrackPosTemp.SetXYZ(srcpos_[0]/10.,
    //			    srcpos_[1]/10.,
    //			    srcpos_[2]/10.);

    thisTrackPosTemp.SetXYZ(oppos_[0]/10.,
			    oppos_[1]/10.,
			    oppos_[2]/10.);
    
    thisTrackPos = TfmDetToGlblCoord(thisTrackPosTemp, true, true);
    
    int tubeNo  = ihPMT_;
    int totalPe = 1;
      
    iterPos = posMap.find(tubeNo);
    iterDir = dirMap.find(tubeNo);

    if (iterPos == posMap.end() || iterDir == dirMap.end()){
      cout << "ERROR PMT POSITION OR DIRECTION NOT FOUND " << tubeNo << endl;
      return;
    }

    
    thisPMTpos = TfmDetToGlblCoord(posMap[tubeNo], true,  true);
    thisPMTdir = TfmDetToGlblCoord(dirMap[tubeNo], false, true);
    
    //    std::cout <<  "PMT " << thisPMTpos.X() << " " << thisPMTpos.Y() << " " << thisPMTpos.Z() << endl;
    //    std::cout <<  "oppos " << thisTrackPos.X() << " " << thisTrackPos.Y() << " " << thisTrackPos.Z() << endl;
    
    relPos = thisTrackPos - thisPMTpos;

    //    cout << "EVT " << evt << "\tR = " << relPos.Mag() << "\tcos(eta) = " << relPos.Dot(thisPMTdir)/relPos.Mag() << "\tthisTrackPos = ( " << thisTrackPos.X() << "," << thisTrackPos.Y() << "," << thisTrackPos.Z() <<" )" << "\tthisPMTpos = ( " << thisPMTpos.X() << "," << thisPMTpos.Y() << "," << thisPMTpos.Z() << ")" <<  endl;
    
    for (int shell = 0; shell < r.size(); shell++){

      //      cout << "Test Shell " << relPos.Mag()  << " " << ( r.at(shell) - dr ) << " " << ( r.at(shell)+dr) << endl;
      if (relPos.Mag() >= ( r.at(shell) - dr ) && relPos.Mag() < ( r.at(shell)+dr)){	// In spherical shell
	
        inFiducialZ = (fabs(thisPMTpos.Z()) + (r.at(shell)+dr)) <= (detLength/2.);
        //	cout << "Test Z " << (fabs(thisPMTpos.Z()) + (r.at(shell)+dr)) << " " << (detLength/2.) << endl;
        inFiducialR = (thisPMTpos.Perp() + (r.at(shell)+dr)) <= detRadius;
        //	cout << "Test R " << (thisPMTpos.Perp() + (r.at(shell)+dr)) << " " <<  detRadius << endl;
        //	cout << "fiducial Z " << inFiducialZ << " fiducial R " << inFiducialR << endl;
        
        if (inFiducialZ && inFiducialR){
          cout << "inFiducialZ and inFiducialR" << endl;
          cout << thisPMTpos.X() << " " <<thisPMTpos.Y() << " "<< thisPMTpos.Z() << endl;
          cout << thisPMTpos.Perp() <<  " " << detRadius << endl;
          return;
	      }
	  
        if (inFiducialZ || inFiducialR){ // Make sure spherical shell is fully contained within "fiducial" volume 

          hAllPe.at(shell)->Fill(relPos.Dot(thisPMTdir)/relPos.Mag(), totalPe);
          if (inFiducialZ) { // Must be a barrel PMT
            hBarrelPe.at(shell)->Fill(relPos.Dot(thisPMTdir)/relPos.Mag(), totalPe);
          } else { // It's an endcap PMT
            hEndCapPe.at(shell)->Fill(relPos.Dot(thisPMTdir)/relPos.Mag(), totalPe);
          }
          break;
        }
        break;
      }
    }
  }
  

  TLegend * leg = new TLegend(0.5, 0.1, 0.9, 0.5);


  for (int i = 0 ; i < r.size(); i++){
    sprintf(buff, "R = %i cm", (int) r.at(i));
    leg->AddEntry(hAllPe.at(i), buff, "f");
    if (hAllPe.at(i)->Integral())    hAllPe.at(i)->   Scale(1/hAllPe.at(i)->   GetBinContent(nBins));
    if (hBarrelPe.at(i)->Integral()) hBarrelPe.at(i)->Scale(1/hBarrelPe.at(i)->GetBinContent(nBins)); 
    if (hEndCapPe.at(i)->Integral()) hEndCapPe.at(i)->Scale(1/hEndCapPe.at(i)->GetBinContent(nBins));
  }

//  TCanvas * c1 = new TCanvas();
//  hAllPe.at(0)->SetFillColor(kBlack);
//  hAllPe.at(0)->Draw("E3");
//  hAllPe.at(1)->SetLineColor(kAzure+2);
//  hAllPe.at(1)->SetFillColor(kAzure+2);
//  hAllPe.at(1)->Draw("SAMEE3");
//  hAllPe.at(2)->SetLineColor(kGreen+2);
//  hAllPe.at(2)->SetFillColor(kGreen+2);
//  hAllPe.at(2)->Draw("SAMEE3");
//
//  leg->Draw("SAME");
//
//  TCanvas * c2 = new TCanvas();
//  hBarrelPe.at(0)->SetFillColor(kBlack);
//  hBarrelPe.at(0)->Draw("E3");
//  hBarrelPe.at(1)->SetLineColor(kAzure+2);
//  hBarrelPe.at(1)->SetFillColor(kAzure+2);
//  hBarrelPe.at(1)->Draw("SAMEE3");
//  hBarrelPe.at(2)->SetLineColor(kGreen+2);
//  hBarrelPe.at(2)->SetFillColor(kGreen+2);
//  hBarrelPe.at(2)->Draw("SAMEE3");
//
//  leg->Draw("SAME");
//  
//  TCanvas * c3 = new TCanvas();
//  hEndCapPe.at(0)->SetFillColor(kBlack);
//  hEndCapPe.at(0)->Draw("E3");
//  hEndCapPe.at(1)->SetLineColor(kAzure+2);
//  hEndCapPe.at(1)->SetFillColor(kAzure+2);
//  hEndCapPe.at(1)->Draw("SAMEE3");
//  hEndCapPe.at(2)->SetLineColor(kGreen+2);
//  hEndCapPe.at(2)->SetFillColor(kGreen+2);
//  hEndCapPe.at(2)->Draw("SAMEE3");
//
//  leg->Draw("SAME");
//
  outRootFile->Write();
  outRootFile->Close();
}

TVector3 TfmDetToGlblCoord(TVector3 tmpVct, bool flgPos, bool flgInverse){
  
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
