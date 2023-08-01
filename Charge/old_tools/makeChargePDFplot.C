#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include "TH1D.h"
#include "TFile.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TTree.h"
#include "/pbs/home/g/gdiazlop/Software/WCSim/install/include/WCSimRootEvent.hh"
#include "/pbs/home/g/gdiazlop/Software/WCSim/install/include/WCSimRootGeom.hh"


void makeChargePDFplot(std::string fNameWCSim, std::string fNameOut){
  
  // WCSim

  gSystem->Load("/pbs/home/g/gdiazlop/Software/WCSim/install/lib/libWCSimRoot.so");
  std::cout << "WCSim ROOT lib loaded" << std::endl;

  TFile * fWCSim = new TFile (fNameWCSim.c_str());
  TTree * tWCSim = (TTree*) fWCSim->Get("wcsimT");
  TTree * tWCSimGeo = (TTree*) fWCSim->Get("wcsimGeoT");

  WCSimRootGeom* WCSimGeom = new WCSimRootGeom();
  TBranch* bWCSimRootGeom = tWCSimGeo->GetBranch("wcsimrootgeom");
  bWCSimRootGeom->SetAddress(&WCSimGeom);
  tWCSimGeo->GetBranch("wcsimrootgeom")->SetAutoDelete(kTRUE);
  
  tWCSimGeo->GetEntry(0);

  int nPMTsGeo = WCSimGeom-> GetWCNumPMT();
  cout << nPMTsGeo << endl;

  WCSimRootEvent* evtWCSim = new WCSimRootEvent();
  
  TBranch * bWCSimRootEvent = tWCSim->GetBranch("wcsimrootevent");
  bWCSimRootEvent->SetAddress(&evtWCSim);
  tWCSim->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  WCSimRootTrigger      * trigWCSim;
  TClonesArray          * digiHitArr;
  //  WCSimRootCherenkovHit * digiHit;

  // Shimpei's binning
  double bEtmp=-1;
  double qbinEdg[500];
  int nqbins;

  ifstream fin("../qbins_wcte.txt");
  int i=0;
  while (1) {
    fin >> qbinEdg[i];

    if (bEtmp>=qbinEdg[i]) {
      cout << "qbins.txt: q must be in ascending order!" << endl;
      exit(-1);
    }
    if ( fin.eof() || qbinEdg[i] == 120 ) break;

    bEtmp=qbinEdg[i];
    i++;
  }
  nqbins=i;
  fin.close();

  

  //  TH1F * hQ  = new TH1F("hQ", "hQ;Q", 100, 0, 10);
  TFile * fOut = new TFile(fNameOut.c_str(), "recreate");

  TH1D * hcpdf2 = new TH1D("hchpdf2","Old PMT Charge PDF",nqbins,qbinEdg);
  hcpdf2->Sumw2();
  TH1D * hcpdf3 = new TH1D("hchpdf3","Old PMT Charge PDF",nqbins,qbinEdg);
  hcpdf3->Sumw2();
  
  TH1D* hctr = new TH1D("hctr","Total # of active PMTs",10,0.5,10.5); // FIX!!
  hctr->Sumw2();

  int nPMTs;

  //  TH1F * hPe = new TH1F("hPe", "hPe;Pe", 100, 0, 10);
  
  vector<int> pmtIDs;
  
  for (int evt = 0; evt < tWCSim->GetEntries(); evt++){
    

    tWCSim->GetEntry(evt);
    int nSubEvts = evtWCSim->GetNumberOfSubEvents();
    
    trigWCSim = evtWCSim->GetTrigger(0);
    //    digiHitArr = trigWCSim->GetCherenkovDigiHits();
    
    int nRawHits  = trigWCSim->GetNcherenkovhits() ;
    int nDigiHits = trigWCSim->GetNcherenkovdigihits();
    
    //    nPMTs = nRawHits;  // Commented out because not all effective PMT has raw hit since WCSim1.5
    nPMTs = nPMTsGeo;
    
    cout << endl << nRawHits << " raw and " << nDigiHits << " digitized hits in event " << evt << " n subevents " << nSubEvts << endl;
    
    for (int hit = 0; hit < nDigiHits; hit++){
      TObject *element = (trigWCSim->GetCherenkovDigiHits())->At(hit);
      
      WCSimRootCherenkovDigiHit * digiHit = dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
      
      
      //      if ( digiHit->GetTubeId() <= 1 || digiHit->GetTubeId() >= 11145 ) cout << digiHit->GetTubeId() << "\t";
      
      //      if (digiHit->GetTubeId() <= 1 ||
      //	  digiHit->GetTubeId() >= 11145){
      //
      //	cout << digiHit->GetTubeId() << "\t" << digiHit->GetQ() << endl;
      //	
      //	bool pmtfound = false;
      //	for (int i = 0; i < nPMTsGeo; i++){
      //	  if (WCSimGeom->GetPMT(i).GetTubeNo() == digiHit->GetTubeId()){
      //	    cout << "tube coords " << WCSimGeom->GetPMT(i).GetPosition(0) << " "<< WCSimGeom->GetPMT(i).GetPosition(1) << " "<< WCSimGeom->GetPMT(i).GetPosition(2) << endl;
      //	    pmtfound = true;
      //	  }
      //	}
      //	if (!pmtfound) cout << "TUBE ID " << digiHit->GetTubeId() << " NOT FOUND!!!!!!!!!!! " << endl;
      //
      //      }
      //      //      cout << "hit " << hit << "\tQ " << digiHit->GetQ() << "\tpe " << endl;
      hcpdf2->Fill(digiHit->GetQ());
      hcpdf3->Fill(digiHit->GetQ());
    }
    //      hQ->Fill(digiHit->GetQ());
    
    //      if ( find(pmtIDs.begin(), pmtIDs.end(), digiHit->GetTubeId()) != pmtIDs.end()){
    //	cout << "Duplicate pmtID " << digiHit->GetTubeId() << " in hit number " << hit << endl;
    //      }
    
    //      if ( digiHit->GetTubeId() < 
    
    
    //      pmtIDs.push_back(digiHit->GetTubeId());
    hctr->Fill(2, nPMTs);    
    hctr->Fill(3, nPMTs);    
  }
  //    pmtIDs.clear();
  
  //    for (int hit = 0; hit < nRawHits; hit++){
  //      TObject *element = (trigWCSim->GetCherenkovHits())->At(hit);
  //	  
  //      WCSimRootCherenkovDigiHit * rawHit =
  //	dynamic_cast<WCSimRootCherenkovDigiHit*>(element);
  //
  //      //      cout << "hit " << hit << "\tQ " << digiHit->GetQ() << "\tpe " << endl;
  //
  //      hQ->Fill(digiHit->GetQ());
  //    }
  //
  
  
  
  
  
  
  //  TCanvas * c1 = new TCanvas();
  //  hPe->Draw();
  //  TCanvas * c2 = new TCanvas();
  //  hQ->Draw();
  
  hcpdf2->Write();
  hcpdf3->Write();
  
  hctr->SetBinContent(10,hctr->GetBinContent(2)+hctr->GetBinContent(3));
  
  hctr->Write();
  
}

