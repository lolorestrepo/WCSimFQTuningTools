// There seems to be a nasty memory leak somewhere in the loop...

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TLine.h>
#include <TCanvas.h>

#include "fiTQun.h"
#include "fiTQun_shared.h"
#include "WCSimWrap.h"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"

int main (int argc, char ** argv){
  
  if( argc != 4 ){
    std::cout << "Arguments are: list of mu, base path to wcsim files, output filename" << std::endl;
    return -2;
  }  
  std::cout << "Opening file " << argv[1] << std::endl;
  
  std::ifstream mus;
  mus.open(argv[1]);
  if (!mus.is_open()) return -1;
  
  std::string mustr;

  char rootFname[500];
  double trueMu;
  
  std::vector<double> trueMuVec;
  std::vector<double> trueMuErrVec;
  std::vector<double> recMuVec;
  std::vector<double> recMuErrVec;
  std::vector<double> recRatioVec;
  std::vector<double> recRatioErrVec;

  std::vector<TH1F*> histVec;

  while ( getline (mus, mustr) ){
    sprintf(rootFname, "%s/Mu_%s/%s.root", argv[2], mustr.c_str(), mustr.c_str());
    trueMu = atof(mustr.c_str());
    std::cout << trueMu << " "<<rootFname << std::endl;
    
    WCSimWrap* wc = WCSimWrap::Get( rootFname );
    
    fiTQun * thefit = new fiTQun( wc->NPMT() );

    int iScat=0;//no scattered light predicted charge
    int fLoadflg[nPID] = {0}; // Don't load any cherenkov profiles?
    
    fiTQun_shared::Get()->SetScatflg(iScat);
    
    thefit->ReadSharedParams(iScat,fLoadflg,true,true,1);
    
    fiTQun_shared::Get()->SetQEEff(0.1);
    fiTQun_shared::Get()->SetPhi0(-1,1.);//remove dependence on tuning const.
    fiTQun_shared::Get()->SetDarkRate(0.);
    
    char histTitle[500];
    char histName[500];
    
    sprintf(histName, "fitMuHist_%f", trueMu);
    sprintf(histTitle, "True #mu = %f p.e.; Fit #mu [p.e.]", trueMu);


    TH1F * hist;
    if (trueMu <=1 ){
      hist = new TH1F(histName, histTitle, 2400, 0.4*trueMu, 1.6*trueMu);
    } else if (trueMu < 100) {
      hist = new TH1F(histName, histTitle, 600, 0.9*trueMu, 1.1*trueMu);
    } else {
      hist = new TH1F(histName, histTitle, 600, 0.98*trueMu, 1.02*trueMu);
    }

    for (int nevt = 0; nevt < wc->NEvt(); nevt++){
      
      wc->LoadEntry(nevt);
      
      thefit->InitEvent(nevt);
      
      thefit->ClrInfo(0);
      
      thefit->SetTimeWindow(0);
      
      hist->Fill(thefit->FitMu());
    }
    
    
    hist->Fit("gaus", "L");
    
    histVec.push_back(hist);

    TF1 * fun = hist->GetFunction("gaus");
    
    trueMuVec.push_back(trueMu);
    trueMuErrVec.push_back(0.);
  
    recMuVec.push_back(fun->GetParameter(1));
    recMuErrVec.push_back(fun->GetParameter(2));

    recRatioVec.push_back(fun->GetParameter(1)/trueMu);
    recRatioErrVec.push_back(fun->GetParameter(2)/trueMu);
    

    delete thefit;
    delete wc;
  }
  
  
  
  gROOT->SetBatch();

  TFile * outFile = new TFile(argv[3], "RECREATE");
  
  for (std::vector<TH1F*>::iterator vecIt = histVec.begin();
       vecIt != histVec.end(); vecIt++) (*vecIt)->Write();
  

  TGraphErrors * gMuBias = new TGraphErrors (trueMuVec.size(), &(trueMuVec[0]), &(recRatioVec[0]), &(recMuErrVec[0]), &(recRatioErrVec[0]));
  gMuBias->SetTitle("Charge PDF bias;#mu_{True} [p.e.];#mu_{Rec} / #mu_{True}");

  gROOT->SetBatch();

  TCanvas * c1 = new TCanvas();
  gMuBias->Draw("AP");
  c1->SetLogx();
  
  TLine * l1 = new TLine( gMuBias->GetHistogram()->GetXaxis()->GetXmin(), 1., gMuBias->GetHistogram()->GetXaxis()->GetXmax(), 1. );

  //TLine * l1 = new TLine( 0.1, 1., 800, 1. );
  l1->SetLineColor(kRed);
  l1->SetLineStyle(2);
  l1->Draw("SAME");
  gMuBias->Draw("SAMEP");

  int lastindex = std::string(argv[3]).find_last_of("."); 
  std::string rawname = std::string(argv[3]).substr(0, lastindex);

  c1->RedrawAxis();
  c1->SaveAs((rawname+".pdf").c_str());

  gMuBias->Write();
  outFile->Write();
  outFile->Close();


}
