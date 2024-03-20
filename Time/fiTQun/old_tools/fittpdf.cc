//script to normalize the cherenkov profiles and put them into one file

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TMath.h>

#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

using namespace std;

const int nthrdpars=10;//number of momentum fit parameters
const int npars=7;//number of logmu fit parameters

bool shiftLeft = true;
bool useEMG = false;
bool useFitSlices = true;

Double_t emg(Double_t *x, Double_t *par){
  double xx = x[0];
  double k = par[0];
  double mu = par[1];
  double sigma = par[2];
  double lambda = par[3];
  
  return k*lambda/2 * TMath::Exp(lambda/2*(2*mu+lambda*pow(sigma,2)-2*xx)) * TMath::Erfc( (mu + lambda*pow(sigma,2) - xx) / (sqrt(2)*sigma) );
}

TString PolyFormula(int nord) {
  TString strPoly="[0]+[1]*x";
  for (int i=2; i<=nord; i++) {
    strPoly+=Form("+[%d]*x^%d",i,i);
  }
  
  return strPoly;
}

void fittpdf(int PID, bool flogfit=false, bool fNoErrorbars=true){

  TFitResultPtr ptrres[2],ptrprres[2];
  const int nmommax=200;
  double momfitrang[2]={0.,7000.};
  double mnpar[npars],sgpar[npars],thdpr[2][npars][nthrdpars];
  
  TString strParNam[2]={"mn","sg"};
  
  Double_t armom[nmommax],arpar[2][npars][nmommax],arparerr[2][npars][nmommax];
  TH1D *hmnsg[2][nmommax];
  TH1D *hmeantmp,*hsigmtmp, *hchi2tmp;
  TH2D *hist2d;
  TGraphErrors *gtpar[2][npars];
  int ntmp;
  
  double pofst=0.;
  if (PID==11) {
    momfitrang[0]=100.;
    momfitrang[1]=1000.;
    pofst=0.;
  }
  else if (PID==13) {
    momfitrang[0]=180.;
    momfitrang[1]=2000.;
    pofst=100.;
  }
  else if (PID==211) {
    momfitrang[0]=220.;
    momfitrang[1]=2000.;
    pofst=140.;
  }
  else if (PID==321) {
    momfitrang[0]=600.;
    momfitrang[1]=5000.;
    pofst=530.;
  }
  else if (PID==2212) {
    momfitrang[0]=1150.;
    momfitrang[1]=5000.;
    pofst=1020.;
  }
  
  if (flogfit) {
    momfitrang[0]=log(momfitrang[0]-pofst);
    momfitrang[1]=log(momfitrang[1]-pofst);
  }
  else {
    pofst=0.;
  }
  
  TFile *fitmp = new TFile(Form("%d_tpdfhist.root",PID));
  if (fitmp->IsZombie()) exit(-1);

  TH3D *h3d = (TH3D*)fitmp->Get("hist_tpdf");
  double tcmin=h3d->GetXaxis()->GetBinLowEdge(1);
  double tcmax=h3d->GetXaxis()->GetBinLowEdge(h3d->GetNbinsX()+1);

  TH2D * hmean_step0 = (TH2D*)h3d->Project3D("yz")->Clone("hmean_step0");
  TH2D * hmean_step1 = (TH2D*)h3d->Project3D("yz")->Clone("hmean_step1");

  TH2D * hsigm_step0 = (TH2D*)h3d->Project3D("yz")->Clone("hsigm_step0");
  TH2D * hsigm_step1 = (TH2D*)h3d->Project3D("yz")->Clone("hsigm_step1");
  
  TFitResultPtr timeResFitPointer;
  TF1 *myGausFit = new TF1("myGausfit","gaus", tcmin, tcmax);
  myGausFit->SetParName(0, "Const");
  myGausFit->SetParName(1, "Mean");
  myGausFit->SetParName(2, "Sigma");
  
  TF1 *fEMG = new TF1("fEMG",emg,tcmin,tcmax,4);
  fEMG->SetParName(0, "Const");
  fEMG->SetParName(1, "Mean");
  fEMG->SetParName(2, "Sigma");
  fEMG->SetParName(3, "Lambda");
  
  
  int nmom=h3d->GetNbinsZ();
  if (PID==13) nmom = 24;
  if (PID==211) nmom = 24;

  // For each momentum: 
  //     Step 0) Gaussian fit of t-residual for each log10(mu)
  //     Step 1) Polynomial fit of mean and sigma vs log10(mu)
  for (int k=0; k<nmom; k++) {
    // int imom=(int)(h3d->GetZaxis()->GetBinLowEdge(k+1)+0.5);
    int imom=(int)(h3d->GetZaxis()->GetBinCenter(k+1));
    cout << "############################" << endl;
    cout << "MOM=" << imom << endl;
    cout << "############################" << endl;
    armom[k]=imom;
    if (flogfit) armom[k]=log(armom[k]-pofst);
    h3d->GetZaxis()->SetRange(k+1,k+1);
    hist2d = (TH2D*)h3d->Project3D("yx");
    hist2d->FitSlicesX(); //Step 0) Gaussian fit of t-residual for each log10(mu)
    hmeantmp = (TH1D*)gDirectory->Get("hist_tpdf_yx_1");
    hsigmtmp = (TH1D*)gDirectory->Get("hist_tpdf_yx_2");
    hchi2tmp = (TH1D*)gDirectory->Get("hist_tpdf_yx_chi2");
    ntmp=hmeantmp->GetNbinsX();

    if(useFitSlices){
      for (int muBin = 1; muBin <= ntmp; muBin++){
        hmean_step0->SetBinContent( k+1,muBin,hmeantmp->GetBinContent(muBin));
        hsigm_step0->SetBinContent( k+1,muBin,hsigmtmp->GetBinContent(muBin));
      }
    }else{ cout << "removed by Gonzalo to simplify code" << endl; exit(-1);}

    for (int ibin=1; ibin<=ntmp; ibin++) {
      double dtmp=hmeantmp->GetBinContent(ibin);
      double chi2tmp=hchi2tmp->GetBinContent(ibin);
      //      if (!(dtmp>tcmin && dtmp<tcmax) || chi2tmp > 20) {
      if (!(dtmp>tcmin && dtmp<tcmax)) {
        hmeantmp->SetBinContent(ibin,0.);
        hmeantmp->SetBinError(ibin,0.);
      }
    }
    hmnsg[0][k] = new TH1D(*hmeantmp);
    hmnsg[1][k] = new TH1D(*hsigmtmp);
    hmnsg[0][k]->SetName(Form("hmean_%d",imom));
    hmnsg[1][k]->SetName(Form("hsigm_%d",imom));

    // Step 1) Polynomial fit of mean and sigma vs log10(mu)
    for (int idx=0; idx<2; idx++) {
      ptrres[idx]=hmnsg[idx][k]->Fit(Form("pol%d", npars-1),"SQ");
      int ires=ptrres[idx];
      // cout << ptrres[idx]->Ndf() << endl;
      if (ires) {
        cout << "Fit abnormally terminated!!! ires: " << ires << " idx " << idx << " imom " << imom <<  endl;
        exit(-1);
      } else if (ptrres[idx]->Chi2()/ptrres[idx]->Ndf() > 25.){
        cout << "Bad fiT! ires: " << ires << " idx " << idx << " imom " << imom << " chi2/ndf " << (ptrres[idx]->Chi2()/ptrres[idx]->Ndf()) <<  endl;
      }
    }
    
    for (int i=0; i<npars; i++) {
      arpar[0][i][k]=ptrres[0]->Parameter(i);//corrected time mean parameters
      arpar[1][i][k]=ptrres[1]->Parameter(i);//corrected time RMS parameters
      arparerr[0][i][k]=ptrres[0]->ParError(i);//corrected time mean parameters
      arparerr[1][i][k]=ptrres[1]->ParError(i);//corrected time RMS parameters
      //        cout << arpar[0][i][k] << ", " << arpar[1][i][k] << endl;
    }
  }

  // Save fit step 1 results
  for (int k = 0; k < nmom; k++){
    for (int muBin = 1; muBin < ntmp; muBin++){
      double testMean = 0.;
      double testSigma = 0.;
      double mupow = 1.;
      for (int i = 0; i < npars; i++){
        testMean += arpar[0][i][k]*mupow;
        testSigma += arpar[1][i][k]*mupow;
        mupow *= hmeantmp->GetXaxis()->GetBinLowEdge(muBin);
      }
      hmean_step1->SetBinContent( k+1,muBin, testMean);
      hsigm_step1->SetBinContent( k+1,muBin, testSigma);
    }
  }
  
  cout << "############################" << endl;
  cout << "FIT FOR EACH MOMENTUM DONE"   << endl;
  cout << "############################" << endl;

  TFile *fout = new TFile(Form("%d_tpdfpar.root",PID),"RECREATE");
  
  TH2D *htpdfpar[2];
  htpdfpar[0] = new TH2D("htpdfparmn","Time pdf mean fit parameters",npars,0,npars,nthrdpars,0,nthrdpars);
  htpdfpar[1] = new TH2D("htpdfparsg","Time pdf sigma fit parameters",npars,0,npars,nthrdpars,0,nthrdpars);
  
  if (!(momfitrang[0]>armom[0])) momfitrang[0]=armom[0];
  if (!(momfitrang[1]<armom[nmom-1])) momfitrang[1]=armom[nmom-1];
  
  for (int i=0; i<npars; i++) {
    cout << "############################" << endl;
    cout << "PAR  " << i << endl;
    cout << "############################" << endl;
    for (int idx=0; idx<2; idx++) {
      if (fNoErrorbars) {
        gtpar[idx][i] = new TGraphErrors(nmom,armom,arpar[idx][i],NULL,NULL);
      }
      else {
        gtpar[idx][i] = new TGraphErrors(nmom,armom,arpar[idx][i],NULL,arparerr[idx][i]);
      }
      gtpar[idx][i]->SetName(Form("gtc%spar_%d",strParNam[idx].Data(),i));
      
      int nthrdtmp=nthrdpars+1;
      int ires;
      do {
        nthrdtmp--;
        
        cout << "nthrd=" << nthrdtmp << endl;
        if (nthrdtmp<6) {
          cout << "Fit abnormally terminated!!!" << endl;
          exit(-1);
        }
        
        ptrprres[idx]=gtpar[idx][i]->Fit(Form("pol%d", nthrdtmp-1),"SQ","",momfitrang[0],momfitrang[1]);
        ires=ptrprres[idx];
        
      } while (ires);
      
      for (int j=0; j<nthrdpars; j++) {
        if (j<nthrdtmp) {
          thdpr[idx][i][j]=ptrprres[idx]->Parameter(j);
        }
        else {
          thdpr[idx][i][j]=0.;
        }
        
        htpdfpar[idx]->SetBinContent(i+1,j+1,thdpr[idx][i][j]);
      }
      
      gtpar[idx][i]->SetMarkerStyle(2);
      gtpar[idx][i]->Write();
      delete gtpar[idx][i];
    }
    
  }
  
  for (int j=0; j<nmom; j++) {
    hmnsg[0][j]->Write();
    hmnsg[1][j]->Write();
    delete hmnsg[0][j];
    delete hmnsg[1][j];
  }
  
  TH1D *htpdfinfo = new TH1D("htpdfinfo","",10,0,10);
  htpdfinfo->SetBinContent(1,npars);
  htpdfinfo->SetBinContent(2,h3d->GetYaxis()->GetBinCenter(1));
  htpdfinfo->SetBinContent(3,h3d->GetYaxis()->GetBinCenter(h3d->GetNbinsY()));
  htpdfinfo->SetBinContent(4,nthrdpars);
  htpdfinfo->SetBinContent(5,momfitrang[0]);
  htpdfinfo->SetBinContent(6,momfitrang[1]);
  htpdfinfo->SetBinContent(7,pofst);
  
  htpdfinfo->Write();
  delete htpdfinfo;

  hmean_step0->Write();
  hsigm_step0->Write();
  hmean_step1->Write();
  hsigm_step1->Write();
  
  htpdfpar[0]->Write();
  htpdfpar[1]->Write();
  delete htpdfpar[0];
  delete htpdfpar[1];
  
  delete fout;
  
  delete fitmp;

}
