//script to normalize the cherenkov profiles and put them into one file

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraphErrors.h>
#include <TFile.h>

using namespace std;
const double pi=3.1415926535898;
const double sqrt2pi=sqrt(2.*pi);

Double_t ftpdf(Double_t *x, Double_t *par) {
  return exp(-(x[0]-par[0])*(x[0]-par[0])/2./par[1]/par[1])/sqrt2pi/par[1];
}

int plottpdf(int PID, int momidx, int muidx){
  
  TFile *fitmp = new TFile(Form("%d_tpdfhist.root",PID));
  if (fitmp->IsZombie()) exit(-1);
  
  TH3D *h3d = (TH3D*)fitmp->Get("hist_tpdf");
  
  int nmu=h3d->GetNbinsY();
  int nmom=h3d->GetNbinsZ();
  
  if (momidx<0 || momidx>=nmom) exit(-1);
  if (muidx<0 || muidx>=nmu) exit(-1);
  
  int imom=(int)(h3d->GetZaxis()->GetBinLowEdge(momidx+1)+0.5);
  double mom=(double)imom;
  double logmu=h3d->GetYaxis()->GetBinCenter(muidx+1);
  cout << "p=" << imom << endl;
  cout << "logmu=" << logmu << endl;
  TH1D *hist1d = (TH1D*)h3d->ProjectionX("hproj",muidx+1,muidx+1,momidx+1,momidx+1);
  
  TString nmParticle="";
  if (PID==11) {
    nmParticle="e^{-}";
  }
  else if (PID==13) {
    nmParticle="#mu^{-}";
  }
  
  hist1d->Scale(1./hist1d->Integral("width"));
  
  gStyle->SetOptStat(0);
  hist1d->SetTitle(Form("PDG=%d, p=%dMeV/c, log#mu=%.3f",PID,imom,logmu));
  hist1d->GetXaxis()->SetTitle("t^{res} [ns]");
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TCanvas *c = new TCanvas("cnvs","",0,0,400,400);
  
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.05);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.15);
  
  hist1d->SetLineWidth(2);
  hist1d->Draw();
  
  hist1d->GetXaxis()->SetRangeUser(-25.,50.);
  
  TFile *fpar = new TFile(Form("%d_tpdfpar.root",PID));
  
  TH2D *htpdfpar[2];
  
  htpdfpar[0] = (TH2D*)fpar->Get("htpdfparmn");
  htpdfpar[1] = (TH2D*)fpar->Get("htpdfparsg");
  
  TH1D *htpdfinfo = (TH1D*)fpar->Get("htpdfinfo");
  double tpdfprang[2];
  double tpdfmurang[2];
  double pofst;
  
  int npars=(int)(htpdfinfo->GetBinContent(1)+0.5);
  tpdfmurang[0]=htpdfinfo->GetBinContent(2);
  tpdfmurang[1]=htpdfinfo->GetBinContent(3);
  int nthrdpars=(int)(htpdfinfo->GetBinContent(4)+0.5);
  tpdfprang[0]=htpdfinfo->GetBinContent(5);
  tpdfprang[1]=htpdfinfo->GetBinContent(6);
  pofst=htpdfinfo->GetBinContent(7);
    
  bool flogfit=false;
  if (fabs(pofst)>0.1) flogfit=true;
  
  TString strParNam[2]={"mn","sg"};
  TGraphErrors *gtpar[2][20];
  
  for (int i=0; i<npars; i++) {
    for (int idx=0; idx<2; idx++) {
      gtpar[idx][i] = (TGraphErrors*)fpar->Get(Form("gtc%spar_%d",strParNam[idx].Data(),i));
    }
  }
  
  if (logmu<tpdfmurang[0]) logmu=tpdfmurang[0];//ensure logmu is in the fitted range
  else if (logmu>tpdfmurang[1]) logmu=tpdfmurang[1];
  
  double momtmp=mom;
  if (flogfit) momtmp=log(momtmp-pofst);
  if (momtmp<tpdfprang[0]) momtmp=tpdfprang[0];
  else if (momtmp>tpdfprang[1]) momtmp=tpdfprang[1];
  
  
  double tpdfpar[2][2][10];
  
  for (int imupar=0; imupar<npars; imupar++) {
    tpdfpar[0][0][imupar]=gtpar[0][imupar]->Eval(momtmp);
    tpdfpar[0][1][imupar]=gtpar[1][imupar]->Eval(momtmp);
    
    tpdfpar[1][0][imupar]=0.;
    tpdfpar[1][1][imupar]=0.;
    double ppow=1.;
    for (int ippar=0; ippar<nthrdpars; ippar++) {
      tpdfpar[1][0][imupar]+=htpdfpar[0]->GetBinContent(imupar+1,ippar+1)*ppow;
      tpdfpar[1][1][imupar]+=htpdfpar[1]->GetBinContent(imupar+1,ippar+1)*ppow;
      ppow*=momtmp;
    }
  }
  
  for (int ifit=1; ifit<2; ifit++) {
    double* mnpar=tpdfpar[ifit][0];
    double* sgpar=tpdfpar[ifit][1];
    
    double tcmean=0.;//corrected time mean
    double tsigm=0.;//corrected time RMS
    double mupow=1.;
    for (int i=0; i<npars; i++) {
      tcmean+=mnpar[i]*mupow;
      tsigm+=sgpar[i]*mupow;
      mupow*=logmu;
    }
    
    TF1 *ffit = new TF1(Form("ffit_%d",ifit),ftpdf,hist1d->GetXaxis()->GetBinLowEdge(1),hist1d->GetXaxis()->GetBinLowEdge(hist1d->GetNbinsX()+1),2);
    ffit->SetParameter(0,tcmean);
    ffit->SetParameter(1,tsigm);
    
    ffit->SetLineColor(1+ifit);
//    ffit->SetLineStyle(1+ifit);
    ffit->SetNpx(1000);
    ffit->Draw("SAME");
  }
  hist1d->Draw("SAME");
  
  hist1d->GetXaxis()->SetTitleSize(0.06);
  hist1d->GetYaxis()->SetTitleSize(0.06);
  hist1d->GetXaxis()->SetTitleOffset(0.8);
  hist1d->GetYaxis()->SetTitleOffset(1.0);
  
  TLatex *ltx = new TLatex;
  ltx->SetTextAlign(11);
  ltx->SetNDC();
  
  ltx->SetTextSize(0.09);
  ltx->DrawLatex(0.54,0.83,Form("%s",nmParticle.Data()));
  ltx->SetTextSize(0.07);
  ltx->DrawLatex(0.55,0.75,Form("p=%dMeV/c",imom));
  ltx->DrawLatex(0.55,0.67,Form("#mu=%5.2f p.e.",pow(10.,logmu)));
  
//  c->Print(Form("tpdf_%d_%d_%d.eps",PID,imom,muidx));
  c->Print(Form("tpdf_%d_%d_%d.pdf",PID,imom,muidx));
  
}
