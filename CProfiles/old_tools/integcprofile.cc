
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TGraph.h>
#include <TFile.h>

using namespace std;

void DoIntegral(TH2D *hg, double R0, double costh0, double *In)
{
  int nsbins=hg->GetNbinsY();//number of s bins
  double sbinw=hg->GetYaxis()->GetBinWidth(1);
  double val;
  double s,sn,costh,Rsx,Rsy,Rs;
  const double cosmx=1.-1e-12;//interpolation fails if costh==1.
  
  double ex=-costh0;
  double ey=sqrt(1.-costh0*costh0);
  
  for (int ns=0; ns<3; ns++) In[ns]=0.;
  for (int i=1; i<=nsbins; i++) {
    s=hg->GetYaxis()->GetBinCenter(i);
    Rsx=R0+ex*s;
    Rsy=ey*s;
    Rs=sqrt(Rsx*Rsx+Rsy*Rsy);
    costh=-(ex*Rsx+ey*Rsy)/Rs;
    if (costh<-1.) costh=-1.;
    if (costh>cosmx) costh=cosmx;
    sn=1.;
    val=hg->Interpolate(costh,s);
    for (int ns=0; ns<3; ns++) {
      In[ns]+=val*sn;
      sn*=s;
    }
  }
  for (int ns=0; ns<3; ns++) In[ns]*=sbinw;
}

double DoIntegralIso(TH1D *hrho, int ns)
{
  int nsbins=hrho->GetNbinsX();//number of s bins
  double s,sn;
  
  double Itg=0.;
  for (int i=1; i<=nsbins; i++) {
    s=hrho->GetXaxis()->GetBinCenter(i);
    sn=1.;
    for (int j=0; j<ns; j++) {
      sn*=s;
    }
    Itg+=hrho->Interpolate(s)*sn;
  }
  Itg*=hrho->GetXaxis()->GetBinWidth(1);
  return Itg;
}

int main(int argc, char* argv[]){
  
  int nR0bin=401;
  int nth0bin=201;
  double R0max=5000.;
  
  int flgHK=0;
  if (argc>2) {
    flgHK=atoi(argv[2]);
  }
  
  if (flgHK==1) {// for HK!
    nR0bin=1201;
    nth0bin=101;
    R0max=30000.;
  }
  
  cout << "R0max=" << R0max << ", nR0bin=" << nR0bin << ", nth0bin=" << nth0bin << endl;

  int PID;
  int imom,i,j,k,nmom;
  const int nmommax=2000;
  int armom[nmommax];
  Double_t mombEdgs[nmommax+1],R0binEdgs[2001],th0binEdgs[2001];
  double R0,costh0;
  TH2D *wtgtmp;
  TH3F *hI3d[3];
  TFile *fin,*fout;
  
  TH1D *rhotmp,*hI1d[2];
  
  PID=atoi(argv[1]);
  cout << "Particle Code: " << PID << endl;
  
  fin = new TFile(Form("%d_wt.root",PID));
  
  if (fin->IsZombie()) {
    cout << "File does not exist!" << endl;
    delete fin;
    exit(-1);
  }
  
  cout << "Available momenta:" << endl;
  nmom=0;
  for (imom=1; imom<=11000; imom++) {
    wtgtmp = (TH2D*)fin->Get(Form("g_%d_%d",PID,imom));
    if (wtgtmp==NULL) continue;
    cout << imom << endl;
    armom[nmom]=imom;
    mombEdgs[nmom]=imom;
    nmom+=1;
  }
  
  mombEdgs[nmom]=mombEdgs[nmom-1]+1.;
  
  for (i=0; i<=nR0bin; i++) {
    R0binEdgs[i]=i*R0max/(nR0bin-1);
  }
  for (i=0; i<=nth0bin; i++) {
    th0binEdgs[i]=i*2./nth0bin-1.;
  }
  
  for (i=0; i<3; i++) {
    hI3d[i] = new TH3F(Form("hI3d_%d",i),Form("I_{%d}",i),nR0bin,R0binEdgs,nth0bin,th0binEdgs,nmom,mombEdgs);
  }
  
  double Inarr[3];
  cout << "Integrating..." << endl;
  for (imom=0; imom<nmom; imom++) {//loop over all available momentum
//    if (imom>2) break;
    cout << "p=" << armom[imom] << "MeV" << endl;
    wtgtmp = (TH2D*)fin->Get(Form("g_%d_%d",PID,armom[imom]));
    for (j=1; j<=nR0bin; j++) {
      R0=hI3d[0]->GetXaxis()->GetBinLowEdge(j);
      for (k=1; k<=nth0bin; k++) {
        costh0=hI3d[0]->GetYaxis()->GetBinLowEdge(k);
        DoIntegral(wtgtmp,R0,costh0,Inarr);
        for (i=0; i<3; i++) hI3d[i]->SetBinContent(j,k,imom+1,Inarr[i]);//Bin contents are evaluated at bin low edges!
      }
    }
  }
  
  for (i=1; i<3; i++) {
    hI1d[i-1] = new TH1D(Form("hI_iso_%d",i),Form("I^{iso}_{%d}",i),nmom,mombEdgs);
  }
  
  cout << "Integrating isotropic source profile..." << endl;
  for (imom=0; imom<nmom; imom++) {//loop over all available momentum
    //    if (imom>2) break;
    cout << "p=" << armom[imom] << "MeV" << endl;
    rhotmp = (TH1D*)fin->Get(Form("rho_%d_%d",PID,armom[imom]));
    for (i=1; i<3; i++) {
      cout << "n=" << i << endl;
      hI1d[i-1]->SetBinContent(imom+1,DoIntegralIso(rhotmp,i));//Bin contents are evaluated at bin low edges!
    }
  }
  
  TGraph *gNphot = (TGraph*)fin->Get("gNphot");
  TGraph *gS = (TGraph*)fin->Get("gsthr");
  
  fout = new TFile(Form("CProf_%d.root",PID),"RECREATE");
  
  for (i=0; i<3; i++) {
    hI3d[i]->Write();
    delete hI3d[i];
  }
  for (i=0; i<2; i++) {
    hI1d[i]->Write();
    delete hI1d[i];
  }
  
  gNphot->Write();
  gS->Write();
  
  delete fout;
  delete fin;
  return 0;
}
