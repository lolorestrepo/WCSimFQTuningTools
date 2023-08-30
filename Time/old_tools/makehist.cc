#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <cmath>
#include <cstdlib>

#include "skheadC.h"
#include "skparmC.h"
#include "geopmtC.h"
#include "sktqC.h"
#include "vcworkC.h"
#include "vcvrtxC.h"

#include "fiTQun.h"
#include "fiTQun_shared.h"

using namespace std;

extern"C" {
  void fortinit_(char*, char*, int, int);
  int fortread_();
  void trginfo_(float *);
  void vcrdvccm_();
}

const double pi=3.1415926535898;
const double c0=29.9792458;//[cm/ns]
double nwtr=1.38;// Water refractive index for TOF subtraction

int main(int argc, char* argv[]){
  
  int ppt = (getenv("PPTRACKING") == NULL) ? 0 : atoi(getenv("PPTRACKING"));
  
  if (argc<2) exit(-1);
  
  if (argc>=3) {//set water refractive index
    nwtr=atof(argv[2]);
  }
  std::cout << "nwtr=" << nwtr << std::endl;
  
  int fPCflcgut=1;
  if (argc>=4) {
    fPCflcgut=atoi(argv[3]);
  }
  std::cout << "PCflg cut: " << fPCflcgut << std::endl;
  
  char fstr[256];

  int nifile_name;
  for (int i=0; i<256; i++) {
    fstr[i]=argv[1][i];
    if (argv[1][i]=='\0') {
      nifile_name=i;
      break;
    }
  }
  std::cout << "Input file: " << argv[1] << ";" << nifile_name << std::endl;
  
  for (int i=nifile_name; i>0; i--) {
    if (fstr[i]=='.') {
      fstr[i]='\0';
      break;
    }
  }
  TString ntplFile = fstr;
  std::cout << "fstr: " << fstr << std::endl;
  
  TTree * pidTree;
  Float_t xPositionFirstInteraction, yPositionFirstInteraction,
    zPositionFirstInteraction, timeOfFirstInteraction;

  std::cout << "ppt: " << ppt << endl;
  
  if (ppt == 1) {
    TString rootFileName = Form("%s.root",fstr);
    TFile * rootFile = new TFile(rootFileName.Data());
    std::cout << "rootFileName: " << rootFileName << std::endl;
    pidTree = (TTree *) rootFile->Get("pid");
    pidTree->SetBranchAddress("xPosFirstInteraction",&xPositionFirstInteraction);
    pidTree->SetBranchAddress("yPosFirstInteraction",&yPositionFirstInteraction);
    pidTree->SetBranchAddress("zPosFirstInteraction",&zPositionFirstInteraction);
    pidTree->SetBranchAddress("timeOfFirstInteraction",&timeOfFirstInteraction);
  }
  
  
  char* dumstr;
  int dumint=0;
  fortinit_(argv[1],dumstr,nifile_name,dumint);// First event is read!
  
  int iret = 0;
  vcrdvccm_();
  
  int iPID;
  int fLoadflg[nPID];
  for (iPID=0; iPID<nPID; iPID++) fLoadflg[iPID]=0;
  for (iPID=0; iPID<nPID; iPID++) {
    if (fiTQun_shared::PIDarr[iPID]==vcwork_.ipvc[0]) {
      fLoadflg[iPID]=1;
      break;
    }
  }
  if (iPID>=nPID) {
    std::cout << "Unknown particle type: " << vcwork_.ipvc[0] << std::endl;
    exit(-1);
  }
  std::cout << "Particle code: " << fiTQun_shared::PIDarr[iPID] << std::endl;
  
  
  fiTQun *thefit = new fiTQun();
  
  int iScat=0;//no scattered light predicted charge
  
  thefit->SetScatflg(iScat);
  
  thefit->ReadSharedParams(iScat,fLoadflg,true,true,1);
  
  thefit->SetQEEff(0.1);
  thefit->SetPhi0(-1,1.);//remove dependence on tuning const.
  
  
  double mom=0.;
  for (int j=0; j<3; j++) mom+=vcwork_.pvc[0][j]*vcwork_.pvc[0][j];
  mom=sqrt(mom);// True momentum
  
  bool fCG=false;
  if (fiTQun_shared::PIDarr[iPID]==48) {
    fCG=true;
    mom/=10.;
  }
  
  double Iiso[2],nphot,smax;
  fiTQun_shared::Get()->SetTypemom(iPID,mom,Iiso,nphot,smax);
  
  double smid=smax/2.;
  std::cout << "p=" << mom << "MeV/c, smid: " << smid << std::endl;
  
  TH2D htimepdf("htimepdf","",400,-100,100,100,-2.,2.);
  
  int nevt=0;
  while (1) {//event loop
    
    if (iret>0) break;
    
    std::cout << "Event#:" << nevt << std::endl;
    
    vcrdvccm_();
    
    float trgofst;
    trginfo_(&trgofst);
//    std::cout << trgofst << std::endl;
    
    double TrkParam[7];
    TrkParam[3]= vcwork_.timvc[0];//time
//    std::cout << TrkParam[3] << std::endl;
    if (ppt == 1) {
      pidTree->GetEntry(nevt);
      TrkParam[0] = xPositionFirstInteraction;
      TrkParam[1] = yPositionFirstInteraction;
      TrkParam[2] = zPositionFirstInteraction;
      cout << "PPT : TrkParam[0-2] = " << TrkParam[0] << ", " << TrkParam[1]
           << ", " << TrkParam[2] << ", nevt = " << nevt << endl << endl;
    }
    
    double momtmp[3];
    for (int j=0; j<3; j++) {
      if (ppt != 1) {
        TrkParam[j]=vcwork_.posivc[0][j];//vertex
      }
      momtmp[j]=vcwork_.pvc[0][j];//momentum (MeV/c)
    }
    
    TrkParam[6]=sqrt(momtmp[0]*momtmp[0]+momtmp[1]*momtmp[1]+momtmp[2]*momtmp[2]);
    double trkdir[3];
    for (int j=0; j<3; j++) trkdir[j]=momtmp[j]/TrkParam[6];
    TrkParam[4]=acos(trkdir[2]);//theta
    TrkParam[5]=atan2(trkdir[1],trkdir[0]);//phi
    
    if (fCG) TrkParam[6]/=10.;
    
    double midpos[3];
    for (int j=0; j<3; j++) midpos[j]=TrkParam[j]+smid*trkdir[j];
    
    double muarr[nPMT_max];
    int PCflg=thefit->Get1Rmudist(iPID,TrkParam,muarr);
    
    if (PCflg==0 || fPCflcgut==0) {//FC event
      for (int ih=0; ih<skq_.nqisk; ih++) {//loop over hit PMTs
        int icab = skchnl_.ihcab[ih] - 1;//cable # in c++
        
        double RmidPMT=sqrt((geopmt_.xyzpm[icab][0]-midpos[0])*(geopmt_.xyzpm[icab][0]-midpos[0])
                     +(geopmt_.xyzpm[icab][1]-midpos[1])*(geopmt_.xyzpm[icab][1]-midpos[1])
                     +(geopmt_.xyzpm[icab][2]-midpos[2])*(geopmt_.xyzpm[icab][2]-midpos[2]));
        
        double tc;// corrected time
        if (ppt == 1) {
          tc = skt_.tisk[icab] - trgofst - (timeOfFirstInteraction * 1.0e9) - RmidPMT*nwtr/c0 - smid/c0;//corrected time
        } else {
          tc = skt_.tisk[icab] - trgofst - RmidPMT*nwtr/c0 - smid/c0;//corrected time
        }
//        std::cout << icab << ", "  << tc << ", " << muarr[icab] << std::endl;
        htimepdf.Fill(tc,log10(muarr[icab]));
      }
    }
    else {
      std::cout << "PC event!!" << std::endl;
    }
    
    iret = fortread_();
    
    nevt++;
  
  }

  TFile *ofile = new TFile(ntplFile + "_hist.root", "RECREATE");
  htimepdf.Write();
  delete ofile;

}
