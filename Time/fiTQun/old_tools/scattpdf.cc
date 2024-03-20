#include <iostream>
#include "fiTQun.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <cmath>

#include "skheadC.h"
#include "skparmC.h"
#include "geopmtC.h"
#include "sktqC.h"
#include "vcworkC.h"
#include "vcvrtxC.h"

extern"C" {
  void fortinit_(char* fname_in, int nfname_in);
  int fortread_();
  void fortconsts_();
  void trginfo_(float *);
  void vcrdvccm_();
  
  extern struct{
    int isk23[MAXPM];
    float qefactor[MAXPM];
  } pmtinfcmn_;
}

const double pi=3.1415926535898;
const double c0=29.9792458;//[cm/ns]
double nwtr=1.36;

int main(int argc, char* argv[]){
  
  bool fcall=false;
  int PCflg=0;
  TString ntplFile;
  char fstr[256];
  
  double phi0;
  const int PIDarr[4]={11,13,211,2212};
  
  if (argc<2) exit(-1);
  
  if (argc==3) {//set water refractive index
    nwtr=atof(argv[2]);
    std::cout << "nwtr=" << nwtr << std::endl;
  }
  
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
  ntplFile=fstr;
  
  fortconsts_();
  
  int PMTtype[nPMT];
  double PMTQE[nPMT];
  for (int icab=0; icab<nPMT; icab++) {
    PMTtype[icab]=pmtinfcmn_.isk23[icab]-2;//set PMT type(0:Old, 1:New)
    if (PMTtype[icab]!=0 && PMTtype[icab]!=1) PMTtype[icab]=-1;//Dead PMT
    PMTQE[icab]=pmtinfcmn_.qefactor[icab];
  }
  
  fortinit_(argv[1],nifile_name);
  
  int iPID;
  int nevt=0;
  int icab;
  double TrkParam[7],muarr[nPMT],musarr[nPMT],momtmp[3],trkdir[3],tc,midpos[3],smid,mom,RmidPMT;
  float trgofst;
  int iret = 0;
  
  TFile *fcprof;
  TGraph *gSmax;
  TString fiTQundir;
  
  
  char *ctmp = getenv("FITQUN_ROOT");
  if (ctmp!=NULL) {
    fiTQundir=ctmp;
    fiTQundir+="/";
  }
  else {
    ctmp = getenv("ATMPD_ROOT");
    if (ctmp!=NULL) {
      fiTQundir=ctmp;
      fiTQundir+="/src/recon/fitqun/";
    }
    else {
      std::cout << "ATMPD_ROOT not set!" << std::endl;
      exit(-1);
    }
  }
  
  fiTQun thefit,scatfit;
  
  thefit.SetPhi0(0,1.);//remove dependence on tuning const.
  thefit.SetScatflg(0);//no scattered light predicted charge
  
  scatfit.SetPhi0(0,1.);//remove dependence on tuning const.
  
  TH2D htimepdf("htimepdf","",100,-25,175,100,-2.,2.);
  while (1) {//event loop
    
    iret = fortread_();
    
    if (iret>0) break;
    
    nevt++;
    
    vcrdvccm_();
    
    TrkParam[3]=vcwork_.timvc[0];//time
    for (int j=0; j<3; j++) {
      TrkParam[j]=vcwork_.posivc[0][j];//vertex
      momtmp[j]=vcwork_.pvc[0][j];//momentum (MeV/c)
    }
    TrkParam[6]=sqrt(momtmp[0]*momtmp[0]+momtmp[1]*momtmp[1]+momtmp[2]*momtmp[2]);
    for (int j=0; j<3; j++) trkdir[j]=momtmp[j]/TrkParam[6];
    TrkParam[4]=acos(trkdir[2]);//theta
    TrkParam[5]=atan2(trkdir[1],trkdir[0]);//phi
    
    if (!fcall) {//identify particle type
      fcall=true;
      for (iPID=0; iPID<4; iPID++) {
        if (PIDarr[iPID]==vcwork_.ipvc[0]) break;
      }
      if (iPID>=4) exit(-1);
      std::cout << "Particle code: " << PIDarr[iPID] << std::endl;
      fcprof = new TFile(Form(fiTQundir+"CProf_%d.root",PIDarr[iPID]));
      gSmax = new TGraph(*(TGraph*)fcprof->Get("gsthr"));
      mom=TrkParam[6];
      smid=gSmax->Eval(mom)/2.;
      std::cout << "p=" << mom << "MeV/c, smid: " << smid << std::endl;
      delete gSmax;
      delete fcprof;
    }
    
    for (int j=0; j<3; j++) midpos[j]=TrkParam[j]+smid*trkdir[j];
    
    PCflg=thefit.Get1Rmudist(iPID,TrkParam,muarr);
    PCflg=scatfit.Get1Rmudist(iPID,TrkParam,musarr);
    
    for (int i=0; i<11146; i++) {
      musarr[i]-=muarr[i];
    }
    
    std::cout << "Event#:" << nevt << std::endl;
//    for (int j=0; j<4; j++) std::cout << TrkParam[j] << ", ";
//    std::cout << std::endl << TrkParam[6] << ": ";
//    std::cout << trkdir[0] << ", " << trkdir[1] << ", " << trkdir[2] << std::endl;
    
    
    trginfo_(&trgofst);
//    std::cout << trgofst << std::endl;
    trgofst+=1000.;
    
    if (PCflg==0) {//FC event
      for (int ih=0; ih<skq_.nqisk; ih++) {//loop over hit PMTs
        icab = skchnl_.ihcab[ih] - 1;//cable # in c++
        
        if (musarr[icab]<10.*muarr[icab]) continue;
        
        RmidPMT=sqrt((geopmt_.xyzpm[icab][0]-midpos[0])*(geopmt_.xyzpm[icab][0]-midpos[0])
                     +(geopmt_.xyzpm[icab][1]-midpos[1])*(geopmt_.xyzpm[icab][1]-midpos[1])
                     +(geopmt_.xyzpm[icab][2]-midpos[2])*(geopmt_.xyzpm[icab][2]-midpos[2]));
        
        tc = skt_.tisk[icab] - trgofst - RmidPMT*nwtr/c0 - smid/c0;//corrected time
//        std::cout << icab << ", "  << tc << ", " << muarr[icab] << std::endl;
        htimepdf.Fill(tc,log10(musarr[icab]));
      }
    }
    else {
      std::cout << "PC event!!" << std::endl;
    }

//    if (nevt>30) break;
  }
  
  TFile *ofile = new TFile(ntplFile + "_hist.root", "RECREATE");
  htimepdf.Write();
  delete ofile;

}
