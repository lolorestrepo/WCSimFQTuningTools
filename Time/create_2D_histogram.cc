#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <sys/stat.h>

#include <TFile.h>
#include <TH2D.h>
#include <TClonesArray.h>

// Load them instad of including?
#include "fiTQun.h"
#include "fiTQun_shared.h"
#include "WCSimWrap.h"
#include "TRuntimeParameters.hxx"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

using namespace std;

const double pi = M_PI;
const double c0 = 29.9792458; // light velocity in vacuum [cm/ns]
// double n_water = 1.38;        // water refractive index       (argument?)


bool fileExists(const string& filePath) {
    struct stat buffer;
    return (stat(filePath.c_str(), &buffer) == 0);
}

int main(int argc, char* argv[]){

  //////////////////// Read arguments
  /////////////////////////////////////////////
  if (argc<3) {
    cout << "Wrong number of arguments, at least 2 arguments are required:" << endl;
    cout << "infile" << endl;
    cout << " " << endl;
    exit(-1);
  }

  cout << argc << endl;

  // Read input filename
  char* infilename = argv[1];
  if (!fileExists(infilename)) { cout << "Input file not found" << endl; exit(-1);}
  cout << "Input file: " << argv[1] << ";" << infilename << endl;

  // Read parameters filename
  char* parameters_filename = argv[2];
  if (!fileExists(parameters_filename)) { cout << "Parameters file not found" << endl; exit(-1);}
  cout << "Parameters file: " << argv[2] << ";" << parameters_filename << endl;

  // Read water refractive index
  double n_water;
  if (argc>3) { n_water = atof(argv[3]); cout << "Water refractive index set to: " << n_water << endl;}
  else        { n_water = 1.38;          cout << "Water refractive index set to default: " << n_water << endl;}

  // Read vertical axis
  int vaxis;
  if (argc>4) { vaxis = atoi(argv[4]); cout << "Rotating detector vertical axis  " << vaxis << endl;}
  else        { vaxis = 2            ; cout << "Detector vertical axis unrotated "          << endl;}
  if ((vaxis<0) & (2<vaxis)){ cout << "Invalid vertical axis" << endl; exit(-1);}

  // define rotation matrix based on vaxis
  TRotation R = TRotation();
  if      (vaxis == 0) { R.RotateY(3.*pi/2.);}
  else if (vaxis == 1) { R.RotateX(pi/2.);}
  
  /////////////////////////////////////////////
  /////////////////////////////////////////////


  // Read WCSim geometry
  WCSimWrap    * wc       = WCSimWrap::Get(argv[1]);
  WCSimRootGeom* geometry = wc->Geo();
  // Read WCSim trigger_offset from WCSim options
  WCSimRootOptions* op = new WCSimRootOptions();
  TFile* f = new TFile(argv[1]);
  TTree* optionsT = (TTree*) f->Get("wcsimRootOptionsT");
  optionsT->SetBranchAddress("wcsimrootoptions",&op);
  optionsT->GetEntry(0);
  double trigger_offset = op->GetTriggerOffset();
  f->Close();

  // Read first event, which will be used to:
  //    - define fiTQun load flag 
  //    - obtain max value of the distance travelled by a particle with its momentum
  wc->LoadEntry(0);

  // Search the track of the parent particle (parenttype = 0, id = 1)
  // The search starts at second track (first 2 tracks special, they correspond to beam and target which in this case are dummy)
  WCSimRootTrigger* trigger = wc->SubEvt(0); // only trigger 0 is considered
  TClonesArray    * tracks  = trigger->GetTracks();
  WCSimRootTrack* track;
  int i = 2;
  while(true){
    track = (WCSimRootTrack*) tracks->At(i);
    i++;
    if ((track->GetParenttype()==0)&(track->GetId()==1)){break;}
  }
  if ((!(track->GetParenttype()==0)&(track->GetId()==1))){
    cout << "Not found track of parent particle for the first event" << endl;
    exit(-1);
  }

  // Find particle index and define load_flag required by fiTQun 
  int PID_index;
  int load_flag[nPID] = {0};
  auto found = find(begin(fiTQun_shared::PIDarr), end(fiTQun_shared::PIDarr), track->GetIpnu());
  if (found != std::end(fiTQun_shared::PIDarr)) {
        PID_index = std::distance(std::begin(fiTQun_shared::PIDarr), found);
        load_flag[PID_index] = 1;
  }
  else {
    std::cout << "Unknown particle type: " << track->GetIpnu() << std::endl;
    exit(-1);
  }
  std::cout << "Particle code: " << fiTQun_shared::PIDarr[PID_index] << std::endl;
  
  // Get a fiTQun instance
  TRuntimeParameters::Get().ReadParamOverrideFile(argv[2]); // (needed here? if not move just above the fiTQun instantiation)
  fiTQun       * fq = new fiTQun(wc->NPMT());
  fiTQun_shared* fqshared = fiTQun_shared::Get(wc->NPMT());

  fq->ReadSharedParams(0,load_flag,true,true,1);

  fqshared->SetScatflg(0);   // don't compute predicted charge by scattered light  
  fqshared->SetWAttL(6800.); // ?
  fqshared->SetQEEff(0.1);   // ?
  fqshared->SetPhi0(-1.,1.); // ?

  // Get momentum
  double momentum = track->GetP();

  // Compute the smid value for this particle and momentum
  double Iiso[2], nphot, smax, smid;
  fiTQun_shared::Get()->SetTypemom(PID_index, momentum, Iiso, nphot, smax);
  smid=smax/2.;
  std::cout << "p=" << momentum << "MeV/c, smid: " << smid << std::endl;
  
  // Define the htimepdf 2D histogram to be filled
  TH2D htimepdf("htimepdf","", 400, -100, 100, 125, -2., 3.);
  htimepdf.Sumw2();

  // Loop through the events in the file
  for (int event = 0; event < wc->NEvt(); event++){

    cout << "Processing event " << (event+1) << "/" << wc->NEvt() << endl;
    wc->LoadEntry(event);
    
    // Repeat the code above to find the parent track, if parent not found, continue processing next event
    trigger = wc->SubEvt(0);
    tracks  = trigger->GetTracks();
    int i = 2;
    while(true){
      track = (WCSimRootTrack*) tracks->At(i);
      i++;
      if ((track->GetParenttype()==0)&(track->GetId()==1)){break;}
    }
    if ((!(track->GetParenttype()==0)&(track->GetId()==1))){continue;}

    // Get parent track parameters (vertex, angles and momentum)
    // Get initial position (vertex) and direction
    double track_params[7];
    TVector3 vertex   (trigger->GetVtx(0), trigger->GetVtx(1), trigger->GetVtx(2));
    TVector3 direction(track  ->GetDir(0), track  ->GetDir(1), track  ->GetDir(2));
    // Rotate vertex and direction if needed
    if (vaxis != 2){
      vertex    = R * vertex;
      direction = R * direction;
    }
    // define track_params as required by fiTQun
    track_params[0] = vertex.X();
    track_params[1] = vertex.Y();
    track_params[2] = vertex.Z();
    track_params[3] = track->GetTime();
    track_params[4] =  acos(direction.Z());
    track_params[5] = atan2(direction.Y(), direction.X());
    track_params[6] = track->GetP();

    for (int a = 0; a <=6; a++){ cout << "track_params " << a << " " << track_params[a] << endl;}
    
    // Calculate the midpoint of the track
    TVector3 midpos = vertex + smid * direction;
    
    // Get the predicted charge at each pmt and the PC flag
    double true_q[nPMT_max];
    int isPC = fq->Get1Rmudist(PID_index, track_params, true_q);

    cout << isPC << endl;

    if (isPC == 0) {
      // Loop through each hit 
      for (int hit = 0; hit < trigger->GetNcherenkovdigihits(); hit++){
        // Get hit
        WCSimRootCherenkovDigiHit* digi_hit = dynamic_cast<WCSimRootCherenkovDigiHit*>(trigger->GetCherenkovDigiHits()->At(hit));

        // Get PMT position and rotate it if needed
        int pmtid = digi_hit->GetTubeId()-1;
        WCSimRootPMT pmt = geometry->GetPMT(pmtid);
        TVector3 pmt_pos(pmt.GetPosition(0), pmt.GetPosition(1), pmt.GetPosition(2));
        if (vaxis != 2){ pmt_pos = R * pmt_pos;}

        // Compute residual time and fill histogram
        double midtrack_pmt_distance = (pmt_pos - midpos).Mag();
        double t = trigger_offset - trigger->GetHeader()->GetDate();
        double t_residual = digi_hit->GetT() - t - smid/c0 - midtrack_pmt_distance*n_water/c0;

        cout << "tres " << t_residual << " logmu " << log10(true_q[pmtid]) << endl;

        htimepdf.Fill(t_residual, log10(true_q[pmtid]));
      }
    }
  }

  const char* ofilename = "tres_trueq_2Dhistogram.root";
  cout << "Writing output to " << ofilename << std::endl;
  TFile* of = new TFile(ofilename, "RECREATE");
  htimepdf.Write();
  of->Close();

}
  
