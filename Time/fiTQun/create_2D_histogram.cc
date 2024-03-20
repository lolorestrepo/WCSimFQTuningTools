#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <regex>
#include <map>
#include <experimental/filesystem>

#include <TFile.h>
#include <TH2D.h>
#include <TClonesArray.h>

#include "fiTQun.h"
#include "fiTQun_shared.h"
#include "WCSimWrap.h"
#include "TRuntimeParameters.hxx"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

using namespace std;
namespace fs = std::experimental::filesystem;

const double pi = M_PI;
const double c0 = 29.9792458; // light velocity in vacuum [cm/ns]


int main(int argc, char* argv[]){

  /////////////////////////////////
  //  Read arguments
  /////////////////////////////////
  if (argc<4) {
    cout << "Wrong number of arguments, at least 3 arguments are required:" << endl;
    cout << "infile" << endl;
    cout << " " << endl;
    exit(-1);
  }

  // Read input directory
  char* indirectory = argv[1];
  if (!fs::exists(indirectory)) { cout << "Input directory does not exists" << endl; exit(-1);}
  cout << "Input directory: " << indirectory << endl;

  // Read particle type
  string particle = argv[2];

  // Read parameters filename
  char* parameters_filename = argv[3];
  if (!fs::exists(parameters_filename)) { cout << "Parameters file not found" << endl; exit(-1);}
  cout << "Parameters file: " << parameters_filename << endl;

  // Read water refractive index
  double n_water;
  if (argc>4) { n_water = atof(argv[4]); cout << "Water refractive index set to: " << n_water << endl;}
  else        { n_water = 1.38;          cout << "Water refractive index set to default: " << n_water << endl;}

  // Read vertical axis and define rotation matrix
  int vaxis;
  TRotation R = TRotation();
  if (argc>5) { vaxis = atoi(argv[5]); cout << "Rotating detector vertical axis  " << vaxis << endl;}
  else        { vaxis = 2            ; cout << "Detector vertical axis unrotated "          << endl;}
  if ((vaxis<0) & (2<vaxis)){ cout << "Invalid vertical axis" << endl; exit(-1);}
  if      (vaxis == 0) { R.RotateY(3.*pi/2.);}
  else if (vaxis == 1) { R.RotateX(pi/2.);}

  // Set contained events statistics (min number of contained events to use, if enough are provided)
  int contained_limit = 10000;
  if (argc>6){contained_limit = atoi(argv[6]);}
  cout << "Contained event statistics set to " << contained_limit << endl;

  // Read histogram binning
  if ((7<argc)&(argc<13)){ 
    cout << "Wrong number of arguments for histogram binning, you must set the four arguments: nx xlow xupp ny ylow yupp" << endl;
    exit(-1);
  }
  // Define default values
  int    n_tres    = 200;
  double tres_low  = -3;
  double tres_upp  = 3;
  int    n_trueq   = 125;
  double trueq_low = -2;
  double trueq_upp = 3;
  if (argc == 13){
    n_tres    = atoi(argv[7]);
    tres_low  = atof(argv[8]);
    tres_upp  = atof(argv[9]);
    n_trueq   = atoi(argv[10]);
    trueq_low = atof(argv[11]);
    trueq_upp = atof(argv[12]);
  }

  /////////////////////////////////////////////
  /////////////////////////////////////////////


  // Group filenames in input directory for each momentum
  map<double, vector<fs::path>> grouped_filenames;
  smatch match;
  double energy;
  int index;
  // Define regular expresion to match filenames
  // first note that if particle charge is positive, the + has to be sustituted by \\+
  string particle_pattern = particle;
  size_t pos = particle_pattern.find("+");
  if (pos != std::string::npos) { particle_pattern.replace(pos, 1, "\\+");}
  regex filename_pattern("out_" + particle_pattern + "_(-?\\d+(?:\\.\\d+)?)_(\\d+)\\.root");

  // Loop through files in directory and save them to "grouped_filenames" 
  // such that the files get grouped for each simulated energy
  for (const auto& entry : fs::directory_iterator(indirectory)) {
    if (fs::is_regular_file(entry)) {
      fs::path full_filename = entry.path().string();
      string filename = full_filename.filename();
      // Find the "energy" and file "index"
      if (std::regex_search(filename, match, filename_pattern)){
        energy = std::stod(match[1]);
        index  = std::stoi(match[2]);
        // Fill "grouped_filenames"
        if (grouped_filenames.find(energy) == grouped_filenames.end()) { 
          grouped_filenames[energy] = std::vector<fs::path>{full_filename};
        } else {grouped_filenames[energy].push_back(full_filename);}
      } else {continue;}
    }
  }

  ////////////////////////////////////////
  //  Create fiTQun instance 
  ////////////////////////////////////////
  TRuntimeParameters::Get().ReadParamOverrideFile(parameters_filename);

  // Load WCSim wrapper and get geometry (needed later)
  char* filename = const_cast<char*> ((grouped_filenames.begin()->second[0]).c_str());

  WCSimWrap    * wc       = WCSimWrap::Get(filename);
  WCSimRootGeom* geometry = wc->Geo();
  // Read first event:
  wc->LoadEntry(0);

  // Search the track of the parent particle (parenttype = 0, id = 1)
  // The search starts at second track (first 2 tracks special, they correspond to beam and target which in this case are dummy)
  WCSimRootTrigger* trigger = wc->SubEvt(0); // only trigger 0 is considered
  TClonesArray    * tracks  = trigger->GetTracks();
  WCSimRootTrack* track;
  for (int i=2; i<tracks->GetEntries(); i++){
    track = (WCSimRootTrack*) tracks->At(i);
    if ((track->GetParenttype()==0)&(track->GetId()==1)){break;}
  }
  if ((!(track->GetParenttype()==0)&(track->GetId()==1))){
    cout << "Not found track of parent particle for the first event" << endl;
    exit(-1);
  }

  // define load_flag and PID_index
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

  // Finally create instance and load parameters
  fiTQun       * fq       = new fiTQun(wc->NPMT());
  fiTQun_shared* fqshared = fiTQun_shared::Get(wc->NPMT());

  fq->ReadSharedParams(1,load_flag,true,true,1);
  // fqshared->SetPhi0(-1.,1.); // ?

  wc->Close();

  ////////////////////////////////////////
  //  Read WCSim trigger_offset from WCSim options
  ////////////////////////////////////////
  WCSimRootOptions* op = new WCSimRootOptions();
  TFile* f = new TFile(filename);
  TTree* optionsT = (TTree*) f->Get("wcsimRootOptionsT");
  optionsT->SetBranchAddress("wcsimrootoptions",&op);
  optionsT->GetEntry(0);
  double trigger_offset = op->GetTriggerOffset();
  f->Close();

  ////////////////////////////////////////
  //  Create output file, to be filled
  ////////////////////////////////////////
  const char* ofilename = "tres_trueq_2Dhistogram.root";
  TFile* of = new TFile(ofilename, "RECREATE");
  of->Close();
  
  ////////////////////////////////////////
  // Loop through each energy value
  ////////////////////////////////////////

  // Save momenta, total processed and PC events
  vector<double> momenta;
  vector<double> ntotal;
  vector<double> contained;
  // Save conversion between nphotons and total charge
  vector<double> nphotons;
  vector<double> totalqs;

  for (auto& entry : grouped_filenames) {

    cout << "--------------------" << endl;
    cout << "Processing E = " << entry.first << endl;
    cout << "--------------------" << endl;
    sort(entry.second.begin(), entry.second.end());

    // Read first file in group and extract the following information:
    //    - trigger offset
    //    - particle momentum
    //    - maximum distance travelled for this particle with this momentum
    //    - instantiate fiTQun and fitqun_shared

    char* filename = const_cast<char*> (entry.second[0].c_str());

    // Load WCSim wrapper
    WCSimWrap* wc = WCSimWrap::Get(filename);

    // Read first event:
    wc->LoadEntry(0);

    // Search the track of the parent particle (parenttype = 0, id = 1)
    // The search starts at second track (first 2 tracks special, they correspond to beam and target which in this case are dummy)
    WCSimRootTrigger* trigger = wc->SubEvt(0); // only trigger 0 is considered
    TClonesArray    * tracks  = trigger->GetTracks();
    WCSimRootTrack  * track;
    for (int i=2; i<tracks->GetEntries(); i++){
      track = (WCSimRootTrack*) tracks->At(i);
      if ((track->GetParenttype()==0)&(track->GetId()==1)){break;}
    }
    if ((!(track->GetParenttype()==0)&(track->GetId()==1))){
      cout << "Not found track of parent particle for the first event" << endl;
      exit(-1);
    }

    // Get momentum
    double momentum = track->GetP();
    momenta.push_back(momentum);

    // Compute the smax value for this particle and momentum
    double Iiso[2], nphot, smax;
    fiTQun_shared::Get()->SetTypemom(PID_index, momentum, Iiso, nphot, smax);
    nphotons.push_back(nphot);

    wc->Close();
    
    // Define the 2D time pdf histograms to be filled
    string hname_direct = "htimepdf_direct_" + to_string(momentum);
    TH2D* hdirect = new TH2D(hname_direct.c_str(),"", n_tres, tres_low, tres_upp, n_trueq, trueq_low, trueq_upp);
    hdirect->Sumw2();

    string hname_indirect = "htimepdf_indirect_" + to_string(momentum);
    TH2D* hindirect = new TH2D(hname_indirect.c_str(),"", n_tres, tres_low, tres_upp, n_trueq, trueq_low, trueq_upp);
    hindirect->Sumw2();

    //////////////////////////////////////////
    // Loop through each file
    //////////////////////////////////////////
    cout << "Looping through files" << endl;
    int ntotal_counter    = 0;
    int contained_counter = 0;
    double totalq         = 0;
    for (const fs::path& path : entry.second) {

      char* filename = const_cast<char*> (path.c_str());
      cout << "    =>" << filename << endl;

      // Read file
      WCSimWrap* wc = WCSimWrap::Get(filename);

      // Loop through each event
      ntotal_counter += wc->NEvt();
      for (int event = 0; event < wc->NEvt(); event++){
        
        // cout << "Processing event " << (event+1) << "/" << wc->NEvt() << endl;
        wc->LoadEntry(event);

        // Repeat the code above to find the parent track, if parent not found, continue processing next event
        trigger = wc->SubEvt(0);
        tracks  = trigger->GetTracks();
        for (int i=2; i<tracks->GetEntries(); i++){
          track = (WCSimRootTrack*) tracks->At(i);
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
        
        // Get the predicted charge at each pmt and the PC flag
        // Direct
        fqshared->SetScatflg(0);
        double true_q_direct[nPMT_max];
        int isPC = fq->Get1Rmudist(PID_index, track_params, true_q_direct);
        // Total
        fqshared->SetScatflg(1);
        double true_q_total[nPMT_max];
        fq->Get1Rmudist(PID_index, track_params, true_q_total);

        // If the event is PC, isPC = 1
        if (isPC == 0) {
          contained_counter += 1;
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
            TVector3 midpos = vertex + direction * (smax/2.);
            double midtrack_pmt_distance = (pmt_pos - midpos).Mag();
            double t = trigger_offset - trigger->GetHeader()->GetDate();
            double t_residual = digi_hit->GetT() - t - (smax/2.)/c0 - midtrack_pmt_distance*n_water/c0;

            hdirect  ->Fill(t_residual, log10(true_q_direct[pmtid]));
            hindirect->Fill(t_residual, log10(true_q_total [pmtid] - true_q_direct[pmtid]));

            totalq += true_q_total [pmtid];
          }

        }
      }

      wc->Close();

      // break if statistics enough
      cout << "# Processed events: " << contained_counter << endl;
      if (contained_counter >= contained_limit){break;}
    }

    // Append number of total and pc events, and mean charge per event
    ntotal   .push_back(ntotal_counter);
    contained.push_back(contained_counter);
    totalqs  .push_back(totalq/contained_counter);

    // Write timepdf to output file
    cout << "Writing to " << ofilename << std::endl;
    TFile* of = new TFile(ofilename, "UPDATE");
    hdirect    ->Write();
    hindirect  ->Write();
    of->Close();
    cout << " " << endl;

    delete hdirect;
    delete hindirect;
  }

  // Define and write graphs for total and PC events, and conversion nphotons-charge
  TGraph* gntotal    = new TGraph( momenta.size(),  momenta.data(), ntotal   .data());
  TGraph* gcontained = new TGraph( momenta.size(),  momenta.data(), contained.data());
  TGraph* gconv      = new TGraph(nphotons.size(), nphotons.data(),   totalqs.data());
  gntotal   ->SetName("Total");
  gcontained->SetName("Contained");
  gconv     ->SetName("Conversion");
  of = new TFile(ofilename, "UPDATE");
  gntotal   ->Write();
  gcontained->Write();
  gconv     ->Write();
  of->Close();

  return 0;
}
