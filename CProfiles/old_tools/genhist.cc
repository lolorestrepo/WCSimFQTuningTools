//script to normalize the cherenkov profiles and put them into one file

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TFile.h>

using namespace std;


int genhist(int PID){//PDG

    /////// MOVE TO ARGS ?
    const char typ[]="wt"; //Track direction: "wt":weighted, "tr":true

    // const char filename[100]="";
    // strcat(filename, getenv("LUSTRE"));
    // strcat(filename, "/CProfiles/out/mu-/");
    // strcat(filename, "cprofile_%dMeV_%d_ccin2p3.root");
    ///////////////
    string fname = getenv("LUSTRE");
    fname += "/CProfiles/full_e-/out/";
    // fname += "cprofile_%d.0MeV_%d.root";
    fname += "out_e-_%d_%d.root";
    const char* filename = fname.c_str();


    /////// DECLARATIONS
    TFile *fin,*fout;
    int imom,j,nmom,nevt;
    double Itmp;
    const int nmommax=2000;
    Double_t armom[nmommax],arnphot[nmommax],thrs[nmommax];
    TH1D *hprj[nmommax];
    TH1D *hnevt;
    TH2D *wtg[nmommax],*wtgtmp;
    ///////////////

    char cPath[256];
    strcpy(cPath,gDirectory->GetPath());

    nmom=0;
    for (imom=0; imom<=5000; imom+=10) {

        if (wtg[nmom]){ 
            // delete wtg[nmom]; //????
            wtg[nmom] = nullptr;
        }

        cout << "Mom: " << imom << endl;
        
        nevt=0;
        
        // open files 1 and 2 file to check if files for this momentum exist
        fin = new TFile(Form(filename,imom,1));
        if (fin->IsZombie()) {
            delete fin;
            fin = new TFile(Form(filename,imom,2));
            if (fin->IsZombie()) {
                delete fin;
                continue;
            }
        }

        for (int j=1; j<100; j++) {

            cout << Form(filename,imom,j) << endl;

            // ifstream file(Form(filename,imom,j));
            // if (!file.good()){break;}

            fin = new TFile(Form(filename,imom,j));
            gDirectory->cd(cPath);
            if (fin->IsZombie()) {
                delete fin;
                continue;
            }
            else if (!wtg[nmom]) {
                wtg[nmom] = new TH2D(*(TH2D*)fin->Get(Form("%sg",typ)));
                wtg[nmom]->Reset();
            }

            hnevt=(TH1D*)fin->Get("nevents");
            nevt+=hnevt->GetEntries();
            wtgtmp=(TH2D*)fin->Get(Form("%sg",typ));
            wtg[nmom]->Add(wtgtmp);
            delete fin;
        }
        if (!wtg[nmom]) continue;//skip this momentum

        armom[nmom]=imom;
        arnphot[nmom]=wtg[nmom]->GetEntries()/nevt;

        // cout << std::fixed << std::setprecision(10) << wtg[nmom]->GetEntries() << " " << nevt << " " << arnphot[nmom] << endl;

        wtg[nmom]->Scale(1./wtg[nmom]->GetEntries(),"width");
        wtg[nmom]->SetName(Form("g_%d_%d",PID,imom));

        //1d cherenkov emission profile
        hprj[nmom] = wtg[nmom]->ProjectionY(Form("rho_%d_%d",PID,imom));
        hprj[nmom]->Scale(wtg[nmom]->GetXaxis()->GetBinWidth(1));
        Itmp=0.;
        for (int i=1; i<=hprj[nmom]->GetNbinsX(); i++) {
            Itmp+=hprj[nmom]->GetBinContent(i)*hprj[nmom]->GetXaxis()->GetBinWidth(i);
            if (Itmp>0.9) {
                thrs[nmom]=hprj[nmom]->GetXaxis()->GetBinCenter(i);
                break;
            }
        }

        nmom+=1;
        if (nmom>=nmommax) {
            cout << "Too many momentum bins!" << endl;
            exit(-1);
        }
    }

    TGraph *gNphot = new TGraph(nmom,armom,arnphot);
    gNphot->SetName("gNphot");
    gNphot->SetTitle("Mean number of photons per event");

    TGraph *gS = new TGraph(nmom,armom,thrs);
    gS->SetName("gsthr");
    gS->SetTitle("Length of track");

    fout = new TFile(Form("%d_%s.root",PID,typ),"RECREATE");
    for (j=0; j<nmom; j++) {
        wtg[j]->Write();
        hprj[j]->Write();
    }

    gNphot->Write();
    gS->Write();

    delete fout;

    return 0;
}
