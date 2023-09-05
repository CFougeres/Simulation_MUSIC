#ifndef MUSICCALIB.H
#define MUSICCALIB.H

#include "main.h"

TCanvas* cEloss = new TCanvas("cEloss","cEloss",800,800);
double Er[2];
double bin_to_E= 0.05;//MeV
int binEr;
TH1F* beamPeak[MaxN];
string nameaxis;
int col;
double FWHM_Si;

TRandom* rESi = new TRandom();

double scattered_theta(int Press) {
    TH1F* srim_scattered;
    TFile* scattered_file;
    string scattered_file_path;
    gRandom->SetSeed(0);
    double theta_p;
    if(Press==0){
        scattered_file_path= pathScattering +Form("%i%iBeamTi_theta_distri.root",int(ZAMBeam[1]),int(ZAMBeam[0]));
        scattered_file = new TFile(scattered_file_path.c_str(), "READ");
        srim_scattered = (TH1F*)scattered_file->Get("theta_distri");
    }
    if (Press>0) {
        scattered_file_path= pathScattering +Form("%i%iBeam%i%iGas%iTorr_theta_distri.root",int(ZAMBeam[1]),int(ZAMBeam[0]), int(ZAMGas[1]),int(ZAMGas[0]), Press);
        scattered_file = new TFile(scattered_file_path.c_str(), "READ");
        srim_scattered = (TH1F*)scattered_file->Get("theta_distri");
    }
    theta_p = srim_scattered->GetRandom();
    scattered_file->Close();
    return theta_p;
}

#endif
