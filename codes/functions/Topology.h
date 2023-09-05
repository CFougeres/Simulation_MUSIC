#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#define MaxN 5000
#define Npad 1024
#define Nsp 125

string pathTrees  = path + "/tree_results/";
string pathSP  = "/Users/chloefougeres/Documents/StoPow/";
//string pathScattering  = path + "/codes/scattering/";

double mu = 931.5;
double pi_val = TMath::Pi();
//MUSIC geometry
double HalfHeight ;
double DeadLayer ;
double PadDim ;
double DetectorDim[2] ;
double AnodeDim[2] ;
double entranceWindow_width=2.9874527*0.001;
double si_dead_layer_width=800*pow(10,-1)*pow(10,-6);//mm
double dxWindow=0.1;//microns
double dxSi=0.5*pow(10,-4);//microns
double dxGas=5;//microns
double pressure;
double lengthDe;
double EbeamEntry;
double ZAMBeam[3];
double ZAMGas[3];
double ZAMRecoil[MaxN][3];
double ZAMultEjectil[MaxN][3];
double ExRecoil[MaxN];

double m4He = 4 * mu + 2.425;
double m3He = 3 * mu + 14.931;
double m2H = 2 * mu + 13.136;
double mn = mu + 8.071;
double mp = mu + 7.289;
double mBeam;
double mGas;
int Nrecoil;
double mRecoil[MaxN];
double mEjectil[MaxN];
double multi_Ejectil[MaxN];
double mRed[MaxN];
double Qval[MaxN];
Double_t masses[MaxN][2];

Double_t EnergySp[MaxN][4][Nsp];
Double_t SPel[MaxN][4][Nsp];
Double_t SPnu[MaxN][4][Nsp];
string fileSP;
int SPwindow;
int SPgas;
int SPsi;
Int_t PressureCalib[MaxN];
int Npressure;
double Pfais;
TGenPhaseSpace event[MaxN];

//TREE structure
TTree* tree = new TTree("tree", "tree");
Float_t EcmReac;
Float_t errEcmReac;
Int_t XReac;
Int_t YReac;
Int_t ReacKind;
Float_t Epad[Npad][Npad]={0};
int NpadX;
int NpadY;


//Variables
int comptPadX;
int comptPadY;
double vect;
double vect1;
double vect2;
double vect3;
double ReacPad[2];
double EntrancePad[2]={0};
int Nreac;
int NEx;
double positionPad=0.0;
double dx;//mm

double elosswindow;
double elossSi;
double Ebeam ;
double vbeamlab;
double ElossBeam;
double thetaRecoil_lab ;
double vRecoil_lab;
double eRecoil_lab;
double ElossRecoil;
double distRecoil;
double thetaEjectil_lab;
double vEjectil_lab;
double eEjectil_lab;
double ElossEjectil;
double distEjectil;
double EtempLossbeam;
double Etempbeam;
double Extemp;
double Exmax;
double ELossBeamRef[Npad][Npad];
double sigmaPad;
int chosenPad[2];
int ReactionChoice;
int  NReaction;
double theta_scat;
double cs[Npad+2];
double errcs[Npad+2];
double thick_Eeff;
int ThickTargetYield;
double incEloss_thick;
double incEloss_thick2;
double incSP;

#endif
