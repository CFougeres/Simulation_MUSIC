#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#define MaxN 5000
#define Nstrip 18
#define Nsp 125

string pathTrees  = path + "/tree_results/";
//string pathSP  = path + "/codes/sp/";
string pathSP  = "/Users/chloefougeres/Documents/StoPow/";
string pathScattering  = path + "/codes/scattering/";

double mu = 931.5;

//MUSIC geometry
double MUSIC_half_height =100.;//mm
double MUSIC_dead_layer =35.9;//mm
double strip_dim = 15.78;// mm
double entranceWindow_width=2.9874527*0.001;
double si_dead_layer_width=800*pow(10,-1)*pow(10,-6);//mm
double MUSIC_length=355.6;//mm
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

Double_t EnergySp[MaxN][2][Nsp];
Double_t SPel[MaxN][2][Nsp];
Double_t SPnu[MaxN][2][Nsp];
string fileSP;
int SPwindow;
int SPgas;
int SPsi;
Int_t PressureCalib[MaxN];
int Npressure;
double Pfais;
TGenPhaseSpace event[MaxN];

//Ecm conversion
Double_t ECMSP[Nstrip];
Double_t AnCM[Nstrip];
Double_t ECM[Nstrip];
Double_t erInfAnCM[Nstrip];
Double_t erSupAnCM[Nstrip];
Double_t erInfECM[Nstrip];
Double_t erSupECM[Nstrip];
TTree* tree = new TTree("tree", "tree");
Int_t stripReac ;
Int_t ReacKind;
Float_t Estrip[Nstrip];
Int_t Istrip[Nstrip];

double strip[18]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};


    //Variables
int Nreac;
int NEx;
double positionMUSIC=0.0;
double positionStrip=0.0;
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
double distRecoil ;
int i_strip_temp;
int remainingStrip;
double thetaExit;
double remainingDist;
double EtempLossbeam;
double Etempbeam;
double Extemp;
double Exmax;
double ELossBeamRef[Nstrip];
double sigmaStrip;
int chosenStrip;
double theta_scat=0;
double cs[Nstrip+2];
double errcs[Nstrip+2];
double thick_Eeff;
int ThickTargetYield;
double incEloss_thick;
double incEloss_thick2;
double incSP;

#endif
