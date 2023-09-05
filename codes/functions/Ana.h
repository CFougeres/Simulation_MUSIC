#ifndef ANA.H
#define ANA.H

#include "main.h"

string filesave;
TChain*  SimTree = new TChain("tree");
int nStat;

TCanvas* creac = new TCanvas("creac","creac",800,800);
TCanvas* ctraces = new TCanvas("ctraces","ctraces",800,800);
TH2F* PadE; TH2F* DimE;
auto legend = new TLegend(0.1, 0.2, 0.4, 0.4);


double dER[2];
double Er[2];
double bin_to_E = 0.1;
int bin_dER;
int binEr;
int max_beam_built;
TGraph* beam[MaxN];
int max_reac_built[MaxN];
TGraph* reac[MaxN][MaxN];
double calibEtot[Npad];
int compteur_beam=0;
int compteur_ev_an[MaxN]={0};
double Eav;
double dEtot=0;
double distX;
double distY;
int compteur_reaction=0;
int find_reaction=0;

#endif
