//********************************************
//
// ROOT macro of Ana.cpp
//
//Author: C. Fougères (2023)
//any remark(s) to be sent at chloe.fougeres@gmail.com
//********************************************
string path  = gSystem->pwd();

#include "functions/Ana.h"

using namespace std;

void Ana()
{
    printf("=====================================\n");
    printf("===             Ana               ===\n");
    printf("=====================================\n");
    
 
    extraction_inputs();
    chosenPad[0]  = param_inputs[20][0];    chosenPad[1]  = param_inputs[20][1];
    Nrecoil = param_inputs[11][0];
    lengthDe = param_inputs[21][0]+1;
    HalfHeight = param_inputs[1][0];
    DeadLayer = param_inputs[2][0];
    PadDim = param_inputs[3][0];
    DetectorDim[0] = param_inputs[5][0]; DetectorDim[1] = param_inputs[5][1];
    AnodeDim[0] = param_inputs[4][0]; AnodeDim[1] = param_inputs[4][1];
    NpadX = AnodeDim[0]/PadDim+1;            NpadY = AnodeDim[1]/PadDim+1;
    ReactionChoice = param_inputs[23][0];
    NReaction = param_inputs[24][0];

    filesave = pathTrees + Form("tree_pad%ix%i.root",int(PadDim),int(PadDim)); cout<<filesave<<endl;
    SimTree->Add(filesave.c_str());     nStat =SimTree->GetEntries();
    
    SimTree->SetBranchAddress("EcmReac", &EcmReac);SimTree->SetBranchAddress("errEcmReac", &errEcmReac);
    SimTree->SetBranchAddress("XReac", &XReac);SimTree->SetBranchAddress("YReac", &YReac);
    SimTree->SetBranchAddress("ReacKind", &ReacKind);// 0 beam 1==(a,a) 2==(a,p) 3==(a,n), 4==(a,2n), 5==(a,g)
    SimTree->SetBranchAddress("Epad", &Epad);
    
    dER[0]=0;dER[1]=150; bin_dER =(dER[1]-dER[0])/bin_to_E;
    Er[0]=0;Er[1]=400; binEr =(Er[1]-Er[0])/bin_to_E;
    PadE = new TH2F("PadE",";Pad X; Pad Y; E (MeV)", NpadX, 0, NpadX, NpadY,0, NpadY);
    DimE = new TH2F(" DimE",";X(mm);Y(mm); E (MeV)",  AnodeDim[0]/(PadDim)+1,0, AnodeDim[0], AnodeDim[1]/(PadDim)+1, 0,  AnodeDim[1]);
    ctraces->Divide(2,2);

    for(int i=0;i<nStat;i++){
        SimTree->GetEntry(i);
        find_reaction=0;
        if(ReacKind== ReactionChoice && XReac==chosenPad[0] && YReac==chosenPad[1]){find_reaction=1;compteur_reaction+=1;}
        if(compteur_reaction<=NReaction && find_reaction==1){
            for(int x=0;x<NpadX;x++){
                distX=double(PadDim*x);
                for(int y=0;y<NpadY;y++){
                    distY=double(PadDim*y);
                    
                    PadE->Fill(x,y,Epad[x][y]);    DimE->Fill(distX,distY,Epad[x][y]);
                }
            }
        }
    }
    gStyle->SetPalette(kThermometer);
    gStyle->SetOptStat("ne");
    ctraces->cd(1);    PadE->Draw("surf3");
    PadE->GetYaxis()->CenterTitle(); PadE->GetYaxis()->SetTitleSize(0.05);PadE->GetYaxis()->SetLabelSize(0.04);
    PadE->GetXaxis()->CenterTitle();PadE->GetXaxis()->SetTitleSize(0.05);PadE->GetXaxis()->SetLabelSize(0.04);
    ctraces->cd(2);    PadE->Draw("colz");
    PadE->GetYaxis()->CenterTitle(); PadE->GetYaxis()->SetTitleSize(0.05);PadE->GetYaxis()->SetLabelSize(0.04);
    PadE->GetXaxis()->CenterTitle();PadE->GetXaxis()->SetTitleSize(0.05);PadE->GetXaxis()->SetLabelSize(0.04);
    ctraces->cd(3);    DimE->Draw("surf3");
    DimE->GetYaxis()->CenterTitle();  DimE->GetYaxis()->SetTitleSize(0.05); DimE->GetYaxis()->SetLabelSize(0.04);
    DimE->GetXaxis()->CenterTitle(); DimE->GetXaxis()->SetTitleSize(0.05); DimE->GetXaxis()->SetLabelSize(0.04);
    ctraces->cd(4);    DimE->Draw("colz");
    DimE->GetYaxis()->CenterTitle();  DimE->GetYaxis()->SetTitleSize(0.05); DimE->GetYaxis()->SetLabelSize(0.04);
    DimE->GetXaxis()->CenterTitle(); DimE->GetXaxis()->SetTitleSize(0.05); DimE->GetXaxis()->SetLabelSize(0.04);
}
