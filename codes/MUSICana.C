//********************************************
//
// ROOT macro of MUSICana.cpp
//
//Author: C. Fougères (2023)
//any remark(s) to be sent at chloe.fougeres@gmail.com
//********************************************
string path  = gSystem->pwd();

#include "functions/MUSICana.h"

using namespace std;

void MUSICana()
{
    printf("=====================================\n");
    printf("===       MUSICana         ===\n");
    printf("=====================================\n");
    
    extraction_inputs();
    chosenStrip  = param_inputs[15][0];    Nrecoil = param_inputs[6][0];
    lengthDe = param_inputs[16][0]+1;
    for(int i=0;i<3;i++){
        ZAMBeam[i]= param_inputs[3][i];
        ZAMGas[i]= param_inputs[2][i];
    }
    for(int r=0;r<Nrecoil ;r++){
        ZAMRecoil[r][0]= param_inputs[7][r];
        ZAMRecoil[r][1]= param_inputs[8][r];
        ZAMRecoil[r][2]= param_inputs[9][r];
        ExRecoil[r]= param_inputs[10][0];
    }
    Elim_plot[0]= param_inputs[24][0];     Elim_plot[1]= param_inputs[24][1];
    filesave = pathTrees + Form("tree_%i%iBeam_%i%iGas.root",int(ZAMBeam[1]),int(ZAMBeam[0]),int(ZAMGas[1]),int(ZAMGas[0])); cout<<filesave<<endl;
    SimTree->Add(filesave.c_str());     nStat =SimTree->GetEntries();

    
    SimTree->Branch("stripReac", &stripReac, "stripReac/I");
    SimTree->SetBranchAddress("stripReac", &stripReac);
    SimTree->Branch("ReacKind", &ReacKind,"ReacKind/I");// 0 beam 1==(a,a) 2==(a,p) 3==(a,n), 4==(a,2n), 5==(a,g)
    SimTree->SetBranchAddress("ReacKind", &ReacKind);
    SimTree->Branch("Estrip", &Estrip,"Estrip[18]/F");
    SimTree->SetBranchAddress("Estrip", &Estrip);
    SimTree->Branch("Istrip", &Istrip, "Istrip[18]/I");
    SimTree->SetBranchAddress("Istrip", &Istrip);
    
    dER[0]=Elim_plot[0];dER[1]=Elim_plot[1]*lengthDe; bin_dER =(dER[1]-dER[0])/bin_to_E;
    Er[0]=Elim_plot[0];Er[1]=Elim_plot[1]*Nstrip; binEr =(Er[1]-Er[0])/bin_to_E;
    dee = new TH2F("dee",Form("Strip %i;#Sigma_{0}^{17}#DeltaE_{i} (MeV);#Sigma_{%i}^{%i}#DeltaE_{i} (MeV)",chosenStrip, chosenStrip+1, int(chosenStrip+lengthDe)), binEr,Er[0],Er[1], bin_dER,dER[0],dER[1]);
    max_beam_built=20;
    
    for(int r=0;r<Nrecoil;r++){max_reac_built[r]=20;}
    for(int i=0;i<nStat;i++){
        SimTree->GetEntry(i);
        dEtot=0.;  Eav=0;
        for(int cc=0;cc<Nstrip;cc++){
            calibEtot[cc]= Estrip[cc];
            dEtot+=calibEtot[cc];
        }
        for(int l=1;l<lengthDe;l++){Eav += calibEtot[chosenStrip+l];}
        if(ReacKind == 0 && compteur_beam<max_beam_built){
            std::cout<<"beam "<<compteur_beam<<std::endl;
            beam[compteur_beam] = new TGraph(Nstrip,strip, calibEtot);
            compteur_beam = compteur_beam+1;
        }
        if(stripReac==chosenStrip && compteur_ev_an[ReacKind-1]<max_reac_built[ReacKind-1]){
            if(ZAMBeam[1]==87 && ZAMGas[1]==4){
                if(ReacKind == 1 || ReacKind == 2||ReacKind == 5){
                    reac[ReacKind-1][compteur_ev_an[ReacKind-1]] = new TGraph(Nstrip,strip, calibEtot);
                    compteur_ev_an[ReacKind-1] = compteur_ev_an[ReacKind-1]+1;
                }
                if(ReacKind == 3|| ReacKind == 4) {//
                    reac[2][compteur_ev_an[2]] = new TGraph(Nstrip,strip, calibEtot);
                    compteur_ev_an[2] = compteur_ev_an[2]+1;
                }
            }
            if(ZAMBeam[1]==14 && ZAMGas[1]==4){
                for(int nc=0;nc< Nrecoil;nc++){
                    if(ReacKind == nc+1) {//
                        reac[nc][compteur_ev_an[nc]] = new TGraph(Nstrip,strip, calibEtot);
                        compteur_ev_an[nc] = compteur_ev_an[nc]+1;
                    }
                }
            }
        }
      /*  if(ReacKind>-1){
            dee->Fill(dEtot,Eav);
        }
        */
        if(stripReac==chosenStrip){
            if(ZAMBeam[1]==87 && ZAMGas[1]==4){
                if( ReacKind == 1||  ReacKind == 3 || ReacKind == 4){                    dee->Fill(dEtot,Eav);                }
            }
            if(ZAMBeam[1]==14 && ZAMGas[1]==4){
                for(int nc=0;nc< Nrecoil;nc++){                    if(ReacKind == nc+1) {dee->Fill(dEtot,Eav);}                }
            }
        }
    }
    gStyle->SetPalette(kThermometer);
    creac->cd();dee->Draw("colz");
    dee->GetYaxis()->CenterTitle(); dee->GetYaxis()->SetTitleSize(0.05);dee->GetYaxis()->SetLabelSize(0.04);
    dee->GetXaxis()->CenterTitle();dee->GetXaxis()->SetTitleSize(0.05);dee->GetXaxis()->SetLabelSize(0.04);
    ctraces->cd();
    beam[0]->Draw("AC");beam[0]->GetYaxis()->SetRangeUser(Elim_plot[0],Elim_plot[1]);
    beam[0]->GetXaxis()->SetRangeUser(0,17);
    beam[0]->SetTitle(";Strip; #DeltaE (MeV)");
    beam[0]->GetYaxis()->CenterTitle(); beam[0]->GetYaxis()->SetTitleSize(0.05);beam[0]->GetYaxis()->SetLabelSize(0.04);
    beam[0]->GetXaxis()->CenterTitle();beam[0]->GetXaxis()->SetTitleSize(0.05);beam[0]->GetXaxis()->SetLabelSize(0.04);
    beam[0]->SetLineColor(1);
    for(int a=1;a<max_beam_built;a++){
        if(a<compteur_beam){
            beam[a]->Draw("same");beam[a]->SetLineColor(1);
        }
        std::cout<<a<<std::endl;
        if(ZAMBeam[1]==87 && ZAMGas[1]==4){
            if(a-1<compteur_ev_an[4]){//(a,g)
                reac[4][a-1]->Draw("same"); reac[4][a-1]->SetLineColor(8);
            }
            if(a-1<compteur_ev_an[2]){//(a,xn)
                reac[2][a-1]->Draw("same"); reac[2][a-1]->SetLineColor(2);
            }
            if(a-1<compteur_ev_an[1]){//(a,p)
                reac[1][a-1]->Draw("same"); reac[1][a-1]->SetLineColor(4);
            }
            if(a-1<compteur_ev_an[0]){//(a,a')
                reac[0][a-1]->Draw("same"); reac[0][a-1]->SetLineColor(2);reac[0][a-1]->SetLineStyle(2);
            }
        }
        if(ZAMBeam[1]==14 && ZAMGas[1]==4){
            for(int nc=0;nc<Nrecoil;nc++){
                if(a-1<compteur_ev_an[nc] && a-1>-1){
                    reac[nc][a-1]->Draw("same");  reac[nc][a-1]->SetLineColor(nc+2);
                }
            }
        }
    }
    if(ZAMBeam[1]==87 && ZAMGas[1]==4){
        legend->AddEntry(beam[0],"^{87}Rb beam");
        legend->AddEntry(reac[2][0],"^{87}Rb(#alpha,xn)Y");
        legend->AddEntry(reac[1][0],"^{87}Rb(#alpha,p)^{90}Sr");
        legend->AddEntry(reac[4][0],"^{87}Rb(#alpha,#gamma)^{91}Y");
        legend->AddEntry(reac[0][0],"^{87}Rb(#alpha,#alpha')^{87}Rb");
    }
    if(ZAMBeam[1]==14 && ZAMGas[1]==4){
        legend->AddEntry(beam[0],"^{14}O beam");
        if(compteur_ev_an[0]>1){   legend->AddEntry(reac[0][0],"^{14}O(#alpha,#alpha')^{14}O");}
        if(compteur_ev_an[1]>1){   legend->AddEntry(reac[1][0],"^{14}O(#alpha,p)^{17}F");}
        if(compteur_ev_an[2]>1){   legend->AddEntry(reac[2][0],"^{14}O(#alpha,2p)^{16}O");}
        if(compteur_ev_an[3]>1){    legend->AddEntry(reac[3][0],"^{14}O(#alpha,n)^{17}Ne");}
        if(compteur_ev_an[4]>1){    legend->AddEntry(reac[4][0],"^{14}O(#alpha,d)^{16}F");}
        if(compteur_ev_an[5]>1){   legend->AddEntry(reac[5][0],"^{14}O(#alpha,#gamma)^{18}Ne");}
    }
    legend->Draw("same");
}
