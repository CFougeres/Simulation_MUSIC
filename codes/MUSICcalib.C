//********************************************
//
// ROOT macro of MUSICcalib.c
//
//Author: C. Fougères (2023)
//any remark(s) to be sent at chloe.fougeres@gmail.com
//********************************************
string path   = gSystem->pwd();

#include "functions/MUSICcalib.h"
using namespace std;
void MUSICcalib()
{
    printf("=====================================\n");
    printf("===           MUSICcalib          ===\n");
    printf("=====================================\n");
    
    //USER INPUTS
    extraction_inputs();
    EbeamEntry = param_inputs[1][0];  Nreac = param_inputs[11][0];
    Npressure = param_inputs[21][0] ;  FWHM_Si  = param_inputs[23][0];
    for(int r=0;r<Npressure;r++){
        PressureCalib[r]= param_inputs[22][r];
    }
    SPwindow = param_inputs[4][0] ; SPgas = param_inputs[4][1]; SPsi = param_inputs[4][2];
    for(int i=0;i<3;i++){
        ZAMBeam[i]= param_inputs[3][i];
        ZAMGas[i]= param_inputs[2][i];
    }
    //SI response
    srand(time(NULL));
    gRandom->SetSeed(0);
    Er[0]=0;  Er[1]=410.; binEr=(Er[1]-Er[0])/bin_to_E;
    for(int a=0;a<Npressure+1;a++){
        if(SPwindow==0){
            if( SPgas==0){nameaxis=Form("Windows_SRIM & Gas_SRIM;E (MeV);Counts/%ikeV",int(1000*bin_to_E));}
            if( SPgas==1){nameaxis=Form("Windows_SRIM & Gas_ATIMA;E (MeV);Counts/%ikeV",int(1000*bin_to_E));}
            if( SPgas>1){nameaxis=Form("Windows_SRIM & Gas_(ATIMA & SRIM);E (MeV);Counts/%ikeV",int(1000*bin_to_E));}
        }
        if(SPwindow==1){
            if( SPgas==0){nameaxis=Form("Windows_ATIMA & Gas_SRIM;E (MeV);Counts/%ikeV",int(1000*bin_to_E));}
            if( SPgas==1){nameaxis=Form("Windows_ATIMA & Gas_ATIMA;E (MeV);Counts/%ikeV",int(1000*bin_to_E));}
            if( SPgas>1){nameaxis=Form("Windows_ATIMA & Gas_(ATIMA & SRIM);E (MeV);Counts/%ikeV",int(1000*bin_to_E));}
        }
        beamPeak[a] =new TH1F(Form("beamPeak%i",a),nameaxis.c_str(),binEr,Er[0],Er[1]);
    }

    
    mBeam = ZAMBeam[1] * mu + ZAMBeam[2] ;   mRed[0] =  mBeam/mu;
    //********************************
    //Simulation at different pressure
    //********************************
    fileSP = pathSP + Form("srim/%i%iBeamTi.dat", int(ZAMBeam[1]),int(ZAMBeam[0]));
    extraction_sp(fileSP, 0, EnergySp[0][0], SPel[0][0],  SPnu[0][0]);
    fileSP = pathSP + Form("lise/%i%iBeamTi.dat",int(ZAMBeam[1]),int(ZAMBeam[0]));
    extraction_sp(fileSP, 1, EnergySp[0][1], SPel[0][1],  SPnu[0][1]);
    fileSP = pathSP + Form("srim/%i%iBeamSi.dat", int(ZAMBeam[1]),int(ZAMBeam[0]));
    extraction_sp(fileSP, 0, EnergySp[2][0], SPel[2][0],  SPnu[2][0]);
    fileSP = pathSP + Form("lise/%i%iBeamSi.dat",int(ZAMBeam[1]),int(ZAMBeam[0]));
    extraction_sp(fileSP, 1, EnergySp[2][1], SPel[2][1],  SPnu[2][1]);
    
    for(int a=0;a<Npressure+1;a++){
            fileSP = pathSP + Form("srim/%i%iBeam%i%iGas%iTorr.dat", int(ZAMBeam[1]),int(ZAMBeam[0]), int(ZAMGas[1]),int(ZAMGas[0]), PressureCalib[a-1]);
            extraction_sp(fileSP, 0, EnergySp[1][0], SPel[1][0],  SPnu[1][0]);
            fileSP = pathSP + Form("lise/%i%iBeam%i%iGas%iTorr.dat", int(ZAMBeam[1]),int(ZAMBeam[0]), int(ZAMGas[1]),int(ZAMGas[0]), PressureCalib[a-1]);
            extraction_sp(fileSP, 1, EnergySp[1][1], SPel[1][1],  SPnu[1][1]);
        
        for(int r=0;r<Nreac;r++){
            Ebeam = EbeamEntry;
            if(a==0){theta_scat=  scattered_theta(0);}
            if(a>0){theta_scat=  scattered_theta(PressureCalib[a-1]);}
            //********************************
            //loss in windows
            //********************************
            if(SPwindow<2){
                elosswindow = loss_E(SPwindow, EbeamEntry, mRed[0], entranceWindow_width, dxWindow , EnergySp[0][SPwindow],SPel[0][SPwindow],  SPnu[0][SPwindow]);
            }
            if(SPwindow>1){
                elosswindow = loss_E(0, EbeamEntry, mRed[0], entranceWindow_width, dxWindow , EnergySp[0][0],SPel[0][0],  SPnu[0][0]);
                elosswindow = (elosswindow + loss_E(1, EbeamEntry,mRed[0],  entranceWindow_width,dxWindow ,EnergySp[0][1],SPel[0][1], SPnu[0][1]))/2.0;
            }
            cout<< " Loss in window "<< elosswindow <<endl;
            Ebeam=Ebeam-elosswindow;
            if(a>0){
                //********************************
                //loss in gas
                //********************************
                if(SPgas<2){
                    ElossBeam = loss_E(SPgas, Ebeam,mRed[0], MUSIC_length/cos(theta_scat*TMath::Pi()/180.), dxGas, EnergySp[1][SPgas], SPel[1][SPgas], SPnu[1][SPgas]);
                }
                if(SPgas>1){
                    ElossBeam = loss_E(0, Ebeam,mRed[0], MUSIC_length/cos(theta_scat*TMath::Pi()/180.), dxGas, EnergySp[1][0],SPel[1][0],  SPnu[1][0]);
                    ElossBeam = (ElossBeam + loss_E(1,Ebeam,mRed[0], MUSIC_length/cos(theta_scat*TMath::Pi()/180.), dxGas,EnergySp[1][1],SPel[1][1], SPnu[1][1]))/2.0;
                }
                cout<<"Loss in gas "<<ElossBeam<<endl;
                Ebeam = Ebeam - ElossBeam;
            }
            //********************************
            //loss in windows
            //********************************
            if(SPwindow<2){
                elosswindow = loss_E(SPwindow, Ebeam, mRed[0], entranceWindow_width/cos(theta_scat*TMath::Pi()/180.), dxWindow , EnergySp[0][SPwindow],SPel[0][SPwindow],  SPnu[0][SPwindow]);
            }
            if(SPwindow>1){
                elosswindow = loss_E(0,  Ebeam, mRed[0], entranceWindow_width/cos(theta_scat*TMath::Pi()/180.), dxWindow , EnergySp[0][0],SPel[0][0],  SPnu[0][0]);
                elosswindow = (elosswindow + loss_E(1, Ebeam,mRed[0],  entranceWindow_width/cos(theta_scat*TMath::Pi()/180.), dxWindow ,EnergySp[0][1],SPel[0][1], SPnu[0][1]))/2.0;
            }
            cout<< " Loss in window "<< elosswindow <<endl;
            Ebeam=Ebeam-elosswindow;
            //********************************
            //loss dead layer in Si
            //********************************
            if(SPsi<2){
                elossSi = loss_E(SPsi, Ebeam, mRed[0],si_dead_layer_width/cos(theta_scat*TMath::Pi()/180.), dxSi , EnergySp[2][SPsi],SPel[2][SPsi],  SPnu[2][SPsi]);
            }
            if(SPsi>1){
                elossSi= loss_E(0,  Ebeam, mRed[0],si_dead_layer_width/cos(theta_scat*TMath::Pi()/180.),  dxSi , EnergySp[2][0],SPel[2][0],  SPnu[2][0]);
                elossSi = (  elossSi + loss_E(1, Ebeam,mRed[0], si_dead_layer_width/cos(theta_scat*TMath::Pi()/180.), dxSi,EnergySp[2][1],SPel[2][1], SPnu[2][1]))/2.0;
            }
            Ebeam=Ebeam-elossSi;
            cout<<"loss in si "<<elossSi<<std::endl;
            beamPeak[a]->Fill(rESi->Gaus(Ebeam,(FWHM_Si*Ebeam/100.)/2.35));
        }
    }
    //********************************
    //PLOTTING
    //********************************
    cEloss->cd(); auto leg = new TLegend(0.1, 0.2, 0.4, 0.4);
    beamPeak[0]->Draw("C");beamPeak[0]->SetLineColor(1);
    beamPeak[0]->GetYaxis()->CenterTitle();beamPeak[0]->GetXaxis()->CenterTitle();
    leg->AddEntry(beamPeak[0],Form("P=%i Torr, entrance wind.",0),"l");
    for(int p=1;p<Npressure+1;p++){
        col = 2 + p;
        beamPeak[p]->Draw("Csame");beamPeak[p]->SetLineColor(col);
        leg->AddEntry(beamPeak[p],Form("P=%i Torr",PressureCalib[p-1],"l"));
    }
    leg->SetNColumns(3);
    leg->Draw("same");
}
