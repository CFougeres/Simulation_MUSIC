//********************************************
//
// ROOT macro of Sim.cpp
//
//Author: C. Fougères (2023)
//any remark(s) to be sent at chloe.fougeres@gmail.com
//********************************************
string path   = gSystem->pwd();

#include "functions/Sim.h"
using namespace std;
void Sim()
{
    printf("=====================================\n");
    printf("===              Sim              ===\n");
    printf("=====================================\n");
    
    //USER INPUTS
    extraction_inputs();
    plotlive= param_inputs[25][0];
    HalfHeight = param_inputs[1][0];
    DeadLayer = param_inputs[2][0];
    PadDim = param_inputs[3][0];
    DetectorDim[0] = param_inputs[5][0]; DetectorDim[1] = param_inputs[5][1];
    AnodeDim[0] = param_inputs[4][0]; AnodeDim[1] = param_inputs[4][1];
    NpadX = AnodeDim[0]/ PadDim+1;            NpadY = AnodeDim[1]/ PadDim+1;
    incSP = param_inputs[22][0];
    sigmaPad= param_inputs[19][0]/2.38;      pressure = param_inputs[10][0] ;
    SPwindow = param_inputs[9][0] ; SPgas = param_inputs[9][1];
    EbeamEntry = param_inputs[6][0]; Nrecoil = param_inputs[11][0];
    Nreac = param_inputs[16][0];  NEx = param_inputs[17][0];
    dx = param_inputs[18][0];
    for(int i=0;i<3;i++){
        ZAMBeam[i]= param_inputs[8][i];
        ZAMGas[i]= param_inputs[7][i];
    }
    for(int r=0;r<Nrecoil ;r++){
        ZAMRecoil[r][0]= param_inputs[12][r];
        ZAMRecoil[r][1]= param_inputs[13][r];
        ZAMRecoil[r][2]= param_inputs[14][r];
        ExRecoil[r]= param_inputs[15][0];
    }
    tree->Branch("EcmReac", &EcmReac,"EcmReac/F");    tree->SetBranchAddress("EcmReac", &EcmReac);
    tree->Branch("errEcmReac", &errEcmReac,"errEcmReac/F");    tree->SetBranchAddress("errEcmReac", &errEcmReac);
    tree->Branch("XReac", &XReac,"XReac/I");    tree->SetBranchAddress("XReac", &XReac);
    tree->Branch("YReac", &YReac,"YReac/I");    tree->SetBranchAddress("YReac", &YReac);
    tree->Branch("ReacKind", &ReacKind,"ReacKind/I");// 0 beam 1==(a,a) 2==(a,p) 3==(a,n), 4==(a,2n), 5==(a,g)
    tree->SetBranchAddress("ReacKind", &ReacKind);
    tree->Branch("Epad", &Epad,"Epad[1024][1024]/F");    tree->SetBranchAddress("Epad", &Epad);
    //********************************
    //KINEMATICS
    //********************************
    Ebeam = EbeamEntry;
    reaction_balance();
    gSystem->Load("libPhysics");
    Pfais = sqrt(EbeamEntry * EbeamEntry + 2. * EbeamEntry * mBeam);
    TLorentzVector target(0.0, 0.0, 0.0, mGas * 0.001);
    TLorentzVector beam(0.0, 0.0, Pfais * 0.001, (EbeamEntry + mBeam) * 0.001);
    TLorentzVector Tr = beam + target;
    TLorentzVector * pRecoil;
    TLorentzVector * pEjectil;
    for(int r=0;r<Nrecoil;r++){event[r].SetDecay(Tr, 2, masses[r]);}
    randEx = new TRandom();     randEloss = new TRandom();  srand(time(NULL));      gRandom->SetSeed(0);
    TrajectoryEjectil = new TH2F("TrajectoryEjectil",";Pad X; Pad Y", NpadX, 0, NpadX, NpadY,0, NpadY);
    TrajectoryRecoil = new TH2F("TrajectoryRecoil",";Pad X; Pad Y", NpadX, 0, NpadX, NpadY,0, NpadY);
    ctraces->Divide(1,2);
    gStyle->SetOptStat("ne");
    gStyle->SetPalette(kThermometer);
    ctraces->cd(1);  TrajectoryRecoil->Draw("scat");
    TrajectoryRecoil->GetYaxis()->CenterTitle(); TrajectoryRecoil->GetYaxis()->SetTitleSize(0.05); TrajectoryRecoil->GetYaxis()->SetLabelSize(0.04);
    TrajectoryRecoil->GetXaxis()->CenterTitle(); TrajectoryRecoil->GetXaxis()->SetTitleSize(0.05); TrajectoryRecoil->GetXaxis()->SetLabelSize(0.04);
    ctraces->cd(2);  TrajectoryEjectil->Draw("scat");
    TrajectoryEjectil->GetYaxis()->CenterTitle(); TrajectoryEjectil->GetYaxis()->SetTitleSize(0.05); TrajectoryEjectil->GetYaxis()->SetLabelSize(0.04);
    TrajectoryEjectil->GetXaxis()->CenterTitle(); TrajectoryEjectil->GetXaxis()->SetTitleSize(0.05); TrajectoryEjectil->GetXaxis()->SetLabelSize(0.04);
  
    
    
    //********************************
    //STOPPING POWERS
    //********************************
    for(int s=0;s<Nrecoil+3;s++){
        if(s==0){
            fileSP = pathSP + Form("srim/%i%iBeamTi.dat", int(ZAMBeam[1]),int(ZAMBeam[0]));
            extraction_sp(fileSP, 0, EnergySp[s][0], SPel[s][0],  SPnu[s][0]);
            fileSP = pathSP + Form("lise/%i%iBeamTi.dat",int(ZAMBeam[1]),int(ZAMBeam[0]));
            extraction_sp(fileSP, 1, EnergySp[s][1], SPel[s][1],  SPnu[s][1]);
        }
        if(s==1){
            fileSP = pathSP + Form("srim/%i%iBeam%i%iGas%iTorr.dat", int(ZAMBeam[1]),int(ZAMBeam[0]), int(ZAMGas[1]),int(ZAMGas[0]), int(pressure));
            extraction_sp(fileSP, 0, EnergySp[s][0], SPel[s][0],  SPnu[s][0]);
            fileSP = pathSP + Form("lise/%i%iBeam%i%iGas%iTorr.dat", int(ZAMBeam[1]),int(ZAMBeam[0]), int(ZAMGas[1]),int(ZAMGas[0]), int(pressure));
            extraction_sp(fileSP, 1, EnergySp[s][1], SPel[s][1],  SPnu[s][1]);
        }
        if(s==2){
            fileSP = pathSP + Form("srim/%i%iBeam%i%iGas%iTorr.dat", int(ZAMGas[1]),int(ZAMGas[0]),int( ZAMGas[1]),int(ZAMGas[0]),int(pressure));
            extraction_sp(fileSP, 0, EnergySp[s][0], SPel[s][0],  SPnu[s][0]);
            fileSP = pathSP + Form("lise/%i%iBeam%i%iGas%iTorr.dat", int(ZAMGas[1]),int(ZAMGas[0]), int(ZAMGas[1]),int(ZAMGas[0]), int(pressure));
            extraction_sp(fileSP, 1, EnergySp[s][1], SPel[s][1],  SPnu[s][1]);
        }
        if(s>2){
            fileSP = pathSP + Form("srim/%i%iBeam%i%iGas%iTorr.dat", int(ZAMRecoil[s-3][1]),int(ZAMRecoil[s-3][0]), int(ZAMGas[1]),int(ZAMGas[0]), int(pressure));
            extraction_sp(fileSP, 0, EnergySp[s][0], SPel[s][0],  SPnu[s][0]);
            fileSP = pathSP + Form("lise/%i%iBeam%i%iGas%iTorr.dat",int(ZAMRecoil[s-3][1]),int(ZAMRecoil[s-3][0]), int(ZAMGas[1]),int(ZAMGas[0]), int(pressure));
            extraction_sp(fileSP, 1, EnergySp[s][1], SPel[s][1],  SPnu[s][1]);
            fileSP = pathSP + Form("srim/%i%iBeam%i%iGas%iTorr.dat",int(ZAMultEjectil[s-3][1]),int(ZAMultEjectil[s-3][0]), int(ZAMGas[1]),int(ZAMGas[0]), int(pressure));
            extraction_sp(fileSP, 0, EnergySp[s][2], SPel[s][2],  SPnu[s][2]);
            fileSP = pathSP + Form("lise/%i%iBeam%i%iGas%iTorr.dat",int(ZAMultEjectil[s-3][1]),int(ZAMultEjectil[s-3][0]), int(ZAMGas[1]),int(ZAMGas[0]), int(pressure));
            extraction_sp(fileSP, 1, EnergySp[s][3], SPel[s][3],  SPnu[s][3]);
        }
    }
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
    //********************************
    //loss dead layer after Ti window
    //********************************
  if(SPgas<2){
        ElossBeam = loss_E(SPgas, Ebeam,mRed[0],  DeadLayer, dxGas, EnergySp[1][SPgas], SPel[1][SPgas], SPnu[1][SPgas]);
    }
    if(SPgas>1){
        ElossBeam = loss_E(0, Ebeam,mRed[0],  DeadLayer, dxGas, EnergySp[1][0],SPel[1][0],  SPnu[1][0]);
        ElossBeam = (ElossBeam + loss_E(1,Ebeam,mRed[0],  DeadLayer, dxGas,EnergySp[1][1],SPel[1][1], SPnu[1][1]))/2.0;
    }
    cout<<"Loss in dead layer "<<ElossBeam<<endl;
    Ebeam = Ebeam - ElossBeam;
    EbeamEntry = Ebeam ;
  
    comptPadX=0; comptPadY=0;
    //********************************
    //Derivation beam trajectory
    //********************************
    XReac=-1; YReac=-1; EcmReac=-1; errEcmReac=-1;

    for(int j=0;j<Nreac;j++){
        Ebeam = EbeamEntry;
        comptPadX=0; comptPadY=int(NpadY/2);
        while(comptPadX<NpadX){
            Etempbeam = Ebeam;
            if(SPgas<2){ElossBeam = loss_E(SPgas, Ebeam,mRed[0], PadDim, dxGas, EnergySp[1][SPgas], SPel[1][SPgas], SPnu[1][SPgas]);}
            if(SPgas>1){
                ElossBeam = loss_E(0, Ebeam,mRed[0],PadDim, dxGas, EnergySp[1][0],SPel[1][0],  SPnu[1][0]);
                ElossBeam = (ElossBeam + loss_E(1,Ebeam,mRed[0],PadDim, dxGas,EnergySp[1][1],SPel[1][1], SPnu[1][1]))/2.0;
            }
            Ebeam = Ebeam - ElossBeam;
            ELossBeamRef[comptPadX][comptPadY]=ElossBeam;
            Epad[comptPadX][comptPadY]=randEloss->Gaus(ElossBeam,sigmaPad);
            comptPadX+=1;
        }
        ReacKind=0;
        tree->Fill();
    }

    //********************************
    //Derivation reaction trajectories
    //********************************
    Ebeam = EbeamEntry;
    comptPadX=0; comptPadY=0;
    XReac=0; YReac=int(NpadY/2);
    cout<<"beam entry anode "<< EbeamEntry<<endl;
    while(XReac<NpadX){
        cout<<"padX reac ongoing "<<XReac<<endl;
        positionPad=0.0;  EtempLossbeam=0.0;    Etempbeam=Ebeam;
        while(positionPad<PadDim-dx){
            cout<<"padX position "<<positionPad<<endl;
            positionPad = positionPad+dx;
            if(SPgas<2){
                EtempLossbeam = loss_E(SPgas,Etempbeam, mRed[0],  positionPad , dxGas, EnergySp[1][SPgas], SPel[1][SPgas], SPnu[1][SPgas]);
                ElossBeam = loss_E(SPgas,Ebeam,mRed[0], dx, dxGas, EnergySp[1][SPgas], SPel[1][SPgas], SPnu[1][SPgas]);
            }
            if(SPgas>1){
                EtempLossbeam = loss_E(0, Etempbeam,mRed[0], positionPad , dxGas, EnergySp[1][0],SPel[1][0],  SPnu[1][0]);
                EtempLossbeam = ( EtempLossbeam + loss_E(1,Etempbeam,mRed[0], positionPad , dxGas,EnergySp[1][1],SPel[1][1], SPnu[1][1]))/2.0;
                ElossBeam  = loss_E(0, Ebeam,mRed[0],dx, dxGas, EnergySp[1][0],SPel[1][0],  SPnu[1][0]);
                ElossBeam  = ( ElossBeam  + loss_E(1,Ebeam,mRed[0],dx, dxGas,EnergySp[1][1],SPel[1][1], SPnu[1][1]))/2.0;
            }
            Ebeam = Ebeam - ElossBeam;
            EcmReac=((mGas)/(mBeam+mGas))*Ebeam;
            errEcmReac= TMath::Abs(EcmReac-((mGas)/(mBeam+mGas))*((1.0-incSP)*Ebeam));
            vbeamlab = sqrt(1 - 1 / pow((Ebeam / mBeam + 1), 2));
            Pfais = sqrt(Ebeam * Ebeam + 2. * Ebeam * mBeam);
            beam.SetPxPyPzE(0, 0, Pfais * 0.001, (Ebeam+mBeam) * 0.001);
            Tr = beam + target;
            for(int r=0;r<5;r++){
                ReacKind=r+1;
                for(int k=0;k<NEx;k++){
                    Exmax =  EcmReac + Qval[r];
                    Extemp =  0.0 + Exmax*randEx->Rndm(); ;
                    if( Exmax<0.0){ Extemp=0;}
                    masses[r][0] = 0.001 *(mRecoil[r]+ Extemp);
                    event[r].SetDecay(Tr, 2, masses[r]);
                    for(int j=0;j<Nreac;j++){
                        for(int jj=0;jj<NpadX;jj++){for(int jjj=0;jjj<NpadY;jjj++){Epad[jj][jjj]= 0;}}
                        for(int jj=0;jj<XReac;jj++){
                            Epad[jj][YReac]= randEloss->Gaus(ELossBeamRef[jj][YReac],sigmaPad);
                        }
                        event[r].Generate();
                        pRecoil = event[r].GetDecay(0); pEjectil = event[r].GetDecay(1);
                        thetaRecoil_lab = pRecoil->Theta();//Rad
                        vRecoil_lab =pRecoil->Beta();
                        eRecoil_lab=(1 / sqrt(1 - pow(vRecoil_lab, 2)) - 1) * mRecoil[r];
                        if(ZAMultEjectil[r][0]>0){
                            thetaEjectil_lab = pEjectil->Theta();//Rad
                            vEjectil_lab =pEjectil->Beta();
                            eEjectil_lab=(1 / sqrt(1 - pow(vEjectil_lab, 2)) - 1) * mEjectil[r];
                        }
                        ElossEjectil=0;ElossRecoil=0;
                        if(eRecoil_lab>0.015){
                            distRecoil = (PadDim-positionPad)/cos(thetaRecoil_lab);
                            if((distRecoil*sin(thetaRecoil_lab))<PadDim/2.){
                                EntrancePad[0]=0.;  EntrancePad[1]=PadDim/2.+distRecoil*sin(thetaRecoil_lab);
                            }
                            if((distRecoil*sin(thetaRecoil_lab))==PadDim/2.){
                                EntrancePad[0]=0.; EntrancePad[1]=0;
                            }
                            if((distRecoil*sin(thetaRecoil_lab))>PadDim/2.){
                                EntrancePad[1]=0.;
                                distRecoil =(PadDim/2.)/sin(thetaRecoil_lab);
                                EntrancePad[0]=positionPad+distRecoil*cos(thetaRecoil_lab);
                            }
                            if(SPgas<2){
                                ElossRecoil = loss_E(SPgas,eRecoil_lab,mRed[r], distRecoil, dxGas, EnergySp[r+3][SPgas], SPel[r+3][SPgas], SPnu[r+3][SPgas]);
                            }
                            if(SPgas>1){
                                ElossRecoil  = loss_E(0, eRecoil_lab,mRed[r],distRecoil, dxGas, EnergySp[r+3][0],SPel[r+3][0],  SPnu[r+3][0]);
                                ElossRecoil  = (ElossRecoil  + loss_E(1,eRecoil_lab,mRed[r], distRecoil, dxGas,EnergySp[r+3][1],SPel[r+3][1], SPnu[r+3][1]))/2.0;
                            }
                            Epad[XReac][YReac]=randEloss->Gaus(ElossRecoil + EtempLossbeam,sigmaPad);
                            eRecoil_lab = eRecoil_lab - ElossRecoil;
                            comptPadY = YReac; comptPadX = XReac;
                            if(ZAMultEjectil[r][0]>0 && plotlive==1){
                                TrajectoryRecoil->Fill(comptPadX,comptPadY);
                                ctraces->cd(1);
                                TrajectoryRecoil->Draw("cont4");
                                ctraces->Update();
                                 gSystem->ProcessEvents();
                                gSystem->Sleep(break_trajectory);
                            }
                            while(comptPadX <NpadX &&  comptPadY<NpadY &&  eRecoil_lab>0.015){
                                distRecoil = (PadDim-EntrancePad[0])/cos(thetaRecoil_lab);
                                if((distRecoil*sin(thetaRecoil_lab))<(PadDim-EntrancePad[1])){
                                    EntrancePad[0]=0.;  EntrancePad[1]=EntrancePad[1]+distRecoil*sin(thetaRecoil_lab);
                                    comptPadX+=1;
                                }else if((distRecoil*sin(thetaRecoil_lab))==(PadDim-EntrancePad[1])){
                                    EntrancePad[0]=0.; EntrancePad[1]=0;
                                    comptPadY+=1; comptPadX+=1;
                                }else{
                                    distRecoil =(PadDim-EntrancePad[1])/sin(thetaRecoil_lab);
                                    EntrancePad[1]=0.;   EntrancePad[0]=EntrancePad[0]+distRecoil*cos(thetaRecoil_lab);
                                    comptPadY+=1;
                                }
                                if(SPgas<2){
                                    ElossRecoil = loss_E(SPgas,eRecoil_lab,mRed[r], distRecoil, dxGas, EnergySp[r+3][SPgas], SPel[r+3][SPgas], SPnu[r+3][SPgas]);
                                }
                                if(SPgas>1){
                                    ElossRecoil  = loss_E(0,eRecoil_lab,mRed[r],distRecoil, dxGas, EnergySp[r+3][0],SPel[r+3][0],  SPnu[r+3][0]);
                                    ElossRecoil  = (ElossRecoil  + loss_E(1,eRecoil_lab,mRed[r], distRecoil, dxGas,EnergySp[r+3][1],SPel[r+3][1], SPnu[r+3][1]))/2.0;
                                }
                           //     cout<<"Epad recoil filling "<<comptPadX<<" "<<NpadX<<" "<<comptPadY<<" "<<NpadY<<endl;
                                Epad[comptPadX][comptPadY]=Epad[comptPadX][comptPadY]+randEloss->Gaus(ElossRecoil,sigmaPad);
                                eRecoil_lab = eRecoil_lab - ElossRecoil;
                                if(ZAMultEjectil[r][0]>0 && plotlive==1){
                                    TrajectoryRecoil->Fill(comptPadX,comptPadY);
                                    ctraces->cd(1);
                                    TrajectoryRecoil->Draw("cont4");
                                    ctraces->Update();
                                     gSystem->ProcessEvents();
                                    gSystem->Sleep(break_trajectory);
                                }
                            }
                            //EJECTIL
                            if(ZAMultEjectil[r][0]>0 && eEjectil_lab>0.015){
                                comptPadY = YReac; comptPadX = XReac;
                                ReacPad[0]=positionPad;ReacPad[1]=PadDim/2.;
                                vect=sqrt(pow(PadDim- ReacPad[0],2)+pow(PadDim- ReacPad[1],2));
                                vect1=sqrt(pow(PadDim- ReacPad[0],2)+pow(ReacPad[1],2));
                                vect2=sqrt(pow(ReacPad[0],2)+pow(PadDim- ReacPad[1],2));
                                vect3==sqrt(pow(ReacPad[0],2)+pow(ReacPad[1],2));
                                //1
                                if(thetaEjectil_lab<acos((PadDim-ReacPad[0])/vect) || thetaEjectil_lab>2*pi_val-acos((PadDim- ReacPad[0])/vect1)){
                                    comptPadX+=1;
                                    EntrancePad[0]=0.;
                                    distEjectil = TMath::Abs((PadDim-ReacPad[0])/cos(thetaEjectil_lab));
                                    EntrancePad[1]= ReacPad[1]+distEjectil*sin(thetaEjectil_lab);
                                }
                                //2
                                if(thetaEjectil_lab==acos((PadDim-ReacPad[0])/vect)){
                                    comptPadX+=1;comptPadY+=1;
                                    EntrancePad[0]=0.;EntrancePad[1]=0.;
                                    distEjectil = TMath::Abs((PadDim-ReacPad[0])/cos(thetaEjectil_lab));
                                }
                                //3
                                if(thetaEjectil_lab>acos((PadDim-ReacPad[0])/vect) && thetaEjectil_lab<(pi_val-acos(ReacPad[0]/vect2))){
                                    comptPadY+=1;
                                    EntrancePad[1]=0.;
                                    distEjectil =TMath::Abs((PadDim-ReacPad[1])/sin(thetaEjectil_lab));
                                    EntrancePad[0] = ReacPad[0] +distEjectil*cos(thetaEjectil_lab);
                                }
                                //4
                                if(thetaEjectil_lab==(pi_val-acos(ReacPad[0]/vect2))){
                                    comptPadX=comptPadX-1;comptPadY+=1;
                                    EntrancePad[0]=PadDim;EntrancePad[1]=0.;
                                    distEjectil = TMath::Abs((PadDim-ReacPad[1])/sin(thetaEjectil_lab));
                                }
                                //5
                                if(thetaEjectil_lab>(pi_val-acos(ReacPad[0]/vect2)) && thetaEjectil_lab<pi_val+acos(ReacPad[0]/vect3)){
                                    comptPadX=comptPadX-1;
                                    EntrancePad[0]=PadDim;
                                    distEjectil = TMath::Abs((ReacPad[0])/cos(thetaEjectil_lab));
                                    EntrancePad[1] = ReacPad[1]+ distEjectil*sin(thetaEjectil_lab);
                                }
                                //6
                                if(thetaEjectil_lab==pi_val+acos(ReacPad[0]/vect3)){
                                    comptPadX=comptPadX-1;  comptPadY=comptPadY-1;
                                    EntrancePad[1]=PadDim;EntrancePad[0]=PadDim;
                                    distEjectil = TMath::Abs((ReacPad[0])/cos(thetaEjectil_lab));
                                }
                                //7
                                if(thetaEjectil_lab>pi_val+acos(ReacPad[0]/vect3) && thetaEjectil_lab<2*pi_val-acos((PadDim- ReacPad[0])/vect1)){
                                    comptPadY=comptPadY-1;
                                    EntrancePad[1]=PadDim;
                                    distEjectil =TMath::Abs(ReacPad[1]/sin(thetaEjectil_lab));
                                    EntrancePad[0]=ReacPad[0] +distEjectil*cos(thetaEjectil_lab);
                                    
                                }
                                //8
                                if(thetaEjectil_lab==2*pi_val-acos((PadDim- ReacPad[0])/vect1)){
                                    comptPadX+=1;comptPadY=comptPadY-1;
                                    EntrancePad[1]=PadDim;EntrancePad[0]=0;
                                    distEjectil =TMath::Abs(ReacPad[1]/sin(thetaEjectil_lab));
                                }
                                if(SPgas<2){
                                    ElossEjectil = loss_E(SPgas, eEjectil_lab, mEjectil[r]/mu, distEjectil, dxGas, EnergySp[r+3][SPgas+2], SPel[r+3][SPgas+2], SPnu[r+3][SPgas+2]);
                                }
                                if(SPgas>1){
                                    ElossEjectil  = loss_E(0, eEjectil_lab, mEjectil[r]/mu,distEjectil, dxGas, EnergySp[r+3][2],SPel[r+3][2],  SPnu[r+3][2]);
                                    ElossEjectil  = (ElossEjectil  + loss_E(1,eEjectil_lab, mEjectil[r]/mu, distEjectil, dxGas,EnergySp[r+3][3],SPel[r+3][3], SPnu[r+3][3]))/2.0;
                                }
                                Epad[XReac][YReac]= Epad[XReac][YReac]+randEloss->Gaus(ElossEjectil,sigmaPad);
                                eEjectil_lab = eEjectil_lab - ElossEjectil;
                                TrajectoryEjectil->Fill(comptPadX,comptPadY);
                                while(comptPadX<NpadX &&  comptPadY<NpadY &&  eEjectil_lab>0.015 &&  comptPadX>-1 && comptPadY>-1){
                                     cout<<comptPadX << " "<< NpadX<< " "<<comptPadY<<" "<<NpadY<<" "<<eEjectil_lab<<" "<<thetaEjectil_lab<<endl;
                                     ReacPad[0]=EntrancePad[0];ReacPad[1]=EntrancePad[1];
                                     vect=sqrt(pow(PadDim- ReacPad[0],2)+pow(PadDim- ReacPad[1],2));
                                     vect1=sqrt(pow(PadDim- ReacPad[0],2)+pow(ReacPad[1],2));
                                     vect2=sqrt(pow(ReacPad[0],2)+pow(PadDim- ReacPad[1],2));
                                     vect3==sqrt(pow(ReacPad[0],2)+pow(ReacPad[1],2));
                                     //1
                                     if(thetaEjectil_lab<acos((PadDim-ReacPad[0])/vect) || thetaEjectil_lab>2*pi_val-acos((PadDim- ReacPad[0])/vect1)){
                                         comptPadX+=1;
                                         EntrancePad[0]=0.;
                                         distEjectil = TMath::Abs((PadDim-ReacPad[0])/cos(thetaEjectil_lab));
                                         EntrancePad[1]= ReacPad[1]+distEjectil*sin(thetaEjectil_lab);
                                     }
                                     //2
                                     if(thetaEjectil_lab==acos((PadDim-ReacPad[0])/vect)){
                                         comptPadX+=1;comptPadY+=1;
                                         EntrancePad[0]=0.;EntrancePad[1]=0.;
                                         distEjectil = TMath::Abs((PadDim-ReacPad[0])/cos(thetaEjectil_lab));
                                     }
                                     //3
                                     if(thetaEjectil_lab>acos((PadDim-ReacPad[0])/vect) && thetaEjectil_lab<(pi_val-acos(ReacPad[0]/vect2))){
                                         comptPadY+=1;
                                         EntrancePad[1]=0.;
                                         distEjectil =TMath::Abs((PadDim-ReacPad[1])/sin(thetaEjectil_lab));
                                         EntrancePad[0] = ReacPad[0] +distEjectil*cos(thetaEjectil_lab);
                                     }
                                     //4
                                     if(thetaEjectil_lab==(pi_val-acos(ReacPad[0]/vect2))){
                                         comptPadX=comptPadX-1;comptPadY+=1;
                                         EntrancePad[0]=PadDim;EntrancePad[1]=0.;
                                         distEjectil = TMath::Abs((PadDim-ReacPad[1])/sin(thetaEjectil_lab));
                                     }
                                     //5
                                     if(thetaEjectil_lab>(pi_val-acos(ReacPad[0]/vect2)) && thetaEjectil_lab<pi_val+acos(ReacPad[0]/vect3)){
                                         comptPadX=comptPadX-1;
                                         EntrancePad[0]=PadDim;
                                         distEjectil = TMath::Abs((ReacPad[0])/cos(thetaEjectil_lab));
                                         EntrancePad[1] = ReacPad[1]+ distEjectil*sin(thetaEjectil_lab);
                                     }
                                     //6
                                     if(thetaEjectil_lab==pi_val+acos(ReacPad[0]/vect3)){
                                         comptPadX=comptPadX-1;  comptPadY=comptPadY-1;
                                         EntrancePad[1]=PadDim;EntrancePad[0]=PadDim;
                                         distEjectil = TMath::Abs((ReacPad[0])/cos(thetaEjectil_lab));
                                     }
                                     //7
                                     if(thetaEjectil_lab>pi_val+acos(ReacPad[0]/vect3) && thetaEjectil_lab<2*pi_val-acos((PadDim- ReacPad[0])/vect1)){
                                         comptPadY=comptPadY-1;
                                         EntrancePad[1]=PadDim;
                                         distEjectil =TMath::Abs(ReacPad[1]/sin(thetaEjectil_lab));
                                         EntrancePad[0]=ReacPad[0] +distEjectil*cos(thetaEjectil_lab);
                                         
                                     }
                                     //8
                                     if(thetaEjectil_lab==2*pi_val-acos((PadDim- ReacPad[0])/vect1)){
                                         comptPadX+=1;comptPadY=comptPadY-1;
                                         EntrancePad[1]=PadDim;EntrancePad[0]=0;
                                         distEjectil =TMath::Abs(ReacPad[1]/sin(thetaEjectil_lab));
                                     }
                                    if(SPgas<2){
                                        ElossEjectil = loss_E(SPgas,eEjectil_lab, mEjectil[r]/mu, distEjectil, dxGas, EnergySp[r+3][SPgas+2], SPel[r+3][SPgas+2], SPnu[r+3][SPgas+2]);
                                    }
                                    if(SPgas>1){
                                        ElossEjectil  = loss_E(0,eEjectil_lab, mEjectil[r]/mu,distEjectil, dxGas, EnergySp[r+3][2],SPel[r+3][2],  SPnu[r+3][2]);
                                        ElossEjectil  = (ElossEjectil  + loss_E(1,eEjectil_lab, mEjectil[r]/mu, distEjectil, dxGas,EnergySp[r+3][3],SPel[r+3][3], SPnu[r+3][3]))/2.0;
                                    }
                                    Epad[comptPadX][comptPadY]= Epad[comptPadX][comptPadY]+randEloss->Gaus(ElossEjectil,sigmaPad);
                                    eEjectil_lab = eEjectil_lab - ElossEjectil;
                                    TrajectoryEjectil->Fill(comptPadX,comptPadY);
                                    if(plotlive==1){
                                        ctraces->cd(2);
                                        TrajectoryEjectil->Draw("cont4");
                                        ctraces->Update();
                                        gSystem->ProcessEvents();
                                        gSystem->Sleep(break_trajectory);
                                    }
                                }
                            }
                            tree->Fill();
                        }
                    }
                }
            }
        }
        XReac+=1;
    }
    //********************************
    //SAVING TREE FILE
    //********************************
    filesave = pathTrees + Form("tree_pad%ix%i.root",int(PadDim),int(PadDim));
    cout<<"simulation done "<< endl;
    TFile fsave(filesave.c_str(),"create");
    tree->Write();
}
