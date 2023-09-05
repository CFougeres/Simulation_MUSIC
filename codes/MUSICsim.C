//********************************************
//
// ROOT macro of MUSICsim.cpp
//
//Author: C. Fougères (2023)
//any remark(s) to be sent at chloe.fougeres@gmail.com
//********************************************
string path   = gSystem->pwd();

#include "functions/MUSICsim.h"
using namespace std;
void MUSICsim()
{
    printf("=====================================\n");
    printf("===            MUSICsim           ===\n");
    printf("=====================================\n");
    
    //USER INPUTS
    extraction_inputs();
    sigmaStrip = param_inputs[14][0]/2.38;      pressure = param_inputs[5][0] ;
    SPwindow = param_inputs[4][0] ; SPgas = param_inputs[4][1];
    EbeamEntry = param_inputs[1][0]; Nrecoil = param_inputs[6][0];
    Nreac = param_inputs[11][0];  NEx = param_inputs[12][0];
    dx = param_inputs[13][0];
    ThickTargetYield= param_inputs[17][0];
    if(ThickTargetYield==1){
        for(int cc=0;cc<Nstrip+2;cc++){
            cs[cc] =param_inputs[18][cc];
            errcs[cc] =param_inputs[19][cc];
        }
    }
    incSP= param_inputs[20][0];
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
    tree->Branch("stripReac", &stripReac, "stripReac/I");
    tree->SetBranchAddress("stripReac", &stripReac);
    tree->Branch("ReacKind", &ReacKind,"ReacKind/I");// 0 beam 1==(a,a) 2==(a,p) 3==(a,n), 4==(a,2n), 5==(a,g)
    tree->SetBranchAddress("ReacKind", &ReacKind);
    tree->Branch("Estrip", &Estrip,"Estrip[18]/F");
    tree->SetBranchAddress("Estrip", &Estrip);
    tree->Branch("Istrip", &Istrip, "Istrip[18]/I");
    tree->SetBranchAddress("Istrip", &Istrip);

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
        }
    }
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
        ElossBeam = loss_E(SPgas, Ebeam,mRed[0], MUSIC_dead_layer, dxGas, EnergySp[1][SPgas], SPel[1][SPgas], SPnu[1][SPgas]);
    }
    if(SPgas>1){
        ElossBeam = loss_E(0, Ebeam,mRed[0], MUSIC_dead_layer, dxGas, EnergySp[1][0],SPel[1][0],  SPnu[1][0]);
        ElossBeam = (ElossBeam + loss_E(1,Ebeam,mRed[0], MUSIC_dead_layer, dxGas,EnergySp[1][1],SPel[1][1], SPnu[1][1]))/2.0;
    }
    cout<<"Loss in dead layer "<<ElossBeam<<endl;
    Ebeam = Ebeam - ElossBeam;
    EbeamEntry = Ebeam ;
    //********************************
    //Derivation beam trajectory
    //********************************
    for(int ee=0;ee<Nstrip;ee++){
            for(int j=0;j<Nreac;j++){
                Ebeam = EbeamEntry;
                for(int eee=0;eee<Nstrip;eee++){
                    Etempbeam = Ebeam;
                    if(SPgas<2){
                        ElossBeam = loss_E(SPgas, Ebeam,mRed[0], strip_dim, dxGas, EnergySp[1][SPgas], SPel[1][SPgas], SPnu[1][SPgas]);
                    }
                    if(SPgas>1){
                        ElossBeam = loss_E(0, Ebeam,mRed[0],strip_dim, dxGas, EnergySp[1][0],SPel[1][0],  SPnu[1][0]);
                        ElossBeam = (ElossBeam + loss_E(1,Ebeam,mRed[0],strip_dim, dxGas,EnergySp[1][1],SPel[1][1], SPnu[1][1]))/2.0;
                    }
                    Ebeam = Ebeam - ElossBeam;
                    ELossBeamRef[eee]=ElossBeam;
                    Estrip[eee]=ELossBeamRef[eee];
                    Estrip[eee]=randEloss->Gaus(Estrip[eee],sigmaStrip);
                    Istrip[eee]=eee;
                    if(ee==0 && j==0 &&   ThickTargetYield==0){
                        AnCM[eee] = Istrip[eee]+0.5;
                        ECM[eee] = ((mGas)/(mBeam+mGas))*((Ebeam+Etempbeam)/2.0);
                        ECMSP[eee] = ((mGas)/(mBeam+mGas))*((1.0-incSP)*(Ebeam+Etempbeam)/2.0);
                        erInfECM[eee] = ECM[eee] - ((mGas)/(mBeam+mGas))*(Ebeam);
                        erSupECM[eee] = ((mGas)/(mBeam+mGas))*(Etempbeam)-ECM[eee];
                        cout<<AnCM[eee]<< "Ecm "<<ECM[eee]<<" err-/+ "<<erInfECM[eee] <<" "<<erSupECM[eee]<<" "<<incSP*100<<"% on SP "<< TMath::Abs(ECMSP[eee] - ECM[eee])<<" ebeam lab "<<((Ebeam+Etempbeam)/2.0)<<" range "<< ((mGas)/(mBeam+mGas))*Etempbeam  <<" "<< ((mGas)/(mBeam+mGas))*Ebeam <<endl;
                    }
                    if(ee==0 && j==0 &&   ThickTargetYield==1){
                        AnCM[eee] = Istrip[eee]+0.5;
                        thick_Eeff =derivation_thick_Eeff(Ebeam, ElossBeam, cs[eee], cs[eee+2]);
                        incEloss_thick=derivation_thick_Eeff(Ebeam,ElossBeam*(1.-incSP), cs[eee], cs[eee+2]);
                        incEloss_thick2=derivation_thick_Eeff(Ebeam,ElossBeam, cs[eee]+errcs[eee],  cs[eee+2]-errcs[eee+2]);
                        ECM[eee] = ((mGas)/(mBeam+mGas))*(thick_Eeff);
                        erInfECM[eee] = ECM[eee] - ((mGas)/(mBeam+mGas))*sqrt(pow(incEloss_thick,2)+pow(incEloss_thick2,2));
                        erSupECM[eee] = ((mGas)/(mBeam+mGas))*sqrt(pow(incEloss_thick,2)+pow(incEloss_thick2,2))-ECM[eee];
                        cout<<AnCM[eee]<< "Ecm "<<ECM[eee]<<" errInf "<<erInfECM[eee] <<" errSup "<<erSupECM[eee]<<" Ucs= "<<TMath::Abs(thick_Eeff-incEloss_thick2)<<" "<<incSP*100<<"% on SP "<< TMath::Abs(thick_Eeff-incEloss_thick)<<" ebeam lab "<<thick_Eeff<<" range "<< ((mGas)/(mBeam+mGas))*Etempbeam  <<" "<< ((mGas)/(mBeam+mGas))*Ebeam <<endl;
                    }
                }
                ReacKind=0;
                stripReac=ee;
                tree->Fill();
            }
        }

    //********************************
    //Derivation reaction trajectories
    //********************************
    Ebeam = EbeamEntry;
    cout<<"beam entry anode "<< EbeamEntry<<endl;
        for(int i_strip=0;i_strip<Nstrip;i_strip++){
            cout<<"strip reac ongoing "<<i_strip<<endl;
            positionStrip=0.0;        EtempLossbeam=0.0;    Etempbeam=Ebeam;
            while(positionStrip<strip_dim-dx){
                positionStrip =positionStrip+dx;
                positionMUSIC=positionMUSIC+dx;
                if(SPgas<2){
                    EtempLossbeam = loss_E(SPgas,Etempbeam, mRed[0], positionStrip, dxGas, EnergySp[1][SPgas], SPel[1][SPgas], SPnu[1][SPgas]);
                    ElossBeam = loss_E(SPgas,Ebeam,mRed[0], dx, dxGas, EnergySp[1][SPgas], SPel[1][SPgas], SPnu[1][SPgas]);
                }
                if(SPgas>1){
                    EtempLossbeam = loss_E(0, Etempbeam,mRed[0],positionStrip, dxGas, EnergySp[1][0],SPel[1][0],  SPnu[1][0]);
                    EtempLossbeam = ( EtempLossbeam + loss_E(1,Etempbeam,mRed[0],positionStrip, dxGas,EnergySp[1][1],SPel[1][1], SPnu[1][1]))/2.0;
                    ElossBeam  = loss_E(0, Ebeam,mRed[0],dx, dxGas, EnergySp[1][0],SPel[1][0],  SPnu[1][0]);
                    ElossBeam  = ( ElossBeam  + loss_E(1,Ebeam,mRed[0],dx, dxGas,EnergySp[1][1],SPel[1][1], SPnu[1][1]))/2.0;
                }
                Ebeam = Ebeam - ElossBeam;
                stripReac = i_strip;
                vbeamlab = sqrt(1 - 1 / pow((Ebeam / mBeam + 1), 2));
                Pfais = sqrt(Ebeam * Ebeam + 2. * Ebeam * mBeam);
                beam.SetPxPyPzE(0, 0, Pfais * 0.001, (Ebeam+mBeam) * 0.001);
                Tr = beam + target;
                for(int r=0;r<5;r++){
                    ReacKind=r+1;
                    for(int k=0;k<NEx;k++){
                        Exmax = ECM[i_strip] + Qval[r];
                        Extemp =  0.0 + Exmax*randEx->Rndm(); ;//Exmax*randEx->Rndm();
                        if( Exmax<0.0){ Extemp=0;}
                        masses[r][0] = 0.001 *(mRecoil[r]+ Extemp);
                        event[r].SetDecay(Tr, 2, masses[r]);
                        for(int j=0;j<Nreac;j++){
                            for(int jj=0;jj<Nstrip;jj++){Estrip[jj]=0;}
                            if(i_strip>0){
                               for(int jj=0;jj<i_strip;jj++){
                                    Estrip[jj]= randEloss->Gaus(ELossBeamRef[jj] ,sigmaStrip);
                                    Istrip[jj]=jj;
                                }
                            }
                            event[r].Generate();
                            pRecoil = event[r].GetDecay(0); pEjectil = event[r].GetDecay(1);
                            thetaRecoil_lab = pRecoil->Theta();//Rad
                            vRecoil_lab =pRecoil->Beta();
                            eRecoil_lab=(1 / sqrt(1 - pow(vRecoil_lab, 2)) - 1) * mRecoil[r];
                            remainingDist = 18.*strip_dim - positionMUSIC;
                            thetaExit = acos(MUSIC_half_height/sqrt(MUSIC_half_height*MUSIC_half_height+remainingDist*remainingDist));
                            if(eRecoil_lab>0){
                                if(thetaRecoil_lab<=thetaExit){
                                    distRecoil = (strip_dim-positionStrip)/cos(thetaRecoil_lab);
                                    if(SPgas<2){
                                        ElossRecoil = loss_E(SPgas,eRecoil_lab,mRed[r], distRecoil, dxGas, EnergySp[r+3][SPgas], SPel[r+3][SPgas], SPnu[r+3][SPgas]);
                                    }
                                    if(SPgas>1){
                                        ElossRecoil  = loss_E(0, eRecoil_lab,mRed[r],distRecoil, dxGas, EnergySp[r+3][0],SPel[r+3][0],  SPnu[r+3][0]);
                                        ElossRecoil  = (ElossRecoil  + loss_E(1,eRecoil_lab,mRed[r], distRecoil, dxGas,EnergySp[r+3][1],SPel[r+3][1], SPnu[r+3][1]))/2.0;
                                    }
                                    
                                    Estrip[i_strip] = ElossRecoil + EtempLossbeam;
                                    Estrip[i_strip]= randEloss->Gaus(Estrip[i_strip] ,sigmaStrip);
                                    eRecoil_lab = eRecoil_lab - ElossRecoil;
                                    Istrip[i_strip]=i_strip;
                                    for(int ii=i_strip+1;ii<Nstrip;ii++){
                                    distRecoil = strip_dim/cos(thetaRecoil_lab);
                                        if(SPgas<2){
                                            ElossRecoil = loss_E(SPgas,eRecoil_lab,mRed[r], distRecoil, dxGas, EnergySp[r+3][SPgas], SPel[r+3][SPgas], SPnu[r+3][SPgas]);
                                        }
                                        if(SPgas>1){
                                            ElossRecoil  = loss_E(0, eRecoil_lab,mRed[r],distRecoil, dxGas, EnergySp[r+3][0],SPel[r+3][0],  SPnu[r+3][0]);
                                            ElossRecoil  = (ElossRecoil  + loss_E(1,eRecoil_lab,mRed[r], distRecoil, dxGas,EnergySp[r+3][1],SPel[r+3][1], SPnu[r+3][1]))/2.0;
                                        }
                                        Estrip[ii]=ElossRecoil;
                                        Estrip[ii]=randEloss->Gaus(Estrip[ii] ,sigmaStrip);
                                        Istrip[ii]=ii;
                                        eRecoil_lab = eRecoil_lab - ElossRecoil;
                                    }
                                        tree->Fill();
                                }
                                if(thetaRecoil_lab>thetaExit){
                                   //cout<<"escape "<< r <<endl;
                                    distRecoil = (MUSIC_half_height)/tan(thetaRecoil_lab);
                                    if(distRecoil<=strip_dim-positionStrip){
                                        if(SPgas<2){
                                            ElossRecoil = loss_E(SPgas,eRecoil_lab,mRed[r], (strip_dim-positionStrip)/cos(thetaRecoil_lab), dxGas, EnergySp[r+3][SPgas], SPel[r+3][SPgas], SPnu[r+3][SPgas]);
                                        }
                                        if(SPgas>1){
                                            ElossRecoil  = loss_E(0, eRecoil_lab,mRed[r],(strip_dim-positionStrip)/cos(thetaRecoil_lab), dxGas, EnergySp[r+3][0],SPel[r+3][0],  SPnu[r+3][0]);
                                            ElossRecoil  = (ElossRecoil  + loss_E(1,eRecoil_lab,mRed[r], (strip_dim-positionStrip)/cos(thetaRecoil_lab), dxGas,EnergySp[r+3][1],SPel[r+3][1], SPnu[r+3][1]))/2.0;
                                        }
                                        Estrip[i_strip]=ElossRecoil+EtempLossbeam;
                                        Estrip[i_strip]= randEloss->Gaus(Estrip[i_strip],sigmaStrip);
                                        Istrip[i_strip]=i_strip;
                                        if(i_strip<Nstrip-1){for(int ii=i_strip+1;ii<Nstrip;ii++){Estrip[ii]=0; Istrip[ii]=ii;}}
                                            tree->Fill();
                                    }
                                    if(distRecoil>strip_dim-positionStrip){
                                        remainingStrip = int(distRecoil/strip_dim);
                                        if(i_strip<Nstrip-1){
                                            distRecoil = (strip_dim-positionStrip)/cos(thetaRecoil_lab);
                                            if(SPgas<2){
                                                ElossRecoil = loss_E(SPgas,eRecoil_lab,mRed[r],distRecoil, dxGas, EnergySp[r+3][SPgas], SPel[r+3][SPgas], SPnu[r+3][SPgas]);
                                            }
                                            if(SPgas>1){
                                                ElossRecoil  = loss_E(0, eRecoil_lab,mRed[r],distRecoil, dxGas, EnergySp[r+3][0],SPel[r+3][0],  SPnu[r+3][0]);
                                                ElossRecoil  = (ElossRecoil  + loss_E(1,eRecoil_lab,mRed[r],distRecoil, dxGas,EnergySp[r+3][1],SPel[r+3][1], SPnu[r+3][1]))/2.0;
                                            }
                                            Estrip[i_strip] = ElossRecoil + EtempLossbeam;
                                            Estrip[i_strip]= randEloss->Gaus(Estrip[i_strip] ,sigmaStrip);
                                            eRecoil_lab = eRecoil_lab - ElossRecoil;
                                            Istrip[i_strip]=i_strip;
                                            for(int ii=i_strip+1;ii<remainingStrip;ii++){
                                                if(SPgas<2){
                                                    ElossRecoil = loss_E(SPgas,eRecoil_lab,mRed[r],strip_dim/cos(thetaRecoil_lab), dxGas, EnergySp[r+3][SPgas], SPel[r+3][SPgas], SPnu[r+3][SPgas]);
                                                }
                                                if(SPgas>1){
                                                    ElossRecoil  = loss_E(0, eRecoil_lab,mRed[r],strip_dim/cos(thetaRecoil_lab), dxGas, EnergySp[r+3][0],SPel[r+3][0],  SPnu[r+3][0]);
                                                    ElossRecoil  = (ElossRecoil  + loss_E(1,eRecoil_lab,mRed[r],strip_dim/cos(thetaRecoil_lab), dxGas,EnergySp[r+3][1],SPel[r+3][1], SPnu[r+3][1]))/2.0;
                                                }
                                                Estrip[ii]=ElossRecoil;
                                                Estrip[ii]= randEloss->Gaus(Estrip[ii],sigmaStrip);
                                                Istrip[ii]=ii;
                                                eRecoil_lab = eRecoil_lab - ElossRecoil;
                                            }
                                            if(SPgas<2){
                                                ElossRecoil = loss_E(SPgas,eRecoil_lab,mRed[r],(distRecoil-remainingStrip*strip_dim)/cos(thetaRecoil_lab), dxGas, EnergySp[r+3][SPgas], SPel[r+3][SPgas], SPnu[r+3][SPgas]);
                                            }
                                            if(SPgas>1){
                                                ElossRecoil  = loss_E(0, eRecoil_lab,mRed[r],(distRecoil-remainingStrip*strip_dim)/cos(thetaRecoil_lab), dxGas, EnergySp[r+3][0],SPel[r+3][0],  SPnu[r+3][0]);
                                                ElossRecoil  = (ElossRecoil  + loss_E(1,eRecoil_lab,mRed[r],(distRecoil-remainingStrip*strip_dim)/cos(thetaRecoil_lab), dxGas,EnergySp[r+3][1],SPel[r+3][1], SPnu[r+3][1]))/2.0;
                                            }
                                            
                                            Estrip[remainingStrip+1]=ElossRecoil;
                                            Estrip[remainingStrip+1]= randEloss->Gaus(Estrip[remainingStrip+1],sigmaStrip);
                                            Istrip[remainingStrip+1]=remainingStrip+1;
                                            if(remainingStrip+1==Nstrip-1){
                                                    tree->Fill();
                                            if(remainingStrip+1<Nstrip-1){
                                                for(int ii=remainingStrip+2;ii<Nstrip;ii++){ Estrip[ii]=0; Istrip[ii]=ii;}
                                                    tree->Fill();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //********************************
    //SAVING TREE FILE
    //********************************
    filesave = pathTrees + Form("tree_%i%iBeam_%i%iGas.root",int(ZAMBeam[1]),int(ZAMBeam[0]),int(ZAMGas[1]),int(ZAMGas[0]));
    cout<<"simulation traces done "<< endl;
    TFile fsave(filesave.c_str(),"create");
    tree->Write();
    
}
