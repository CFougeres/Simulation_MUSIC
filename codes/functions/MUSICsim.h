#ifndef MUSICSIM.H
#define MUSICSIM.H

#include "main.h"

string filesave;
TRandom* randEx;
TRandom* randEloss;

void reaction_balance()
{
    mBeam = ZAMBeam[1] * mu + ZAMBeam[2] ;
    mGas = ZAMGas[1] * mu + ZAMGas[2];
    int dZ, dA;
    for(int r=0;r<Nrecoil;r++){
        mRecoil[r] = ZAMRecoil[r][1] * mu +ZAMRecoil[r][2];
        mRed[r] = mRecoil[r] /mu;
        dZ = ZAMGas[0]+ZAMBeam[0]-ZAMRecoil[r][0];
        dA = ZAMGas[1]+ZAMBeam[1]-ZAMRecoil[r][1];
        if(dZ==1){
            ZAMultEjectil[r][0]=dZ;
            ZAMultEjectil[r][1]=dA;
            multi_Ejectil[r]=1;
            if(dA==1){  mEjectil[r]=mp; }
            if(dA==2){  mEjectil[r]=m2H;}
        }
        if(dZ==0){
            ZAMultEjectil[r][0]=dZ;
            if(dA>0){
                ZAMultEjectil[r][1]=1;
                mEjectil[r] = mn;
                multi_Ejectil[r]=dA;
            }
            if(dA==0){
                ZAMultEjectil[r][1]=0;
                mEjectil[r] = ExRecoil[r]; multi_Ejectil[r]=1;
            }
        }
        if(dZ==2){
            if(dA==2){ ZAMultEjectil[r][0]=1; ZAMultEjectil[r][1]=1; mEjectil[r]=mp;  multi_Ejectil[r]=2;}
            if(dA>2){
                ZAMultEjectil[r][0]=dZ;
                ZAMultEjectil[r][1]=dA;
                multi_Ejectil[r]=1;
                if(dA==4){ mEjectil[r]=m4He;  }
                if(dA==3){ mEjectil[r]=m3He;   }
            }
        }
        ZAMultEjectil[r][2]= multi_Ejectil[r];
        Qval[r] = mBeam+mGas - mRecoil[r] - multi_Ejectil[r]*mEjectil[r];
        cout<< " Qval A,Z "<<ZAMRecoil[r][0] <<"," <<ZAMRecoil[r][1] <<" = "<<Qval[r] <<" multiEjec "<<multi_Ejectil[r]<<endl;
        masses[r][0] = 0.001 *   mRecoil[r];   masses[r][1] = 0.001 * mEjectil[r]*multi_Ejectil[r];
    }
}
#endif
