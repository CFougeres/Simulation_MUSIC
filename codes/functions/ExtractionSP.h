#ifndef EXTRACTIONSP_H
#define EXTRACTIONSP_H

#include "main.h"

using namespace std;

void extraction_sp(string file, int SPref,Double_t EnergySp[Nsp], Double_t SPel[Nsp], Double_t SPnu[Nsp]){
    ifstream in;
    cout<<file<<endl;
    in.open(file.c_str());
    Int_t nlines = 0;
    while (1) {
        if(SPref==0){
            in >> EnergySp[nlines] >> SPel[nlines] >>  SPnu[nlines] ; // SP in keV/microns, E in keV, Range in Angstrom
        //    cout<<EnergySp[nlines] <<endl;
        }
        if(SPref==1){
            in >> EnergySp[nlines] >> SPel[nlines] ; // SP in MeV/microns, E in MeV/u
           // cout<<EnergySp[nlines] <<endl;
        }
        if (!in.good()) break;
        nlines++;
    }
    in.close();
    return 1;
}

#endif
