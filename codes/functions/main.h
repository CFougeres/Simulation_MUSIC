#ifndef MAIN_H
#define MAIN_H

#include "Topology.h"
#include "Tools.h"
#include "ExtractionSP.h"
#include <iostream>
#include <stdlib.h>

Double_t index_inputs[500]={0};
Double_t param_inputs[500][100]={0};

void parseFloats(const string& s , vector<double>& res) {
    istringstream isf(s);
    int comp=0;
    while ( isf.good() ) {
        if(comp==0 || comp>1){
            double value;
            isf >> value;
            res.push_back(value);
        }
        if(comp==1){
            char value[200];
            isf >> value;
        }
        comp++;
    }
}
void extraction_inputs(){
    cout<<path <<endl;    string path_input = path  +"/inputs.dat";     cout<<path_input<<endl;
    char variable[500][200]={0};
    vector<vector<double>> list;
    vector<double> list_temp;
    ifstream in;
    in.open(path_input.c_str());
    Int_t nlines = 0;   string line; int nlist = 0;
    while(1){
        if(nlines==0){
            in >> index_inputs[nlines] >> variable[nlines] ;
            cout<<index_inputs[nlines]<<" "<<variable[nlines] <<endl;
        }
       if(nlines==7 || nlines ==8 || nlines ==9 || nlines ==10 || nlines ==22){
           ws(in);
           getline(in,line);
           parseFloats(line,list_temp);
           list.push_back(list_temp);
           nlist++;
           list_temp.clear();
        }
        if(nlines==4){
            in >> index_inputs[nlines] >> variable[nlines] >>   param_inputs[nlines][0]  >>  param_inputs[nlines][1] >>  param_inputs[nlines][2];
            cout<<index_inputs[nlines]<<" "<<variable[nlines] <<" "<<  param_inputs[nlines][0]<<" "<<  param_inputs[nlines][1]<<" "<< param_inputs[nlines][2] <<endl;
        }
        if(nlines==2 || nlines==3 ){
            in >> index_inputs[nlines] >> variable[nlines] >>  param_inputs[nlines][0]  >>  param_inputs[nlines][1] >>  param_inputs[nlines][2] ;
            cout<<index_inputs[nlines]<<" "<<variable[nlines] <<" "<<  param_inputs[nlines][0]<<" "<< param_inputs[nlines][1] <<" "<<  param_inputs[nlines][2] <<endl;
        }
        if(nlines!=0 && nlines!=2 && nlines!=3 && nlines!=4 && nlines!=7 && nlines!=8 && nlines!=9 && nlines!=10 && nlines!=18 && nlines!=19 && nlines!=22 && nlines!=24){
            in >> index_inputs[nlines] >> variable[nlines] >>  param_inputs[nlines][0];
            cout<<index_inputs[nlines]<<" "<<variable[nlines] <<" "<<  param_inputs[nlines][0] <<endl;
        }
        if(nlines==18 || nlines==19){
            in >> index_inputs[nlines] >> variable[nlines] >>  param_inputs[nlines][0]  >>  param_inputs[nlines][1] >>  param_inputs[nlines][2] >>param_inputs[nlines][3]  >>  param_inputs[nlines][4] >>  param_inputs[nlines][5] >>param_inputs[nlines][6]  >>  param_inputs[nlines][7] >>  param_inputs[nlines][8]>> param_inputs[nlines][9]  >>  param_inputs[nlines][10] >>  param_inputs[nlines][11]>>param_inputs[nlines][12]  >>  param_inputs[nlines][13] >>  param_inputs[nlines][14] >> param_inputs[nlines][15]  >>  param_inputs[nlines][16] >>  param_inputs[nlines][17] >> param_inputs[nlines][18]  >>  param_inputs[nlines][19] ;
            cout<<index_inputs[nlines]<<" "<<variable[nlines] <<" "<<  param_inputs[nlines][0]<<" "<< param_inputs[nlines][1] <<" "<<  param_inputs[nlines][2] <<endl;
        }
        if(nlines==24){
            in >> index_inputs[nlines] >> variable[nlines] >>  param_inputs[nlines][0]  >>  param_inputs[nlines][1];
            cout<<index_inputs[nlines]<<" "<<variable[nlines] <<" "<<  param_inputs[nlines][0] <<" "<< param_inputs[nlines][1] <<endl;
        }
        if (!in.good()) break;
        nlines++;
    }
    in.close();
    cout<<"\n Recoils" <<endl;
    for(int r=0;r<param_inputs[6][0];r++){
        param_inputs[7][r]= int(list[0][r+1]);
        param_inputs[8][r]=int(list[1][r+1]);
        param_inputs[9][r]=list[2][r+1];
        param_inputs[10][r]=list[3][r+1];
        cout<< param_inputs[7][r]<<" "<<param_inputs[8][r]<<" "<<param_inputs[9][r]<<" "<<param_inputs[10][r]<<endl;
    }
    for(int r=0;r<param_inputs[21][0];r++){
        param_inputs[22][r]= int(list[4][r+1]);
        cout<< param_inputs[22][r]<<endl;
    }
    return 1;
}

#endif
