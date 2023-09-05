#ifndef TOOLS_H
#define TOOLS_H
#include "Topology.h"

double interpol(double x[], double y[], int const n, double p)
{
    int i = 1;
    double val, b1, c1, d1, c2, d2;

    // b1 false if something is wrong
    bool bb1 = true;

    // loop breaking
    bool bb2 = true;

    if (p <= x[1])
    {
        bb1 = false;
    }
    else
    {
        do
        {
            i++;
            if (p <= x[i])bb2 = false;
        } while (bb2 && (i < n - 2));
    }

    if (p > x[n - 2])bb1 = false;

    if (bb1 == true)
    {
        b1 = (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]);
        b1 = b1 * p - b1 * x[i - 2] + y[i - 2];

        c1 = (y[i] - y[i - 2]) / (x[i] - x[i - 2]);
        c1 = c1 * p - c1 * x[i - 2] + y[i - 2];

        d1 = (y[i + 1] - y[i - 2]) / (x[i + 1] - x[i - 2]);
        d1 = d1 * p - d1 * x[i - 2] + y[i - 2];

        c2 = (c1 - b1) / (x[i] - x[i - 1]);
        c2 = c2 * p - c2 * x[i - 1] + b1;

        d2 = (d1 - b1) / (x[i + 1] - x[i - 1]);
        d2 = d2 * p - d2 * x[i - 1] + b1;

        val = (d2 - c2) / (x[i] - x[i - 1]);
        val = val * p - val * x[i] + c2;
    }
    else
    {
        val = 0.;
    }

    return val;

}
double loss_E(int SPref, double  EEntry, double red_mass, double path, double dpath, Double_t EnergyPart[], Double_t SPel_Part[], Double_t SPnu_Part[]){
    double Etemp = EEntry;
    double dpart = 0;
    double ElossEl, ElossNuc, EtempLoss;
    while (dpart <= path*1000) {
        EtempLoss = Etemp *1000.*double(1-SPref) + (Etemp/red_mass)*double(SPref);
        ElossEl =  0.001 * interpol(EnergyPart, SPel_Part, Nsp, EtempLoss) * dpath;
        ElossNuc = 0.001* interpol(EnergyPart, SPnu_Part, Nsp, EtempLoss) * dpath;
        Etemp = Etemp - (ElossEl*double(1-SPref)+double(SPref)*1000.*ElossEl) - ElossNuc*double(1-SPref);
        dpart = dpart + dpath;
    }
    return EEntry - Etemp;
}
double derivation_thick_Eeff(double Ei,double Eloss, double cs1, double cs2){
    double coeffcs= (cs1*cs1+cs2*cs2)/(2*pow(cs1-cs2,2));
    double Eeff = Ei - Eloss + Eloss*((-cs2)/(cs1-cs2)+sqrt(coeffcs));
    return Eeff;
}
#endif
