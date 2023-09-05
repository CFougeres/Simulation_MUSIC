#include <iostream>
void distribution(){
	int Nion = 1600;
	double thetai = 0;
	TH1F* theta_distri = new TH1F("theta_distri", Form("P= %.2f Torr; #theta_{beam} (deg)",thetai),45*4,0,45);
	double theta,x,y,z, i; double Pi = TMath::Pi();
	Double_t X[Nion], Y[Nion], Z[Nion];
	std::ifstream in;
	in.open("/Users/chloefougeres/Documents/post_doc/argonne/work/experiment/MUSIC@ATLAS/88Sr_disc_time/sim/scattering/Ti/range.dat");
	Int_t nlines = 0;
	while (1) {
		in >> i >> z >> x >> y;
		Z[nlines] = z / 10000000.; X[nlines] = x / 10000000.;  Y[nlines] = y / 10000000.;
		if (!in.good()) break;
		theta = acos(Z[nlines] / sqrt(X[nlines] * X[nlines] + Y[nlines] * Y[nlines] + Z[nlines] * Z[nlines])) * 180. / Pi;
		theta_distri->Fill(theta);
		nlines++;
	}
	printf(" found %d points\n", nlines);
	in.close();
	TGraph* g = new TGraph(Nion, X, Y);
	g->SetTitle(Form("SRIM scattering P= %.2f Torr;X (mm); Y (mm))",thetai));
	TCanvas* c = new TCanvas("c", "", 800, 800); c->Divide(2, 1);
	c->cd(1);	g->Draw("A*");
	c->cd(2);	theta_distri->Draw();
}
