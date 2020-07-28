void makeDrawDeckLines(){
	gStyle->SetOptStat(0);
	TFile* infile_dat = TFile::Open("degALL_data_2017_newProspectusGraphs/degALL_data_2017_hists_DSelector.root");
	TH2F *any2DHist;
	infile_dat->GetObject("tetaVsMpi0eta",any2DHist);
	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);

	int num_tBins=14;
	double tMin=0;
	double tMax=2.8;
	const int num_massBins=12;
	double mMin=1.6;
	double mMax=2.8;
	double tStep=(tMax-tMin)/num_tBins;
	double mStep=(mMax-mMin)/num_massBins;

	any2DHist->Draw("COLZ");
	cout << "Drew teta vs Mpi0eta" << endl;
	TLine *line = new TLine();	
	line->SetLineColor(kRed);
	double mass;
	double t;
	for (int iMass=0; iMass<num_massBins+1; ++iMass){
		cout << "Drawing mass line" << endl;
		mass = mMin+iMass*mStep;
		line->DrawLine(mass,tMin,mass,tMax);
	}
	for (int it=0; it<num_tBins+1; ++it){
		cout << "Drawing t line" << endl;
		t = tMin+it*tStep;
		line->DrawLine(mMin,t,mMax,t);
	}
	allCanvases->SaveAs("newGraphs/overlayDeck.png");
}
