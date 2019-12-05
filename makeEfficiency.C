void makeEfficiency(){
	TFile* genFile = TFile::Open("v20_flat_gen_hists_DSelector_pi0eta.root");
	TFile* recFile = TFile::Open("pi0eta_flat8GeVPlus_hists_DSelector.root");
	TCanvas *allCanvases_yields = new TCanvas("anyHists_yields","",1440,900);
	allCanvases_yields->Divide(3,4,0,0);

	TH1F *tetaVsMpi0eta_genCounts;
	TH1F *tpi0VsMpi0eta_genCounts;
	TH1F *tetaVsMpi0eta_recCounts;
	TH1F *tpi0VsMpi0eta_recCounts;
	
	genFile->GetObject("tetaVsMpi0eta_genCounts",tetaVsMpi0eta_genCounts);
	genFile->GetObject("tpi0VsMpi0eta_genCounts",tpi0VsMpi0eta_genCounts);
	recFile->GetObject("tetaVsMpi0eta_recCounts",tetaVsMpi0eta_recCounts);
	recFile->GetObject("tpi0VsMpi0eta_recCounts",tpi0VsMpi0eta_recCounts);

	int num_tBins=14;
	int num_massBins=12;
	const int numHists = num_tBins*num_massBins;

	// bin0 = underflow
	// bin1 = -1 which we used to hold all the data that wasn't in a bin
	
	double c_teta_genCounts;
	double c_tpi0_genCounts;
	double c_teta_recCounts;
	double c_tpi0_recCounts;
	double efficiencies_pi0[numHists];
	double efficiencies_eta[numHists];
	int stopIter=0;
	for (int i=2; i<numHists+2; ++i) { 
		if ( stopIter == 0 ) { 
			c_teta_genCounts = tetaVsMpi0eta_genCounts->GetBinContent(i);
			c_tpi0_genCounts = tpi0VsMpi0eta_genCounts->GetBinContent(i);
			c_teta_recCounts = tetaVsMpi0eta_recCounts->GetBinContent(i);
			c_tpi0_recCounts = tpi0VsMpi0eta_recCounts->GetBinContent(i);
			
			if ( c_teta_genCounts == 0 ) { stopIter=i; continue; } 
			if ( c_tpi0_genCounts == 0 ) { stopIter=i; continue; } 
			if ( c_teta_recCounts == 0 ) { stopIter=i; continue; } 
			if ( c_tpi0_recCounts == 0 ) { stopIter=i; continue; } 

			cout << "c_teta_genCounts, c_tpi0_genCounts, c_teta_recCounts, c_tpi0_recCounts: " << c_teta_genCounts << ", " << c_tpi0_genCounts << ", " << c_teta_recCounts << ", " << c_tpi0_recCounts << endl;

			efficiencies_eta[numHists-2] = c_teta_recCounts/c_teta_genCounts; 
			efficiencies_pi0[numHists-2] = c_tpi0_recCounts/c_tpi0_genCounts; 
			cout << efficiencies_eta[numHists-2] << endl;
		}
	}
}
