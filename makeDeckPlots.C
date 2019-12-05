int numDOFsig=8;
Double_t g2(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	Double_t r2 = Double_t((x[1]-par[3])/par[4]);
	return par[5]+par[6]*(x[0]-par[1])+par[7]*(x[1]-par[3])+par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}

//int numDOF_1D=9;
//Double_t g1_eta(Double_t *x, Double_t *par) {
//	Double_t r1 = Double_t((par[8]-par[1])/par[2]);
//	Double_t r2 = Double_t((x[0]-par[3])/par[4]);
//	return par[5]+par[6]*(par[8]-par[1])+par[7]*(x[0]-par[3])+par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
//}
//Double_t g1_pi0(Double_t *x, Double_t *par) {
//	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
//	Double_t r2 = Double_t((par[8]-par[3])/par[4]);
//	return par[5]+par[6]*(x[0]-par[1])+par[7]*(par[8]-par[3])+par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
//}

int numDOF_1D=11;
Double_t g1_pi0(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	//return par[5]*par[8]+par[8]*par[6]*(x[0]-par[1])+par[7]*par[8]*(par[8]/2-par[3])+par[0]*par[8]/sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-0.5*(r1*r1));
	return ( par[5]+par[6]*(x[0]-par[1])-par[7]*par[3])*(par[10]-par[9])*par[8]+par[7]*(0.5*par[8]*(par[10]*par[10]-par[9]*par[9]))+par[0]*par[8]/sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-0.5*(r1*r1));
}

Double_t g1_eta(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[3])/par[4]);
	//return par[5]*par[8]+par[8]*par[7]*(x[0]-par[3])+par[6]*par[8]*(par[8]/2-par[1])+par[0]*par[8]/sqrt(2*TMath::Pi())/par[4]*TMath::Exp(-0.5*(r1*r1));
	return ( par[5]+par[7]*(x[0]-par[3])-par[6]*par[1])*(par[10]-par[9])*par[8]+par[6]*(0.5*par[8]*(par[10]*par[10]-par[9]*par[9]))+par[0]*par[8]/sqrt(2*TMath::Pi())/par[4]*TMath::Exp(-0.5*(r1*r1));
}


int numDOFsigFlat=6;
Double_t gflat2(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	Double_t r2 = Double_t((x[1]-par[3])/par[4]);
	return par[5]+par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}

int numDOFsig1D=5;
Double_t gaus(Double_t *x, Double_t *par){
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	return par[3]+par[4]*(x[0]-par[1])+par[0]/sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-0.5*r1*r1 );
}


void makeDeckPlot(string selectionString){
   	TFile *deckDiagnosticFile = new TFile("deckDiagnosticHists.root", "RECREATE");
	// *********** PREPARING THE CODE ***************
	// *********************************************
	//gStyle->SetErrorX(0.000001); // remove the x-error bars
    	ofstream logFile;
    	logFile.open(("deckPlots/"+selectionString+"/failedFittingPlotIDs.txt").c_str());

	
	TFile* dataFile = TFile::Open("pi0eta_datatreeFlat_DSelector.root");
	TFile* dataHists = TFile::Open("pi0eta_data_hists_DSelector.root");
	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	TCanvas *allCanvases_tSlope = new TCanvas("anyHists_tSlope","Blue=Eta Red=Pi0 DarkGray=teta Gray=tpi0",1440,900);
	TCanvas *allCanvases_yields = new TCanvas("anyHists_yields","Blue=Eta Red=Pi0 DarkGray=teta Gray=tpi0",1440,900);
	TCanvas *allCanvases_unscaledYields = new TCanvas("anyHists_unscaledYields","Blue=Eta Red=Pi0 DarkGray=teta Gray=tpi0",1440,900);

	allCanvases_tSlope->SetLeftMargin(0.15);
	allCanvases_yields->SetLeftMargin(0.15);
	allCanvases_unscaledYields->SetLeftMargin(0.15);
	allCanvases_tSlope->SetBottomMargin(0.2);
	allCanvases_yields->SetBottomMargin(0.2);
	allCanvases_unscaledYields->SetBottomMargin(0.2);

	allCanvases_yields->Divide(3,4,0,0);
	allCanvases_unscaledYields->Divide(3,4,0,0);
	allCanvases_tSlope->Divide(3,4,0,0);

	TTree *dataTree;
	dataFile->GetObject("pi0eta_datatree_flat",dataTree);

	// THese will define the bins
	int num_tBins=14;
	double tMin=0;
	double tMax=2.8;
	int num_massBins=12;
	double mMin=1.7;
	double mMax=2.9;
	double tStep=(tMax-tMin)/num_tBins;
	double mStep=(mMax-mMin)/num_massBins;
	const int numHists = (const int)num_tBins*num_massBins;
	cout << "numHists: " << numHists << endl;
	
	// Match the addresses from the flat tree
	double new_t;
	double new_m;
	double new_accWeight;
	Bool_t new_isUnique12B_1234B;
	Bool_t new_isUnique34B_1234B;
	Bool_t isUnique;
	double new_meta;
	double new_mpi0;
	Bool_t new_ptGT1;
	Bool_t new_ptLT05;
	Bool_t new_ptGT05LT1;
        dataTree->SetBranchAddress("AccWeight",&new_accWeight);
        dataTree->SetBranchAddress("isNotRepeated_pi0_pi0eta",&new_isUnique12B_1234B);
        dataTree->SetBranchAddress("isNotRepeated_eta_pi0eta",&new_isUnique34B_1234B);
	dataTree->SetBranchAddress("Meta_meas",&new_meta);
	dataTree->SetBranchAddress("Mpi0_meas",&new_mpi0);
	dataTree->SetBranchAddress("Mpi0eta_meas",&new_m);
	dataTree->SetBranchAddress("ptGT1",&new_ptGT1);
	dataTree->SetBranchAddress("ptLT05",&new_ptLT05);
	dataTree->SetBranchAddress("ptGT05LT1",&new_ptGT05LT1);
	cout << "Loaded addresses" << endl;


	// *********** CALC EFFIENCY FIRST ***************
	// *********************************************
	 // For the efficiency plots
	TFile* genFile = TFile::Open("flatUpTo3GeVResMass_gen_hists_DSelector_pi0eta.root");
	TFile* recFile = TFile::Open("pi0eta_flat8GeVPlus_hists_DSelector.root");
	TH1F *tetaVsMpi0eta_genCounts;
	TH1F *tpi0VsMpi0eta_genCounts;
	TH1F *tetaVsMpi0eta_recCounts;
	TH1F *tpi0VsMpi0eta_recCounts;
	genFile->GetObject( ("tetaVsMpi0eta_genCounts_"+selectionString).c_str(),tetaVsMpi0eta_genCounts);
	genFile->GetObject( ("tpi0VsMpi0eta_genCounts_"+selectionString).c_str(),tpi0VsMpi0eta_genCounts);
	recFile->GetObject( ("tetaVsMpi0eta_recCounts_"+selectionString).c_str(),tetaVsMpi0eta_recCounts);
	recFile->GetObject( ("tpi0VsMpi0eta_recCounts_"+selectionString).c_str(),tpi0VsMpi0eta_recCounts);
	allCanvases->cd(); allCanvases->Clear(); allCanvases->SetLogy(); 
	tetaVsMpi0eta_genCounts->Draw();
	allCanvases->SaveAs(("deckPlots/"+selectionString+"/tetaVsMpi0eta_genCounts.png").c_str());
	allCanvases->Clear();allCanvases->SetLogy();
	tpi0VsMpi0eta_genCounts->Draw();
	allCanvases->SaveAs(("deckPlots/"+selectionString+"/tpi0VsMpi0eta_genCounts.png").c_str());
	allCanvases->Clear();allCanvases->SetLogy();
	tetaVsMpi0eta_recCounts->Draw();
	allCanvases->SaveAs(("deckPlots/"+selectionString+"/tetaVsMpi0eta_recCounts.png").c_str());
	allCanvases->Clear();allCanvases->SetLogy();
	tpi0VsMpi0eta_recCounts->Draw();
	allCanvases->SaveAs(("deckPlots/"+selectionString+"/tpi0VsMpi0eta_recCounts.png").c_str());

	allCanvases->SetLogy(0);

	TH1F *hist_efficiencies_pi0[num_massBins];
	TH1F *hist_efficiencies_eta[num_massBins];
	for (iHist=0; iHist < num_massBins; ++iHist){
		hist_efficiencies_pi0[iHist] = new TH1F("","",num_tBins,tMin,tMax);
		hist_efficiencies_eta[iHist] = new TH1F("","",num_tBins,tMin,tMax);
	}
	double c_teta_genCounts;
	double c_tpi0_genCounts;
	double c_teta_recCounts;
	double c_tpi0_recCounts;
	std::vector<double > efficiencies_pi0; efficiencies_pi0.reserve(numHists);
	std::vector<double > efficiencies_eta; efficiencies_eta.reserve(numHists);
	std::vector< std::vector<double> > efficiencies; 
	std::vector<double > efficiencies_error_pi0; efficiencies_error_pi0.reserve(numHists);
	std::vector<double > efficiencies_error_eta; efficiencies_error_eta.reserve(numHists);
	std::vector< std::vector<double> > efficiencies_error; 
	double maxEfficiency=DBL_MIN;
	bool skipCalc;
	// bin0 = underflow
	// bin1 = -1 which we used to hold all the data that wasn't in a bin
	for (int i=2; i<numHists+2; ++i) { 
		// To make the code more understandable we will just calculate the efficiencies first then fill the histogram. There will be too much going on with the array indicies if
		//     we have to consider the overflow bin and the shifting of the histograms we read in since they also have an overflow and an extra bin that holds non-bin events. 
		//     ALSO: we have to find the maximum value first
		int j=i-2; // shifted i
		c_teta_genCounts = tetaVsMpi0eta_genCounts->GetBinContent(i);
		c_tpi0_genCounts = tpi0VsMpi0eta_genCounts->GetBinContent(i);
		c_teta_recCounts = tetaVsMpi0eta_recCounts->GetBinContent(i);
		c_tpi0_recCounts = tpi0VsMpi0eta_recCounts->GetBinContent(i);
		
		skipCalc=false;
		cout << " ** IF COUNTS <=0 THEN WE SET THEM =0 **" << endl;
		if ( c_teta_genCounts <= 0 ) { efficiencies_eta[j]=0; efficiencies_error_eta[j]=0; skipCalc=true; cout << "SETTNG TO ZERO TO c_teta_genCounts ZERO" << endl;} 
		if ( c_tpi0_genCounts <= 0 ) { efficiencies_pi0[j]=0; efficiencies_error_pi0[j]=0; skipCalc=true; cout << "SETTNG TO ZERO TO c_tpi0_genCounts ZERO" << endl;} 
		if ( c_teta_recCounts <= 0 ) { efficiencies_eta[j]=0; efficiencies_error_eta[j]=0; skipCalc=true; cout << "SETTNG TO ZERO TO c_teta_recCounts ZERO" << endl;} 
		if ( c_tpi0_recCounts <= 0 ) { efficiencies_pi0[j]=0; efficiencies_error_pi0[j]=0; skipCalc=true; cout << "SETTNG TO ZERO TO c_tpi0_recCounts ZERO" << endl;} 

		cout << "c_teta_genCounts, c_tpi0_genCounts, c_teta_recCounts, c_tpi0_recCounts: " << c_teta_genCounts << ", " << c_tpi0_genCounts << ", " << c_teta_recCounts << ", " << c_tpi0_recCounts << endl;
		if (!skipCalc) {
			efficiencies_eta.push_back( c_teta_recCounts/c_teta_genCounts); 
			efficiencies_pi0.push_back( c_tpi0_recCounts/c_tpi0_genCounts); 
			efficiencies_error_eta.push_back( efficiencies_eta[j]*sqrt(1/c_teta_recCounts+1/c_teta_genCounts) ); 
			efficiencies_error_pi0.push_back( efficiencies_pi0[j]*sqrt(1/c_tpi0_recCounts+1/c_tpi0_genCounts) ); 

			if ( efficiencies_eta[j] > maxEfficiency ) { maxEfficiency = efficiencies_eta[j]; }
			if ( efficiencies_pi0[j] > maxEfficiency ) { maxEfficiency = efficiencies_pi0[j]; }
		}
		
		int massBin = j/num_tBins;
		int tBin = j%num_tBins;
		cout << "\tETA - Filling Mass Bin: " << massBin << " at t Bin: " << tBin << " with efficiency,error: " << efficiencies_eta[j] << "," << efficiencies_error_eta[j] << endl;
		cout << "\tPI0 -Filling Mass Bin: " << massBin << " at t Bin: " << tBin << " with efficiency,error: " << efficiencies_pi0[j] << "," << efficiencies_error_pi0[j] << endl;
		hist_efficiencies_pi0[massBin]->SetBinContent( tBin+1, efficiencies_pi0[j]);
		hist_efficiencies_pi0[massBin]->SetBinError( tBin+1, efficiencies_error_pi0[j]);
		hist_efficiencies_eta[massBin]->SetBinContent( tBin+1,  efficiencies_eta[j]);
		hist_efficiencies_eta[massBin]->SetBinError( tBin+1, efficiencies_error_eta[j]);
	}
	efficiencies.push_back(efficiencies_eta);
	efficiencies.push_back(efficiencies_pi0);
	efficiencies_error.push_back(efficiencies_error_eta);
	efficiencies_error.push_back(efficiencies_error_pi0);

	//for (int massBin=0; massBin < num_massBins; ++massBin){
	//	allCanvases_yields->cd(massBin+1);
	//	hist_efficiencies_pi0[massBin]->SetMarkerStyle(kFullCircle);
	//	hist_efficiencies_pi0[massBin]->SetMarkerSize(0.5);
	//	hist_efficiencies_pi0[massBin]->SetMarkerColor(kRed);
	//	hist_efficiencies_pi0[massBin]->Draw("E1 PMC");
	//	hist_efficiencies_eta[massBin]->SetMarkerStyle(kFullCircle);
	//	hist_efficiencies_eta[massBin]->SetMarkerSize(0.5);
	//	hist_efficiencies_eta[massBin]->SetMarkerColor(kBlue);
	//	hist_efficiencies_eta[massBin]->Draw("E1 SAME");
	//}
	//allCanvases_yields->SaveAs("effs.png");
	

	// -----------------------PROBABLY WONT NEED THIS ANYMORE -------------------------------
	// -------------------------------------------------------------------------------------
	// *********** LOAD THE teta/tpi0 HISTOGRAMS TO OVERLAY IN THE END  ***************
	// *********************************************
	//allCanvases->SetLogy(0);
	//TH1F *teta_binnedHists[num_massBins];
	//TH1F *tpi0_binnedHists[num_massBins];
	//TF1 *linFit_eta = new TF1("linFit_eta","pol2",tMin,tMax);
	//TF1 *linFit_pi0 = new TF1("linFit_pi0","pol2",tMin,tMax);
	//double maxT=DBL_MIN;
	//for (int iMass=0; iMass<num_massBins; ++iMass){
	//	dataHists->GetObject(("tetaMassBinned"+to_string(iMass)).c_str(),teta_binnedHists[iMass]);
	//	dataHists->GetObject(("tpi0MassBinned"+to_string(iMass)).c_str(),tpi0_binnedHists[iMass]);
	//	if ( maxT < teta_binnedHists[iMass]->GetMaximum() ) { 
	//		maxT = teta_binnedHists[iMass]->GetMaximum(); 
	//	}
	//	if ( maxT < tpi0_binnedHists[iMass]->GetMaximum() ) { 
	//		maxT = tpi0_binnedHists[iMass]->GetMaximum(); 
	//	}
	//	teta_binnedHists[iMass]->Divide(hist_efficiencies_eta[iMass]);
	//	teta_binnedHists[iMass]->Fit("linFit_eta","Q");
	//	allCanvases->SaveAs(("deckPlots/teta_binnedHists_"+to_string(iMass)+".png").c_str());
	//	allCanvases->Clear();
	//	tpi0_binnedHists[iMass]->Divide(hist_efficiencies_pi0[iMass]);
	//	tpi0_binnedHists[iMass]->Fit("linFit_pi0","Q");
	//	tpi0_binnedHists[iMass]->Draw();
	//	allCanvases->SaveAs(("deckPlots/tpi0_binnedHists_"+to_string(iMass)+".png").c_str());
	//}



	// *********** CALCULATE YIELDS IN BINS  ***************
	// *********************************************
	//string branchNames[2]={"mandelstam_tpi0_meas","mandelstam_teta_meas"};
	string branchNames[2]={"mandelstam_teta_meas","mandelstam_tpi0_meas"};
	int branchIdx=-1;
	std::vector<TH1F> massHists_eta;
	std::vector<TH1F> massHists_pi0;
	std::vector<TF1> expPlusLinFits_eta;
	std::vector<TF1> expPlusLinFits_pi0;
	std::vector<double> tSlopes_eta;
	std::vector<double> tSlopes_pi0;
	std::vector<double> tSlopesError_eta;
	std::vector<double> tSlopesError_pi0;

	double maxYield=DBL_MIN;
	double minYield=DBL_MAX;
	double maxUnscaledYield=DBL_MIN;
	double yields[2][numHists];
	double yieldErrors[2][numHists];
	double unscaledYields[2][numHists];
	double unscaledYieldErrors[2][numHists];
	for ( string branchName: branchNames ) {
		++branchIdx;
		gStyle->SetErrorX(0.000001); // remove the x-error bars
		gStyle->SetTitleSize(0.06,"t");
		logFile << branchName << endl;
		allCanvases->cd();
		dataTree->SetBranchAddress( branchName.c_str(),&new_t);

		// So the following indicies would describe the bin for a specific variable. We have to multiplex them into a single array. 
		// The obvious way to do that is to use the array index = num_tBins*m+t
		// This would order the array where the first [0,num_tBins] would belong in the smallest mass bin, [num_tBins,2*num_tBins] belongs to the second smallest.
		int idx_t;
		int idx_m;
		TH2F* hists_meta_mpi0[numHists];
		//TH1F* hists_meta[numHists];
		//TH1F* hists_mpi0[numHists];
		TH2F* hists_mpi0eta_t[numHists];
		TH2F* full_meta_mpi0 = new TH2F(("full_meta_mpi0_"+std::to_string(branchIdx)).c_str(), "Cuts=GeneralCuts;M(#pi_{0}) GeV;M(#eta)(GeV)", 100,0.05,0.25,100,0.25,0.85 );

		double xfitMin=0.085;
		double xfitMax=0.185;

		double yfitMin;
		double yfitMax;

		if (selectionString=="tLT1" || selectionString=="tAll" || selectionString=="tGT05LT1"){
			//  These were fitRange we used for the analysis presentation
			yfitMin=0.4;
			yfitMax=0.64;
		}
		else if (selectionString=="tLT05") { 
			//  Since tLT05 doesnt fit easily we can shorten the range in the eta
			yfitMin=0.46;
			yfitMax=0.625;
		}

		double chiSq;
		double dof;

		double yBins[2] = {0.25,0.85};
		double xBins[2] = {0.05,0.25};
		int nxBins=100;
		int nyBins=100;
		double xWidth = (xBins[1]-xBins[0])/nxBins;
		double yWidth = (yBins[1]-yBins[0])/nyBins;
		double differentialArea = 1/(xWidth*yWidth);
		cout << "Differential Area:"<< differentialArea << endl;
		TH1F* chiSqPerDOF = new TH1F(("chiSqPerDOF_"+std::to_string(branchIdx)).c_str(), "#chi^{2}/n_{dof}; #chi^{2}/n_{dof}; Events / 0.1", 30,0,3 );
		TH1F* full_meta = new TH1F(("full_meta_"+std::to_string(branchIdx)).c_str(), ";M(#eta)(GeV)", nxBins,yBins[1],yBins[0] );
		TH1F* full_mpi0 = new TH1F(("full_mpi0_"+std::to_string(branchIdx)).c_str(), ";M(#pi_0)(GeV)",nyBins,xBins[1],xBins[0] );
		
		for (int i=0; i< numHists; ++i){
			hists_meta_mpi0[i] = new TH2F(("meta_mpi0"+to_string(i)+"_"+std::to_string(branchIdx)).c_str(), "Cuts=GeneralCuts;M(#pi_{0}) GeV;M(#eta)(GeV)", 100,0.05,0.25,100,0.25,0.85 ); 
			//hists_meta[i] = new TH1F(("meta"+to_string(i)+"_"+std::to_string(branchIdx)).c_str(), "Cuts=GeneralCuts;M(#eta)(GeV)", 100,0.25,0.85 ); 
			//hists_mpi0[i] = new TH1F(("mpi0"+to_string(i)+"_"+std::to_string(branchIdx)).c_str(), "Cuts=GeneralCuts;M(#pi_{0}) GeV", 100,0.05,0.25 ); 
			hists_mpi0eta_t[i] = new TH2F(("mpi0eta_t"+to_string(i)+"_"+std::to_string(branchIdx)).c_str(), "Cuts=mMandelstamT_eta;M(#pi_{0}#eta) (GeV);t_{#eta} (GeV^2)", 260, 0.6, 3.2, 80,0,8);
		}
		int histIdx;
		cout << "Initialized histograms" << endl;

		Long64_t nentries = dataTree->GetEntries();

		for (int ientry=0; ientry<nentries; ++ientry){
			dataTree->GetEntry(ientry);
			if ( (selectionString == "tGT1" && new_ptGT1) || (selectionString == "tLT05" && new_ptLT05) || (selectionString == "tGT05LT1" && new_ptGT05LT1) || 
					selectionString == "tAll" ) {  
				if (branchIdx==0){ isUnique = new_isUnique34B_1234B; }
				else { isUnique = new_isUnique12B_1234B; } 
				if ( new_t < tMin || new_m < mMin || new_t>tMax || new_m>mMax ) {continue; }
				idx_t = (int)( (new_t-tMin)/tStep ); 
				idx_m = (int)( (new_m-mMin)/mStep );
				//cout << idx_t << ", " << idx_m << endl;
				histIdx=num_tBins*idx_m+idx_t;
				hists_meta_mpi0[histIdx]->Fill(new_mpi0,new_meta,new_accWeight);
				//hists_meta[histIdx]->Fill(new_meta,new_accWeight);
				//hists_mpi0[histIdx]->Fill(new_mpi0,new_accWeight);
				hists_mpi0eta_t[histIdx]->Fill(new_m,new_t,new_accWeight);
				full_meta_mpi0->Fill(new_mpi0,new_meta,new_accWeight);
				full_meta->Fill(new_meta,new_accWeight);
				full_mpi0->Fill(new_mpi0,new_accWeight);
			}
		}
		cout << "Loaded all the data" << endl;

		logFile << "THESE ARE PROBABLY THE MOST IMPORTANT FITS TO CHECK SINCE THEY ARE THE INITIALIZATIONS OF THE BINNED FITS" << endl;
		logFile << " -------------------------------------------------------------------------------------------------------" << endl;
		logFile << "Distribution, fitStatus, chiSqPerDOF" << endl;

		///////////////////////////////////////////
		// DOING 1D FIT TO USE AS INITIALIZATION FOR 2D
		///////////////////////////////////////////
		cout << "-----------------------\nSTARTING FULL ETA FIT!\n---------------------------- " << endl;
		double fitPars_eta[5];
		int fit_nentries=full_meta_mpi0->GetEntries();
		cout << "TOTAL ENTRIES: " << fit_nentries << endl;
		//gStyle->SetOptFit();
		TF1 * feta = new TF1("meta_fit",gaus,yfitMin,yfitMax,numDOFsig1D); 
		feta->SetParameters(2000,0.55,0.03,0, 0);
		feta->SetParLimits(0,0,5000);
		feta->SetParLimits(1,0.52,0.585);
		feta->SetParLimits(2,0,0.1);
		feta->SetParLimits(4,0,1000); 
		Int_t fitStatus_eta = full_meta->Fit(feta,"RLV");
		full_meta->SetTitle(("ChiSqPerDOF: "+to_string(feta->GetChisquare()/feta->GetNDF())).c_str());
		full_meta->Draw();
		feta->GetParameters(fitPars_eta);
		allCanvases->SaveAs(("deckPlots/"+selectionString+"/"+branchName+"/project_full_meta.png").c_str());

		logFile << "M(eta), " << fitStatus_eta << ", " << feta->GetChisquare()/feta->GetNDF() << endl;

		cout << "-----------------------\nSTARTING FULL PI0 FIT!\n---------------------------- " << endl;
		double fitPars_pi0[5];
		allCanvases->Clear();
		TF1 * fpi0 = new TF1("mpi0_fit",gaus,xfitMin,xfitMax,numDOFsig1D); 
		fpi0->SetParameters(500,0.135,0.015,0, 0);
		fpi0->SetParLimits(0,0,1000);
		fpi0->SetParLimits(1,0.125,0.145);
		fpi0->SetParLimits(2,0.005,0.02);
		fpi0->SetParLimits(4,0,1000); 
		Int_t fitStatus_pi0  = full_mpi0->Fit(fpi0,"RLV");
		full_mpi0->SetTitle(("ChiSqPerDOF: "+to_string(fpi0->GetChisquare()/fpi0->GetNDF())).c_str());
		full_mpi0->Draw();
		fpi0->GetParameters(fitPars_pi0);
		allCanvases->SaveAs(("deckPlots/"+selectionString+"/"+branchName+"/project_full_mpi0.png").c_str());
		logFile << "M(pi0), " << fitStatus_pi0 << ", " << fpi0->GetChisquare()/fpi0->GetNDF() << endl;

		cout << "-----------------------\nSTARTING FULL PI0ETA FIT!\n---------------------------- " << endl;
		TF2 * f2 = new TF2("f2",g2,xfitMin,xfitMax,yfitMin,yfitMax,numDOFsig); 
		TF2 * fflat2 = new TF2("fflat2",gflat2,xfitMin,xfitMax,yfitMin,yfitMax,numDOFsigFlat); 
		double peakBin=1200;
		// We will try to estimate the initialization of the 2D gaussian using the 1D gaussians
		// we get the peak from looking at the histogram. Take the peak*2*pi*sigmaX*sigmaY
		// take directly the mean, sigma from the pi0 and eta 1D fits since they dont scale
		// for the flat bkg we can take some avg between the two const parameters, and for the linear we can take the values from the individual fits. These 3 have to be divided by 100 since
		// 	it might make sense the scale grows with the number of bins we project onto it. In this case the number of bins is 100.
		double full_initParams[8]={peakBin*2*TMath::Pi()*fitPars_pi0[2]*fitPars_eta[2],fitPars_pi0[1], fitPars_pi0[2], fitPars_eta[1], fitPars_eta[2], (fitPars_eta[3]+fitPars_pi0[3])/2/100, fitPars_eta[4]/100, fitPars_pi0[4]/100};
		double full_convParams[8];
		double fitPars[8];
		double fitParErrors[8];
		const double *pointParError;
		f2->SetParameters(full_initParams[0],full_initParams[1],full_initParams[2],full_initParams[3],full_initParams[4],full_initParams[5],full_initParams[6],full_initParams[7]); 
		f2->SetParLimits(0,0,full_initParams[0]*10);
		f2->SetParLimits(1,0.13,0.14);
		f2->SetParLimits(2,0.005,0.015);
		f2->SetParLimits(3,0.52,0.585);
		f2->SetParLimits(4,0,0.1);
		allCanvases->Clear();
		f2->SetLineColorAlpha(kRed, 1);
		f2->SetContour(8);
		f2->SetLineWidth(2);
		Int_t fitStatus_pi0eta = full_meta_mpi0->Fit(f2,"RLV");
		full_meta_mpi0->SetTitle(("ChiSqPerDOF: "+to_string(f2->GetChisquare()/f2->GetNDF())).c_str());
		f2->GetParameters(full_convParams);
		full_meta_mpi0->Draw("COLZ");
		//full_meta_mpi0->Draw("SURF2");
		//f2->Draw("SURF SAME");
		allCanvases->SaveAs(("deckPlots/"+selectionString+"/"+branchName+"/full_meta_mpi0.png").c_str());	
		allCanvases->Clear();
		logFile << "M(pi0eta), " << fitStatus_pi0eta << ", " << f2->GetChisquare()/f2->GetNDF() << endl;


		cout << "-----------------------\nSTARTING BINNED FITS!\n---------------------------- " << endl;
		logFile << " -------------------------------------------------------------------------------------------------------" << endl;
		logFile << "id, fitStatus, chiSqPerDOF, numEntries"  << endl;


		///////////////////////////////////////////
		// STARTING THE INDIVIDUAL FITS IN EACH BIN
		///////////////////////////////////////////

		TF1 * feta_projected = new TF1("meta_projected",g1_eta,yfitMin,yfitMax,numDOF_1D); 
		TF1 * fpi0_projected = new TF1("mpi0_projected",g1_pi0,xfitMin,xfitMax,numDOF_1D); 
		double scaled_initParams[8];
		double scaleFactor;
		TRandom randomGenerator;
		int maxResamplings=10;
		double scales[3]={0.1,0.2,0.3};
		int scaleThresholds[3]={2,4,6};
		for(int i=0; i<numHists; ++i){
        		scaleFactor = hists_meta_mpi0[i]->GetEntries()/fit_nentries;
			for (int j=0; j<8; ++j){
				if ( j==0 || j==5 || j==6 || j==7){
					// the only paramters that scale are the massHist of  the 2D gaus and the paramters for the planar bkg
					scaled_initParams[j] = full_convParams[j]*scaleFactor;
				}
				else {
					// We cant scale everything since the mean and variance shouldn't change
					scaled_initParams[j] = full_convParams[j];
				}
			}
			f2->SetParameters(scaled_initParams[0],scaled_initParams[1],scaled_initParams[2],scaled_initParams[3],scaled_initParams[4],scaled_initParams[5],scaled_initParams[6],scaled_initParams[7]); 
			f2->SetParLimits(0,0,hists_meta_mpi0[i]->GetEntries()); // just restrict to be positive
			f2->SetParLimits(1,0.13,0.14);
			f2->SetParLimits(2,0.005,0.015);
			f2->SetParLimits(3,0.52,0.585);
			f2->SetParLimits(4,0,0.1);
			//f2->SetParLimits(5,0,scaled_initParams[5]*10);
			//f2->SetParLimits(6,-scaled_initParams[5]*10,scaled_initParams[5]*10);
			//f2->SetParLimits(7,-scaled_initParams[6]*10,scaled_initParams[6]*10);
			cout << "\n\n*******************\nFITTING NEXT BIN: " << i << " ************************"<<endl;
			TCanvas *allCanvases_binChecks = new TCanvas("anyHists_binChecks","",1440,900);
			allCanvases_binChecks->Divide(3,1);
			allCanvases_binChecks->cd();
			TPad *pad1 = new TPad("pad1", "",0.3,0.5,0.7,1.0);
			TPad *pad2 = new TPad("pad2", "",0.0,0.0,0.5,0.5);
			TPad *pad3 = new TPad("pad3", "",0.5,0.0,1.0,0.5);
			TPad *pad4 = new TPad("pad4", "",0.1,0.65,0.25,0.85);
			pad1->Draw();
			pad2->Draw();
			pad3->Draw();
			pad4->Draw();

			pad1->cd();
			gStyle->SetOptStat(0);

			Int_t fitStatus = hists_meta_mpi0[i]->Fit(f2,"RLE");

			if ( fitStatus != 0 ) {  
				f2->SetParameters(scaled_initParams[0],scaled_initParams[1],scaled_initParams[2],scaled_initParams[3],scaled_initParams[4],scaled_initParams[5],0,0); 
				f2->FixParameter(6,0);
				f2->FixParameter(7,0);
				fitStatus = hists_meta_mpi0[i]->Fit(f2,"RLE");
				if (fitStatus != 0) { cout << "FUDGE WE STILL NEED TO FIX THIS.... EXITING..." << endl; }
			}
			
			dof=f2->GetNDF();
			chiSq=f2->GetChisquare();
			chiSqPerDOF->Fill(chiSq/dof);
			hists_meta_mpi0[i]->SetTitle(("chiSqPerDOF: "+to_string(chiSq/dof)).c_str());
			hists_meta_mpi0[i]->Draw("COLZ");
			logFile << i << ", " << fitStatus << ", " << chiSq/dof << ", " << hists_meta_mpi0[i]->GetEntries() << endl; 
			//f2->SetParameter(0,randomGenerator.Uniform( scaled_initParams[0]*(1-scale),scaled_initParams[0]*(1+scale) ) );
			//f2->SetParameter(5,randomGenerator.Uniform( scaled_initParams[5]*(1-scale),scaled_initParams[5]*(1+scale) ) );
			
			// just resample the amplitudes and copy over the peak,widths
			//double scale;
			//for ( int maxIters=0; maxIters<maxResamplings; ++maxIters ){
			//	if (fitStatus == 0 ) { break; }
			//	if (maxIters==scaleThresholds[0]) { scale = scales[0]; } 
			//	if (maxIters==scaleThresholds[1]) { scale = scales[1]; } 
			//	if (maxIters==scaleThresholds[2]) { scale = scales[2]; } 
			//	cout << "                  *NEXT RESAMPLE*"<<endl;
			//	cout << randomGenerator.Uniform( scaled_initParams[0]*(1-scale),scaled_initParams[0]*(1+scale)) << endl;
			//	f2->SetParameter(0,randomGenerator.Uniform( scaled_initParams[0]*(1-scale),scaled_initParams[0]*(1+scale) ) );
			//	f2->SetParameter(5,randomGenerator.Uniform( scaled_initParams[5]*(1-scale),scaled_initParams[5]*(1+scale) ) );
			//	f2->SetParameter(6,randomGenerator.Uniform( scaled_initParams[6]*(1-scale),scaled_initParams[6]*(1+scale) ) );
			//	f2->SetParameter(7,randomGenerator.Uniform( scaled_initParams[7]*(1-scale),scaled_initParams[7]*(1+scale) ) );
			//	fitStatus = hists_meta_mpi0[i]->Fit(f2,"RLEV");
			//}
			//if (fitStatus != 0) { cout << "FUDGE WE STILL NEED TO FIX THIS.... EXITING..." << endl; exit(0); }
			
			f2->GetParameters(fitPars);

			feta_projected->SetParameters(fitPars[0],fitPars[1],fitPars[2],fitPars[3],fitPars[4],fitPars[5],fitPars[6],fitPars[7],1/0.002,xfitMin,xfitMax);
			fpi0_projected->SetParameters(fitPars[0],fitPars[1],fitPars[2],fitPars[3],fitPars[4],fitPars[5],fitPars[6],fitPars[7],1/0.006,yfitMin,yfitMax);// yfitMax-yfitMin);
			TH1D * hists_meta = hists_meta_mpi0[i]->ProjectionY("proj_eta", 0, -1);
			TH1D * hists_mpi0 = hists_meta_mpi0[i]->ProjectionX("proj_pi0", 0, -1);
			hists_meta->SetTitle("");
			hists_mpi0->SetTitle("");

			pointParError = f2->GetParErrors();

			pad2->cd();
			gStyle->SetOptStat(0);
			hists_meta->Draw("E1");
			feta_projected->SetLineColor(kRed);
			feta_projected->Draw("SAME");

			pad3->cd();
			gStyle->SetOptStat(0);
			hists_mpi0->Draw("E1");
			fpi0_projected->SetLineColor(kRed);
			fpi0_projected->Draw("SAME");

			pad4->cd();
			TPaveText *histID_pt = new TPaveText(0,0,1,1);
			histID_pt->AddText(("hist id: "+std::to_string(i)).c_str());
			histID_pt->Draw();

			if ( i == 0 ) {
				//allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/fitDiagnostic.pdf(").c_str(),"pdf");	
				allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/fitDiagnostic.pdf(").c_str(),"pdf");
				cout << "Starting to fill the pdf" << endl;
			}
			else if (i == numHists-1) { 
				//allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/fitDiagnostic.pdf)").c_str(),"pdf");	
				allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/fitDiagnostic.pdf)").c_str(),"pdf");
				cout << "Finished filling the pdf" << endl;
			}
			else {
				//allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/fitDiagnostic.pdf").c_str(),"pdf");	
				allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/fitDiagnostic.pdf").c_str(),"pdf");
				cout << "Adding to the pdf" << endl;
			}

			//allCanvases->Clear();
			//hists_mpi0eta_t[i]->Draw("COLZ");
			//allCanvases->SaveAs(("deckPlots/mpi0eta_t-"+to_string(i)+".png").c_str());	


			double unscaledYield = differentialArea*fitPars[0];
			double unscaledYieldError = differentialArea*pointParError[0];
			double yieldErrorTerm=(unscaledYieldError/unscaledYield)*(unscaledYieldError/unscaledYield);
			double efficErrorTerm=(efficiencies_error[branchIdx][i]/efficiencies[branchIdx][i])*(efficiencies_error[branchIdx][i]/efficiencies[branchIdx][i]);
			yields[branchIdx][i]=unscaledYield/efficiencies[branchIdx][i];
			yieldErrors[branchIdx][i]=abs(yields[branchIdx][i])*sqrt(yieldErrorTerm+efficErrorTerm);
			unscaledYields[branchIdx][i]=unscaledYield;
			unscaledYieldErrors[branchIdx][i]=unscaledYieldError;

			// put this here so we can still properly save the diagnostic plots. Even though we will waste a fit calculation....
			if ( hists_meta_mpi0[i]->GetEntries() < 100 ) { 
				yields[branchIdx][i]=0;
				yieldErrors[branchIdx][i]=0;
				unscaledYields[branchIdx][i]=0;
				unscaledYieldErrors[branchIdx][i]=0;
				logFile << i << ", " << "Not enough events, only: " << hists_meta_mpi0[i]->GetEntries() << endl; 
				continue;
			}

			cout << "Yield, error = " << fitPars[0] << ", " << pointParError[0];
			cout << "   Scaled Yield, Scaled error = " << unscaledYield << ", " << unscaledYieldError;
			cout << "   EffCorrected Yield, Scaled error = " << yields[branchIdx][i] << ", " << yieldErrors[branchIdx][i];
			if (maxYield<yields[branchIdx][i]){
				maxYield=yields[branchIdx][i];
			}
			if (minYield>yields[branchIdx][i]){
				minYield=yields[branchIdx][i];
			}
			if (maxUnscaledYield<unscaledYields[branchIdx][i]){
				maxUnscaledYield=unscaledYields[branchIdx][i];
			}


		}
		cout << "    maxYield = " << maxYield << endl;
		cout << "    minYield = " << minYield << endl;
		if ( minYield < 0 ) {
			cout << "NEGATIVE MIN YIELD NOT GOOD. FIX ME" << endl;
		}
		cout << "    maxUnscaledYield = " << maxUnscaledYield << endl;

		allCanvases->Clear();
		chiSqPerDOF->Draw();
		allCanvases->SaveAs(("deckPlots/"+selectionString+"/"+branchName+"/chiSq.png").c_str());
	}
	
	// **************************************************************************************
	// ********* NOW WE CAN START BUILDING OUT teta/tpi0 MASS HISTOGRAMS *********************
	// **************************************************************************************
	branchIdx=-1;
	int skipLastN=0;
	for ( string branchName: branchNames ) {
		++branchIdx;
		std::vector<string> massBinTitles; 
		for (Int_t i=0; i<num_massBins; ++i){
			massBinTitles.push_back( std::to_string( mMin+i*mStep )+" < M(#pi_{0}#eta) < " + std::to_string( mMin+(i+1)*mStep ) );
		}
		TH1F* massHist;
		TH1F* unscaledMassHist;
		TF1* expFit_t;
		gStyle->SetErrorX(0.000001); // remove the x-error bars


		if ( selectionString == "tLT05" ){
			skipLastN = 0;
		}	


		for (Int_t iMass=0; iMass<massBinTitles.size()-skipLastN; ++iMass){
			cout << "Plotting histogram for mass bin: " << iMass << endl;
			massHist = new TH1F(("massHist_"+std::to_string(iMass)+"_"+std::to_string(branchIdx)).c_str(),(massBinTitles[iMass]+";t_{#eta}/t_{#pi_{0}} GeV^{2};").c_str(),num_tBins,tMin,tMax);	
			unscaledMassHist = new TH1F(("unscaledMassHist_"+std::to_string(iMass)+"_"+std::to_string(branchIdx)).c_str(),(massBinTitles[iMass]+";t_{#eta}/t_{#pi_{0}} GeV^{2};").c_str(),num_tBins,tMin,tMax);	

			if (minYield<1) { 
				massHist->SetAxisRange(1,maxYield*1.1,"Y");
			}
			else { 
				massHist->SetAxisRange(minYield*0.9,maxYield*1.1,"Y");
			}
			gStyle->SetTitleSize(0.1,"t");
			massHist->GetXaxis()->SetLabelSize(0.06);
			massHist->GetYaxis()->SetLabelSize(0.06);
			unscaledMassHist->SetAxisRange(1,maxUnscaledYield*1.1,"Y");
			unscaledMassHist->SetTitleSize(1,"t");
			unscaledMassHist->GetXaxis()->SetLabelSize(0.06);
			unscaledMassHist->GetYaxis()->SetLabelSize(0.06);
			expFit_t = new TF1("expFit_t","expo",tMin+tStep,tMax-2*tStep);//+pol0(2)) ,tMin+tStep,tMax);
			expFit_t->SetLineStyle(2);
			for (Int_t i=0; i<num_tBins; ++i){
				massHist->SetBinContent( i+1, yields[branchIdx][i+num_tBins*iMass]  );
				massHist->SetBinError( i+1, yieldErrors[branchIdx][i+num_tBins*iMass]  );
				unscaledMassHist->SetBinContent( i+1, unscaledYields[branchIdx][i+num_tBins*iMass]  );
				unscaledMassHist->SetBinError( i+1, unscaledYieldErrors[branchIdx][i+num_tBins*iMass]  );
				//cout << "yield in bin: " << i+num_tBins*iMass << " is " << yields[i+num_tBins*iMass] << endl;
			}
			massHist->Fit("expFit_t","RBN");

			cout << "Plotting on pad : " << iMass+1 << endl;
			massHist->SetMarkerStyle(kFullCircle);
			massHist->SetMarkerSize(0.75);
			unscaledMassHist->SetMarkerStyle(kFullCircle);
			unscaledMassHist->SetMarkerSize(0.75);
			double scaleAxis;
			double scaleUnscaledAxis;
			double scaleAxisT;
			if (branchIdx==0){
				allCanvases_unscaledYields->cd( iMass+1 );	
				gStyle->SetOptStat(0);
				unscaledMassHist->SetMarkerColor(kBlue);
				unscaledMassHist->Draw("E1 PMC");
				unscaledMassHist->GetXaxis()->SetTitleSize(0.08);
				unscaledMassHist->GetYaxis()->SetTitleSize(0.08);
				allCanvases_unscaledYields->Update();

				allCanvases_tSlope->cd( iMass+1 );
				gPad->SetLogy();
				gStyle->SetOptStat(0);
				massHist->SetMarkerColor(kBlue);
				massHist->Draw("E1 PMC");
				massHist->GetXaxis()->SetTitleSize(0.08);
				massHist->GetYaxis()->SetTitleSize(0.08);
				allCanvases_tSlope->Update();
				expFit_t->SetLineColor(kBlue);
				expFit_t->Draw("SAME");

				allCanvases_yields->cd( iMass+1 );	
				gStyle->SetOptStat(0);
				massHist->SetMarkerColor(kBlue);
				massHist->Draw("E1 PMC");
				massHist->GetXaxis()->SetTitleSize(0.08);
				massHist->GetYaxis()->SetTitleSize(0.08);
				allCanvases_yields->Update();
				scaleAxis = gPad->GetUymax()/(maxEfficiency*1.1); // this should be the same for all pads since I set the axisRange above
				hist_efficiencies_eta[iMass]->Scale(scaleAxis);  
				hist_efficiencies_eta[iMass]->SetLineColor(kBlue);
				hist_efficiencies_eta[iMass]->Draw("E SAME");
				hist_efficiencies_eta[iMass]->Draw("SAME Lhist");

				massHists_eta.push_back(*massHist);
				expPlusLinFits_eta.push_back(*expFit_t);
				tSlopes_eta.push_back( abs(expFit_t->GetParameter(1)) );
				tSlopesError_eta.push_back( abs(expFit_t->GetParError(1)) );
			}
			else { 
				allCanvases_unscaledYields->cd( iMass+1 );	
				unscaledMassHist->SetMarkerColor(kRed);
				unscaledMassHist->Draw("E1 PMC SAME");
				allCanvases_unscaledYields->Update();

				allCanvases_tSlope->cd( iMass+1 );
				massHist->SetMarkerColor(kRed);
				massHist->Draw("E1 PMC SAME");
				allCanvases_tSlope->Update();
				expFit_t->SetLineColor(kRed);
				expFit_t->Draw("SAME");

				allCanvases_yields->cd( iMass+1 );	
				massHist->SetMarkerColor(kRed);
				massHist->Draw("E1 PMC SAME");
				hist_efficiencies_pi0[iMass]->Scale(scaleAxis);  
				hist_efficiencies_pi0[iMass]->SetLineColor(kRed);
				hist_efficiencies_pi0[iMass]->Draw("E SAME");
				hist_efficiencies_pi0[iMass]->Draw("SAME Lhist");

				massHists_pi0.push_back(*massHist);
				expPlusLinFits_pi0.push_back(*expFit_t);
				tSlopes_pi0.push_back( abs(expFit_t->GetParameter(1)) );
				tSlopesError_pi0.push_back( abs(expFit_t->GetParError(1)) );
			}

			TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
					            gPad->GetUxmax(), gPad->GetUymax(),0,maxEfficiency*1.1,510,"+L");
			axis->SetLineColor(kRed);
			axis->SetLabelColor(kRed);
			axis->SetLabelSize(0.06);
			axis->Draw();
		}
		//cout << "Initial parameter values: " << endl;
		//for ( auto parValue: full_initParams ) {
		//	cout << parValue << endl;
		//}

		//cout << "ALL THE YIELDS: " <<endl;
		//for (auto yield : yields ) {
		//	cout << yield << endl;
		//}
	logFile << "\n\n\n======================= END OF THIS ITERATION ===============================\n" << endl;
    	}

	 TLegend *leg;
	 for (int iMass=0; iMass<num_massBins-skipLastN; ++iMass){
	 	allCanvases_tSlope->cd(iMass+1);
	 	leg = new TLegend(0.1,0.11,0.35,0.35);
	 	leg->AddEntry(&expPlusLinFits_eta[iMass], ("#eta tSlope = "+to_string(tSlopes_eta[iMass])).c_str(), "l");
	 	leg->AddEntry(&expPlusLinFits_pi0[iMass], ("#pi^{0} tSlope = "+to_string(tSlopes_pi0[iMass])).c_str(), "l");
	 	leg->SetBorderSize(0);
	 	leg->SetFillStyle(0);
	 	leg->SetTextSize(0.07);
	 	leg->Draw();
	 }

	 allCanvases_yields->SaveAs(("deckPlots/"+selectionString+"/yields.png").c_str());
	 allCanvases_unscaledYields->SaveAs(("deckPlots/"+selectionString+"/unscaledYields.png").c_str());
	 allCanvases_tSlope->SaveAs(("deckPlots/"+selectionString+"/yields_tSlope.png").c_str());
	 allCanvases_yields->SaveAs(("deckPlots/"+selectionString+"/yields.C").c_str());
	 allCanvases_unscaledYields->SaveAs(("deckPlots/"+selectionString+"/unscaledYields.C").c_str());
	 allCanvases_tSlope->SaveAs(("deckPlots/"+selectionString+"/yields_tSlope.C").c_str());


	 allCanvases->Clear();
	 allCanvases->cd();
	 TH1F* hist_tSlopes_eta = new TH1F("hist_tSlope_eta","; M(#pi^{0}#eta) (GeV); t-slope",num_massBins,mMin,mMax);
	 TH1F* hist_tSlopes_pi0 = new TH1F("hist_tSlope_pi0","; M(#pi^{0}#eta) (GeV); t-slope",num_massBins,mMin,mMax);
	 for ( int iMass=0; iMass < num_massBins-skipLastN; ++iMass ){
	 	hist_tSlopes_eta->SetBinContent( iMass+1, tSlopes_eta[iMass]);
	 	hist_tSlopes_pi0->SetBinContent( iMass+1, tSlopes_pi0[iMass]);
	 	hist_tSlopes_eta->SetBinError( iMass+1, tSlopesError_eta[iMass]);
	 	hist_tSlopes_pi0->SetBinError( iMass+1, tSlopesError_pi0[iMass]);
	 }

	 hist_tSlopes_eta->SetMarkerStyle(kFullCircle);
	 hist_tSlopes_eta->SetMarkerSize(1.5);
	 hist_tSlopes_eta->SetMarkerColor(kBlue);
	 hist_tSlopes_eta->GetXaxis()->SetLabelSize(0.05);
	 hist_tSlopes_eta->GetYaxis()->SetLabelSize(0.05);
	 hist_tSlopes_eta->Draw("E1 PMC");

	 hist_tSlopes_pi0->SetMarkerStyle(kFullCircle);
	 hist_tSlopes_pi0->SetMarkerSize(1.5);
	 hist_tSlopes_pi0->SetMarkerColor(kRed);
	 hist_tSlopes_pi0->Draw("E1 PMC SAME");
	 allCanvases->SaveAs(("deckPlots/"+selectionString+"/tSlopesVsMpi0eta.png").c_str());
	 allCanvases->SaveAs(("deckPlots/"+selectionString+"/tSlopesVsMpi0eta.C").c_str());
	 

	 //// Output the tslope for the eta
	 //allCanvases->Clear(); allCanvases->cd();
	 //allCanvases->SetLogy(0);
	 //deckDiagnosticFile->cd();
	 //for (int iMass=0; iMass<num_massBins-skipLastN; ++iMass){
	 //	allCanvases->Clear();
	 //	massHists_eta[iMass].Draw();
	 //	expPlusLinFits_eta[iMass].SetLineColor(kRed);
	 //	expPlusLinFits_eta[iMass].SetLineStyle(1);
	 //	expPlusLinFits_eta[iMass].Draw("SAME");
	 //	allCanvases->SaveAs(("deckPlots/"+selectionString+"/tSlope/eta_tSlope_"+to_string(iMass)+".png").c_str());
	 //	massHists_eta[iMass].Write();
	 //}
	 //// Output the tslope for the pi0
	 //for (int iMass=0; iMass<num_massBins; ++iMass){
	 //	allCanvases->Clear();
	 //	massHists_pi0[iMass].Draw();
	 //	expPlusLinFits_pi0[iMass].SetLineColor(kRed);
	 //	expPlusLinFits_pi0[iMass].SetLineStyle(1);
	 //	expPlusLinFits_pi0[iMass].Draw("SAME");
	 //	allCanvases->SaveAs(("deckPlots/"+selectionString+"/tSlope/pi0_tSlope_"+to_string(iMass)+".png").c_str());
	 //	massHists_pi0[iMass].Write();
	 //}

}


void makeDeckPlots(){
	//gSystem->Exec("rm -rf deckPlots");
	//gSystem->Exec("mkdir -p deckPlots/tAll/mandelstam_teta_meas");
	//gSystem->Exec("mkdir deckPlots/tAll/mandelstam_tpi0_meas");
	//gSystem->Exec("mkdir deckPlots/tAll/tSlope");
	//makeDeckPlot("tAll");

	//gSystem->Exec("mkdir -p deckPlots/tGT1/mandelstam_teta_meas");
	//gSystem->Exec("mkdir deckPlots/tGT1/mandelstam_tpi0_meas");
	//gSystem->Exec("mkdir deckPlots/tGT1/tSlope");
	//makeDeckPlot("tGT1");

	gSystem->Exec("mkdir -p deckPlots/tLT05/mandelstam_teta_meas");
	gSystem->Exec("mkdir deckPlots/tLT05/mandelstam_tpi0_meas");
	gSystem->Exec("mkdir deckPlots/tLT05/tSlope");
	makeDeckPlot("tLT05");

	//gSystem->Exec("mkdir -p deckPlots/tGT05LT1/mandelstam_teta_meas");
	//gSystem->Exec("mkdir deckPlots/tGT05LT1/mandelstam_tpi0_meas");
	//gSystem->Exec("mkdir deckPlots/tGT05LT1/tSlope");
	//makeDeckPlot("tGT05LT1");
}
