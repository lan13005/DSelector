int numDOFsig=8;
// 2d gaussian with linear bkg
Double_t gaus2D(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	Double_t r2 = Double_t((x[1]-par[3])/par[4]);
	return par[5]+par[6]*(x[0]-par[1])+par[7]*(x[1]-par[3])+par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}

// 2d gaussian with linear bkg
int numDOFsigFlat=6;
Double_t gaus2D_flat(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	Double_t r2 = Double_t((x[1]-par[3])/par[4]);
	return par[5]+par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}

// 2d double gaussian with linear bkg
Double_t doubleGaus2D(Double_t *x, Double_t *par){
	return par[9]+par[10]*(x[0]-par[1])+par[11]*(x[1]-par[5])+par[0]/2/TMath::Pi()*(1/par[2]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]))
	     + par[3]/(par[4]*par[2])*exp(-0.5*((x[0]-par[1])/(par[4]*par[2]))*((x[0]-par[1])/(par[4]*par[2])))
	     )*(1/par[6]*exp(-0.5*((x[1]-par[5])/par[6])*((x[1]-par[5])/par[6]))
	     + par[7]/(par[8]*par[6])*exp(-0.5*((x[1]-par[5])/(par[8]*par[6]))*((x[1]-par[5])/(par[8]*par[6])))
	     );
}

// 2d double gaussian with flat bkg
Double_t doubleGaus2D_flat(Double_t *x, Double_t *par){
	return par[9]+par[0]*(1/par[2]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]))
	     + par[3]/(par[4]*par[2])*exp(-0.5*((x[0]-par[1])/(par[4]*par[2]))*((x[0]-par[1])/(par[4]*par[2])))
	     )*(1/par[6]*exp(-0.5*((x[1]-par[5])/par[6])*((x[1]-par[5])/par[6]))
	     + par[7]/(par[8]*par[6])*exp(-0.5*((x[1]-par[5])/(par[8]*par[6]))*((x[1]-par[5])/(par[8]*par[6])))
	     );
}


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


// 1d double gaussian
Double_t doubleGaus(Double_t *x, Double_t *par){
	return par[5]+par[6]*(x[0]-par[1])+par[0]/par[2]/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]))
	     + par[3]*par[0]/(par[4]*par[2])/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/(par[4]*par[2]))*((x[0]-par[1])/(par[4]*par[2])));

}

// skew gaussian
Double_t skewGaus(Double_t *x, Double_t *par){
	return par[4]+par[5]*(x[0]-par[1])+par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2));
}

// 1d gaussian to fit the M(eta) and M(pi0) distribution
int numDOFsig1D_eta=5;
int numDOFsig1D_pi0=5;
Double_t gaus(Double_t *x, Double_t *par){
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	return par[3]+par[4]*(x[0]-par[1])+par[0]/sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-0.5*r1*r1 );
}

// asymmetric gaus
Double_t asymGaus(Double_t *x, Double_t *par){
	double sigma;
	if (x[0] < par[1] ) { sigma = par[2]; }
	else { sigma = par[3]; }
	Double_t r1 = Double_t((x[0]-par[1])/sigma);
	return par[4]+par[5]*(x[0]-par[1])+par[0]/sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-0.5*r1*r1 );
}

void makeDeckPlot(string selectionString){
	cout << " --------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << " --------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << "                                             BEGIN " << selectionString << " FIT CALCULATION" << endl;
	cout << " --------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << " --------------------------------------------------------------------------------------------------------------------------" << endl;



	// *********************************************
	// *********** PREPARING THE CODE ***************
	// *********************************************
	//gStyle->SetErrorX(0.000001); // remove the x-error bars
   	TFile *deckDiagnosticFile = new TFile("deckDiagnosticHists.root", "RECREATE");
    	ofstream logFile;
    	logFile.open(("deckPlots/"+selectionString+"/failedFittingPlotIDs.txt").c_str());

	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	TCanvas *allCanvases_tSlope = new TCanvas("anyHists_tSlope","Blue=Eta Red=Pi0 DarkGray=teta Gray=tpi0",1440,900);
	TCanvas *allCanvases_yields = new TCanvas("anyHists_yields","Blue=Eta Red=Pi0 DarkGray=teta Gray=tpi0",1440,900);
	TCanvas *allCanvases_unscaledYields = new TCanvas("anyHists_unscaledYields","Blue=Eta Red=Pi0 DarkGray=teta Gray=tpi0",1440,900);

	allCanvases_yields->SetLeftMargin(0.15);
	allCanvases_unscaledYields->SetLeftMargin(0.15);
	allCanvases_tSlope->SetBottomMargin(0.2);
	allCanvases_yields->SetBottomMargin(0.2);
	allCanvases_unscaledYields->SetBottomMargin(0.2);

	allCanvases_yields->Divide(3,4,0,0);
	allCanvases_unscaledYields->Divide(3,4,0,0);
	allCanvases_tSlope->Divide(3,4,0,0);

	// THese will define the bins
	int num_tBins=14;
	double tMin=0;
	double tMax=2.8;
	const int num_massBins=12;
	double mMin=1.6;
	double mMax=2.8;
	double tStep=(tMax-tMin)/num_tBins;
	double mStep=(mMax-mMin)/num_massBins;
	const int numHists = (const int)num_tBins*num_massBins;
	cout << "numHists: " << numHists << endl;



	// *********************************************
	// *********** CALCULATE YIELDS IN BINS  ***************
	// *********************************************
	//string dataFileLoc[2] = {"pi0eta_data_treeFlat_DSelector.root","pi0eta_flat_21t_treeFlat_DSelector.root"}; 
	//string dataTreeName[2] = {"pi0eta_datatree_flat","pi0eta_flat_21ttree_flat"}; 
	string dataFileLoc[2] = {"degALL_data_2017_treeFlat_DSelector.root","degALL_acc_2017_treeFlat_DSelector.root"}; 
	string dataTreeName[2] = {"degALL_data_2017_tree_flat","degALL_acc_2017_tree_flat"}; 
	string dataTypes[2]={"data","reco"};
	string branchNames[2]={"mandelstam_teta_meas","mandelstam_tpi0_meas"};
	std::vector<std::vector<double>> tSlopes_etas;
	std::vector<std::vector<double>> tSlopes_pi0s;
	std::vector<std::vector<double>> tSlopesError_etas;
	std::vector<std::vector<double>> tSlopesError_pi0s;
	double yields[2][numHists];
	double yieldErrors[2][numHists];
	double unscaledYields[2][2][numHists]; // unscaledYields will hold data and reco fitted yields. yields on the other hand will just hold the efficiency correct data.
	double unscaledYieldErrors[2][2][numHists];
	double entriesInMetaMpi0[2][2][numHists]; // will give us an idea if the fit might have been good. If there were like 10 events probably a fit would probably be noisy
	// Create some variabls to load the data from the tree
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
	// Some variables to setup histograms with certain bins and widths
	double yBins[2] = {0.25,0.85};
	double xBins[2] = {0.05,0.25};
	int nxBins=100;
	int nyBins=100;
	double xWidth = (xBins[1]-xBins[0])/nxBins;
	double yWidth = (yBins[1]-yBins[0])/nyBins;
	double inverseDifferentialArea = 1/(xWidth*yWidth);
	cout << "Inverse Differential Area:"<< inverseDifferentialArea << endl;
	// First we loop over the dataTypes: data or acc. 
	for ( int dType=0; dType < sizeof(dataTypes)/sizeof(string); ++dType ){	
		logFile << "\n\n" << dataTypes[dType] << "\n-----------------------------------\n-------------------------------" << endl;
		TFile* dataFile = TFile::Open(dataFileLoc[dType].c_str());
		cout << "Loaded root file: " << dataFileLoc[dType] << endl;
		TTree *dataTree;
		dataFile->GetObject(dataTreeName[dType].c_str(),dataTree);
		cout << "Loaded tree: " << dataTreeName[dType] << endl;

		// Match the addresses from the flat tree
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

		int branchIdx=-1;
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
			TH2F* hists_mpi0eta_t[numHists];

			double xfitMin=0.085;
			double xfitMax=0.185;

			double yfitMin;
			double yfitMax;

			if (selectionString=="tGT1" || selectionString=="tAll" || selectionString=="tGT05LT1"){
				//  These were fitRange we used for the analysis presentation
				yfitMin=0.45;
				yfitMax=0.65;
			}
			else if (selectionString=="tLT05") { 
				//  Since tLT05 doesnt fit easily we can shorten the range in the eta
				yfitMin=0.46;
				yfitMax=0.625;
			}
			cout << "Current t selection will use [" << yfitMin << "," << yfitMax << "] as the fit range" << endl;

			double chiSq;
			double dof;

			TH1F* chiSqPerDOF = new TH1F(("chiSqPerDOF_"+std::to_string(branchIdx)).c_str(), "#chi^{2}/n_{dof}; #chi^{2}/n_{dof}; Events / 0.1", 30,0,3 );
			TH1F* full_meta = new TH1F(("full_meta_"+std::to_string(branchIdx)).c_str(), ";M(#eta)(GeV)", nxBins,yBins[0],yBins[1] );
			TH1F* full_mpi0 = new TH1F(("full_mpi0_"+std::to_string(branchIdx)).c_str(), ";M(#pi^0)(GeV)",nyBins,xBins[0],xBins[1] );
			TH2F* full_meta_mpi0 = new TH2F(("full_meta_mpi0_"+std::to_string(branchIdx)).c_str(), "Cuts=GeneralCuts;M(#pi^{0}) GeV;M(#eta)(GeV)", 100,0.05,0.25,100,0.25,0.85 );
			
			for (int i=0; i< numHists; ++i){
				hists_meta_mpi0[i] = new TH2F(("meta_mpi0"+to_string(i)+"_"+std::to_string(branchIdx)).c_str(), "Cuts=GeneralCuts;M(#pi^{0}) GeV;M(#eta)(GeV)", nxBins,xBins[0],xBins[1],nyBins,yBins[0],yBins[1] ); 
				hists_mpi0eta_t[i] = new TH2F(("mpi0eta_t"+to_string(i)+"_"+std::to_string(branchIdx)).c_str(), "Cuts=mMandelstamT_eta;M(#pi^{0}#eta) (GeV);t_{#eta} (GeV^2)", 260, 0.6, 3.2, 80,0,8);
			}
			int histIdx;
			cout << "Initialized histograms" << endl;

			Long64_t nentries = dataTree->GetEntries();

			for (int ientry=0; ientry<nentries; ++ientry){
				dataTree->GetEntry(ientry);
				// here is where we decide whether or not we are going to fill a histogram. If it matches the selection string and satisfies the corresponding bool then fill.
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
			//
			cout << "-----------------------\nSTARTING FULL ETA FIT!\n---------------------------- " << endl;
			double fitPars_eta[numDOFsig1D_eta];
			int fit_nentries=full_meta_mpi0->GetEntries();
			cout << "TOTAL ENTRIES: " << fit_nentries << endl;
			//gStyle->SetOptFit();
			TF1 * feta = new TF1("meta_fit",asymGaus,yfitMin,yfitMax,numDOFsig1D_eta); 
			TF1 * feta_before = new TF1("meta_fit_before",asymGaus,yfitMin,yfitMax,numDOFsig1D_eta); 
			// FIT PARAMS FOR SKEW GAUSSIAN
			// 15299.2, 0.5506, 0.0123107, -0.142242, 513.223, 1421.38, 0
			//feta->SetParameters(19000,0.55,0.01,-0.2,300,1000);
			//feta_before->SetParameters(19000,0.55,0.01,-0.2,300,1000);
			// FIT PARAMS FOR DOUBLE GAUSSIAN
			//feta->SetParameters(800,0.541,0.02,0.4,2.5,300,0);
			// FIT PARAMS FOR GAUSSIAN
			feta->SetParameters(2000, 0.55,0.01, 0.01, 0, 0);
			feta->SetParLimits(0, 0, 5000);
			feta->SetParLimits(1, 0.52, 0.585);
			feta->SetParLimits(2, 0, 0.1);
			feta->SetParLimits(3, 0, 0.1);
			feta->SetParLimits(4, 0, 1000); 
			feta->SetLineColor(kRed);
			Int_t fitStatus_eta = full_meta->Fit(feta,"RLQ");
			double chiPerDOF = feta->GetChisquare()/feta->GetNDF();
			for(int iFit=0; iFit < 25; ++iFit){
				if (chiPerDOF > 200){
					if (iFit==49){ 
						cout << "Meta full fit never converged... exiting" << endl; 
						feta->Draw();
						allCanvases->SaveAs(("deckPlots/"+selectionString+"/"+branchName+"/"+dataTypes[dType]+"/bad-project_full_meta.png").c_str());
						exit(0); 
					}
					double rand_par0 = rand() % 2000;
					double rand_par4 = rand() % 1000;
					feta->SetParameters(2000, 0.55,0.01, 0.01, 0, 0);
					feta->SetParLimits(0, 0, 10000);
					feta->SetParLimits(1, 0.52, 0.585);
					feta->SetParLimits(2, 0, 0.1);
					feta->SetParLimits(4, 0, 5000); 
					feta->SetLineColor(kRed);
					Int_t fitStatus_eta = full_meta->Fit(feta,"RLQ");
				}
			}
			full_meta->SetTitle(("ChiSqPerDOF: "+to_string(chiPerDOF)).c_str());
			full_meta->Draw();
			feta_before->SetLineColor(kBlue);
			//feta_before->Draw("SAME");
			feta->GetParameters(fitPars_eta);
			cout << "full_eta fitted parameters: " << endl;
			for ( auto full_initParam : fitPars_eta ) {
				cout << full_initParam << ", ";
			}
			cout << endl;
			allCanvases->SaveAs(("deckPlots/"+selectionString+"/"+branchName+"/"+dataTypes[dType]+"/project_full_meta.png").c_str());

			logFile << "M(eta), " << fitStatus_eta << ", " << feta->GetChisquare()/feta->GetNDF() << endl;

			cout << "-----------------------\nSTARTING FULL PI0 FIT!\n---------------------------- " << endl;
			double fitPars_pi0[numDOFsig1D_pi0];
			allCanvases->Clear();
			TF1 * fpi0 = new TF1("mpi0_fit",gaus,xfitMin,xfitMax,numDOFsig1D_pi0); 
			// FIT PARAMS FOR DOUBLE GAUSSIAN
			//fpi0->SetParameters(100,0.134,0.005,2.5,1.5,300,0);
			// FIT PARAMS FOR GAUSSIAN
			fpi0->SetParameters(500,0.135,0.015,0, 0);
			fpi0->SetParLimits(0,0,1000);
			fpi0->SetParLimits(1,0.125,0.145);
			fpi0->SetParLimits(2,0.005,0.02);
			fpi0->SetParLimits(4,0,1000); 
			Int_t fitStatus_pi0  = full_mpi0->Fit(fpi0,"RLQ");
			full_mpi0->SetTitle(("ChiSqPerDOF: "+to_string(fpi0->GetChisquare()/fpi0->GetNDF())).c_str());
			full_mpi0->Draw();
			fpi0->GetParameters(fitPars_pi0);
			cout << "full_pi0 fitted parameters: " << endl;
			for ( auto full_initParam : fitPars_pi0 ) {
				cout << full_initParam << ", ";
			}
			cout << endl;
			allCanvases->SaveAs(("deckPlots/"+selectionString+"/"+branchName+"/"+dataTypes[dType]+"/project_full_mpi0.png").c_str());
			logFile << "M(pi0), " << fitStatus_pi0 << ", " << fpi0->GetChisquare()/fpi0->GetNDF() << endl;

			cout << "-----------------------\nSTARTING FULL PI0ETA FIT!\n---------------------------- " << endl;
			// 2D DOUBLE GAUS
			//TF2 * f2 = new TF2("f2",doubleGaus2D,xfitMin,xfitMax,yfitMin,yfitMax,numDOFsig); 
			//TF2 * fflat2 = new TF2("fflat2",doubleGaus2D_flat,xfitMin,xfitMax,yfitMin,yfitMax,numDOFsigFlat); 
			// 2D GAUS
			TF2 * f2 = new TF2("f2",gaus2D,xfitMin,xfitMax,yfitMin,yfitMax,numDOFsig); 
			TF2 * fflat2 = new TF2("fflat2",gaus2D_flat,xfitMin,xfitMax,yfitMin,yfitMax,numDOFsigFlat); 
			// We will try to estimate the initialization of the 2D gaussian using the 1D gaussians
			// we get the peak from looking at the histogram. Take the peak*2*pi*sigmaX*sigmaY
			// take directly the mean, sigma from the pi0 and eta 1D fits since they dont scale
			// for the flat bkg we can take some avg between the two const parameters, and for the linear we can take the values from the individual fits. These 3 have to be divided by 100 since
			// 	it might make sense the scale grows with the number of bins we project onto it. In this case the number of bins is 100.
			// initial parameters for 2D double gaussian
			//double full_initParams[12]={(fitPars_pi0[0]/nxBins+fitPars_eta[0]/nyBins)/2,fitPars_pi0[1], fitPars_pi0[2], fitPars_pi0[3], fitPars_pi0[4], fitPars_eta[1], fitPars_eta[2], fitPars_eta[3], fitPars_eta[4], (fitPars_eta[5]/nyBins+fitPars_pi0[5]/nxBins)/2, fitPars_pi0[6]/nxBins, fitPars_eta[6]/nyBins};

			// initial parameters for 2D gaussian
			double peakBin=1200;
			double full_initParams[8]={peakBin*2*TMath::Pi()*fitPars_pi0[2]*fitPars_eta[2],fitPars_pi0[1], fitPars_pi0[2], fitPars_eta[1], fitPars_eta[2], (fitPars_eta[3]+fitPars_pi0[3])/2/100, fitPars_eta[4]/100, fitPars_pi0[4]/100};


			cout << "full_pi0eta fitted parameters: " << endl;
			for ( auto full_initParam : full_initParams ) {
				cout << full_initParam << ", ";
			}
			cout << endl;

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
			Int_t fitStatus_pi0eta = full_meta_mpi0->Fit(f2,"RLQ");
			full_meta_mpi0->SetTitle(("ChiSqPerDOF: "+to_string(f2->GetChisquare()/f2->GetNDF())).c_str());
			f2->GetParameters(full_convParams);
			full_meta_mpi0->Draw("COLZ");
			//full_meta_mpi0->Draw("SURF2");
			//f2->Draw("SURF SAME");
			allCanvases->SaveAs(("deckPlots/"+selectionString+"/"+branchName+"/"+dataTypes[dType]+"full_meta_mpi0.png").c_str());	
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
			double randomInitParams[8];
			int maxItersToResampleRandomly = 10;
			TRandom* randomScale = new TRandom();
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
				f2->SetParLimits(1,0.125,0.145);
				f2->SetParLimits(2,0.0025,0.0175);
				f2->SetParLimits(3,0.45,0.65);
				f2->SetParLimits(4,0,0.15);
				f2->SetParLimits(5,0,scaled_initParams[5]*10);
				f2->SetParLimits(6,-scaled_initParams[5]*10,scaled_initParams[5]*10);
				f2->SetParLimits(7,-scaled_initParams[6]*10,scaled_initParams[6]*10);
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

				Int_t fitStatus = hists_meta_mpi0[i]->Fit(f2,"RLQE");
				
				int currentIter = 0;
				while ( fitStatus != 0 ) {  
					if( currentIter < maxItersToResampleRandomly ) {
						cout << "RANDOMLY RESAMPLING THE INIT PARAMETERS TO TRY TO GET CONVERGENCE" << endl;
						for (int j=0; j<8; ++j){
							randomInitParams[j] = scaled_initParams[j]*randomScale->Uniform(0.9,1.1);
						}
						f2->SetParameters(randomInitParams[0],randomInitParams[1],randomInitParams[2],randomInitParams[3],randomInitParams[4],randomInitParams[5],0,0); 
						f2->FixParameter(6,0);
						f2->FixParameter(7,0);
						fitStatus = hists_meta_mpi0[i]->Fit(f2,"RLQE");
						if (fitStatus != 0) { cout << "FUDGE WE STILL NEED TO FIX THIS.... EXITING..." << endl; }
						++currentIter;
					}
					else{
						break;
					}
				}
				
				dof=f2->GetNDF();
				chiSq=f2->GetChisquare();
				chiSqPerDOF->Fill(chiSq/dof);
				hists_meta_mpi0[i]->SetTitle(("chiSqPerDOF: "+to_string(chiSq/dof)+" -- FitStatus: "+to_string(fitStatus)).c_str());
				hists_meta_mpi0[i]->Draw("COLZ");
				logFile << i << ", " << fitStatus << ", " << chiSq/dof << ", " << hists_meta_mpi0[i]->GetEntries() << endl; 

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
					allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/"+dataTypes[dType]+"/fitDiagnostic.pdf(").c_str(),"pdf");
					cout << "Starting to fill the pdf" << endl;
				}
				else if (i == numHists-1) { 
					//allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/fitDiagnostic.pdf)").c_str(),"pdf");	
					allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/"+dataTypes[dType]+"/fitDiagnostic.pdf)").c_str(),"pdf");
					cout << "Finished filling the pdf" << endl;
				}
				else {
					//allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/fitDiagnostic.pdf").c_str(),"pdf");	
					allCanvases_binChecks->Print(("deckPlots/"+selectionString+"/"+branchName+"/"+dataTypes[dType]+"/fitDiagnostic.pdf").c_str(),"pdf");
					cout << "Adding to the pdf" << endl;
				}

				unscaledYields[dType][branchIdx][i] = inverseDifferentialArea*fitPars[0];
				unscaledYieldErrors[dType][branchIdx][i] = inverseDifferentialArea*pointParError[0];
				// we will keep track of the number of entries in our binned 2D fit. Will allow us to make a threshold on when the fit should be considered
				entriesInMetaMpi0[dType][branchIdx][i] = hists_meta_mpi0[i]->GetEntries();

				cout << "Unscaled Yield: " << inverseDifferentialArea*fitPars[0] << endl;
				cout << "Unscaled Yield Errors: " << inverseDifferentialArea*pointParError[0] << endl;
				cout << "Entries in fit: " << hists_meta_mpi0[i]->GetEntries() << endl;

			} // closes the hist for loop
			allCanvases->Clear();
			chiSqPerDOF->Draw();
			allCanvases->SaveAs(("deckPlots/"+selectionString+"/"+branchName+"/"+dataTypes[dType]+"/chiSq.png").c_str());
      		} // closes the teta or tpi0 for loop
      } // closes the dataTypes for loop
	

	// *********************************************
	// *********** CALCULATE EFFICIENCIES ***************
	// *********************************************
	//TFile* genFile = TFile::Open("flat_21t_gen_hists_DSelector_pi0eta.root");
	TFile* genFile = TFile::Open("degALL_gen_2017_hists_DSelector.root");
	TH1F *tetaVsMpi0eta_genCounts;
	TH1F *tpi0VsMpi0eta_genCounts;
	genFile->GetObject( ("tetaVsMpi0eta_genCounts_"+selectionString).c_str(),tetaVsMpi0eta_genCounts);
	genFile->GetObject( ("tpi0VsMpi0eta_genCounts_"+selectionString).c_str(),tpi0VsMpi0eta_genCounts);

	 // Save the histograms for easy access
	allCanvases->cd(); allCanvases->Clear(); allCanvases->SetLogy(); 
	tetaVsMpi0eta_genCounts->Draw();
	allCanvases->SaveAs(("deckPlots/"+selectionString+"/tetaVsMpi0eta.png").c_str());
	allCanvases->Clear();allCanvases->SetLogy();
	tpi0VsMpi0eta_genCounts->Draw();
	allCanvases->SaveAs(("deckPlots/"+selectionString+"/tpi0VsMpi0eta.png").c_str());

	allCanvases->SetLogy(0);

	TH1F *hist_efficiencies_pi0[num_massBins];
	TH1F *hist_efficiencies_eta[num_massBins];
	for (int iHist=0; iHist < num_massBins; ++iHist){
		hist_efficiencies_pi0[iHist] = new TH1F("","",num_tBins,tMin,tMax);
		hist_efficiencies_eta[iHist] = new TH1F("","",num_tBins,tMin,tMax);
	}
	double c_teta_genCounts;
	double c_tpi0_genCounts;
	double c_teta_recCounts;
	double c_tpi0_recCounts;
	double c_teta_recCountErrors;
	double c_tpi0_recCountErrors;
	std::vector<double > efficiencies_pi0; efficiencies_pi0.reserve(numHists);
	std::vector<double > efficiencies_eta; efficiencies_eta.reserve(numHists);
	std::vector< std::vector<double> > efficiencies; 
	std::vector<double > efficiencies_error_pi0; efficiencies_error_pi0.reserve(numHists);
	std::vector<double > efficiencies_error_eta; efficiencies_error_eta.reserve(numHists);
	std::vector< std::vector<double> > efficiencies_error; 
	double maxEfficiency=DBL_MIN;
	bool skipCalc;

	// I SET THESE BINS THIS WAY
	// bin0 = underflow
	// bin1 = -1 which we used to hold all the data that wasn't in a bin
	for (int i=2; i<numHists+2; ++i) { 
		// To make the code more understandable we will just calculate the efficiencies first then fill the histogram. There will be too much going on with the array indicies if
		//     we have to consider the overflow bin and the shifting of the histograms we read in since they also have an overflow and an extra bin that holds non-bin events. 
		//     ALSO: we have to find the maximum value first
		int j=i-2; // shifted i
		c_teta_genCounts = tetaVsMpi0eta_genCounts->GetBinContent(i);
		c_tpi0_genCounts = tpi0VsMpi0eta_genCounts->GetBinContent(i);
		c_teta_recCounts = unscaledYields[1][0][j]; // [data/reco][eta/pi0][MassBins]
		c_tpi0_recCounts = unscaledYields[1][1][j]; // [data/reco][eta/pi0][MassBins]
		c_teta_recCountErrors = unscaledYieldErrors[1][0][j]; // [data/reco][eta/pi0][MassBins]
		c_tpi0_recCountErrors = unscaledYieldErrors[1][1][j]; // [data/reco][eta/pi0][MassBins]
		
		skipCalc=false;
		cout << " ** IF COUNTS <=0 THEN WE SET THEM =0 **" << endl;
		if ( c_teta_genCounts <= 0 || c_teta_recCounts <= 0 ) { efficiencies_eta.push_back(1); efficiencies_error_eta.push_back(0); skipCalc=true; cout << "SETTNG TO EFFICIENCY=1 FOR ETA" << endl;} 
		if ( c_tpi0_genCounts <= 0 || c_tpi0_recCounts <= 0 ) { efficiencies_pi0.push_back(1); efficiencies_error_pi0.push_back(0); skipCalc=true; cout << "SETTNG TO EFFICIENCY=1 FOR PI0" << endl;} 

		cout << "c_teta_genCounts, c_tpi0_genCounts, c_teta_recCounts, c_tpi0_recCounts: " << c_teta_genCounts << ", " << c_tpi0_genCounts << ", " << c_teta_recCounts << ", " << c_tpi0_recCounts << endl;
		if (!skipCalc) {
			efficiencies_eta.push_back( c_teta_recCounts/c_teta_genCounts); 
			efficiencies_pi0.push_back( c_tpi0_recCounts/c_tpi0_genCounts); 
			// the error of poisson counting is 1/sqrt(N) which will be used for the gen counts. For the acc we will use the fitted error. Error propagating
			efficiencies_error_eta.push_back( abs(efficiencies_eta[j])*sqrt(TMath::Power(c_teta_recCountErrors/c_teta_recCounts,2)+1/c_teta_genCounts) ); 
			efficiencies_error_pi0.push_back( abs(efficiencies_pi0[j])*sqrt(TMath::Power(c_tpi0_recCountErrors/c_tpi0_recCounts,2)+1/c_tpi0_genCounts) ); 

			if ( efficiencies_eta[j] > maxEfficiency ) { maxEfficiency = efficiencies_eta[j]; }
			if ( efficiencies_pi0[j] > maxEfficiency ) { maxEfficiency = efficiencies_pi0[j]; }
		}
		
		int massBin = j/num_tBins;
		int tBin = j%num_tBins;
		cout << "\tETA - Filling Mass Bin: " << massBin << " at t Bin: " << tBin << " with efficiency,error: " << efficiencies_eta[j] << "," << efficiencies_error_eta[j] << endl;
		cout << "\tPI0 - Filling Mass Bin: " << massBin << " at t Bin: " << tBin << " with efficiency,error: " << efficiencies_pi0[j] << "," << efficiencies_error_pi0[j] << endl;
		hist_efficiencies_pi0[massBin]->SetBinContent( tBin+1, efficiencies_pi0[j]);
		hist_efficiencies_pi0[massBin]->SetBinError( tBin+1, efficiencies_error_pi0[j]);
		hist_efficiencies_eta[massBin]->SetBinContent( tBin+1,  efficiencies_eta[j]);
		hist_efficiencies_eta[massBin]->SetBinError( tBin+1, efficiencies_error_eta[j]);
		cout << "Max Efficiency is: " << maxEfficiency << endl;
	} // closes the numHist loop
	// to hold both teta and tpi0 efficiency vectors
	efficiencies.push_back(efficiencies_eta);
	efficiencies.push_back(efficiencies_pi0);
	efficiencies_error.push_back(efficiencies_error_eta);
	efficiencies_error.push_back(efficiencies_error_pi0);

	// *********************************************
	// *********** EFFICIENCY CORRECT THE YIELDS ***************
	// *********************************************
	double maxYield=DBL_MIN;
	double minYield=DBL_MAX;
	double maxUnscaledYield=DBL_MIN;
	int branchIdx=-1;
	int minNumberOfEvents=20;
	for ( string branchName: branchNames ) {
		++branchIdx;
		for(int i=0; i<numHists; ++i){
			double unscaledYield = unscaledYields[0][branchIdx][i];
			double unscaledYieldError = unscaledYieldErrors[0][branchIdx][i];
			// Note! yieldErrorTerm and efficErrorTerm are just TERMS of an expression. Not an actual error.....
			double yieldErrorTerm=(unscaledYieldError/unscaledYield)*(unscaledYieldError/unscaledYield);
			double efficErrorTerm=(efficiencies_error[branchIdx][i]/efficiencies[branchIdx][i])*(efficiencies_error[branchIdx][i]/efficiencies[branchIdx][i]);
			yields[branchIdx][i]=unscaledYield/efficiencies[branchIdx][i];
			yieldErrors[branchIdx][i]=abs(yields[branchIdx][i])*sqrt(yieldErrorTerm+efficErrorTerm);

			// put this here so we can still properly save the diagnostic plots. Even though we will waste a fit calculation....
			if ( entriesInMetaMpi0[0][branchIdx][i] < minNumberOfEvents || entriesInMetaMpi0[1][branchIdx][i] < minNumberOfEvents ) { 
				yields[branchIdx][i]=0;
				yieldErrors[branchIdx][i]=0;
				// we only set the "data" part of unscaledYields to zero. Lets leave the "reco" part for now.
				unscaledYields[0][branchIdx][i]=0;
				unscaledYieldErrors[0][branchIdx][i]=0;
				logFile << i << ", " << "Not enough events, only {#eventsData, #eventsReco}: {" << entriesInMetaMpi0[0][branchIdx][i] << ", " << entriesInMetaMpi0[1][branchIdx][i] << "} with a minimum of " << minNumberOfEvents <<  endl; 
				continue;
			}

			cout << "Yield, error = " << unscaledYield/inverseDifferentialArea << ", " << unscaledYieldError/inverseDifferentialArea;
			cout << "   Scaled Yield, Scaled error = " << unscaledYield << ", " << unscaledYieldError;
			cout << "   EffCorrected Yield, Scaled error = " << yields[branchIdx][i] << ", " << yieldErrors[branchIdx][i];
			cout << "   Efficiency was = " << efficiencies[branchIdx][i] << endl;
			if (maxYield<yields[branchIdx][i]){
				maxYield=yields[branchIdx][i];
			}
			if (minYield>yields[branchIdx][i]){
				minYield=yields[branchIdx][i];
			}
			if (maxUnscaledYield<unscaledYields[0][branchIdx][i]){
				maxUnscaledYield=unscaledYields[0][branchIdx][i];
			}
		}
		cout << "    maxYield = " << maxYield << endl;
		cout << "    minYield = " << minYield << endl;
		if ( minYield < 0 ) {
			cout << "NEGATIVE MIN YIELD NOT GOOD. FIX ME" << endl;
		}
		cout << "    maxUnscaledYield = " << maxUnscaledYield << endl;
	}
	

	// **************************************************************************************
	// ********* NOW WE CAN START BUILDING OUT teta/tpi0 MASS HISTOGRAMS *********************
	// **************************************************************************************
	//std::vector<double> etaFitMin;
	//std::vector<double> etaFitMax;
	//std::vector<double> pi0FitMin;
	//std::vector<double> pi0FitMax; 

	//// ranges for t<0.5
	//if (selectionString == "tLT05") {
	//double etaFitMin[num_massBins] = {tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin};
	//double etaFitMax[num_massBins] = {tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-6*tStep,tMax-4*tStep,tMax-5*tStep,tMax-7*tStep,tMax-8*tStep,tMax-9*tStep,tMax-10*tStep,tMax-11*tStep,tMin};
	//double pi0FitMin[num_massBins] = {tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin,tMin};
	//double pi0FitMax[num_massBins] = {tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax,tMax,tMax-4*tStep,tMax-7*tStep,tMax-8*tStep,tMax-7*tStep,tMax-9*tStep,tMax-11*tStep,tMin};
	//}
	//// ranges for 0.5<t<1
	//else if (selectionString == "tGT05LT1") { 
	//double etaFitMin[num_massBins] = {tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep};
	//double etaFitMax[num_massBins] = {tMax,tMax,tMax,tMax,tMax,tMax-4*tStep,tMax-5*tStep,tMax-5*tStep,tMax-5*tStep,tMax-7*tStep,tMax-8*tStep,tMax-8*tStep};
	//double pi0FitMin[num_massBins] = {tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep};
	//double pi0FitMax[num_massBins] = {tMax,tMax,tMax-2*tStep,tMax,tMax,tMax-2*tStep,tMax-3*tStep,tMax-3*tStep,tMax-5*tStep,tMax-5*tStep,tMax-5*tStep,tMax-6*tStep};
	//}
	//// ranges for all
	//else if ( selectionString == "tAll"){
	//	etaFitMin = {tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep};
	//	etaFitMax = {tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep};
	//	pi0FitMin = {tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep};
	//	pi0FitMax = {tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep};
	//}
	double etaFitMin[num_massBins] = {tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep};
	double etaFitMax[num_massBins] = {tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep};
	double pi0FitMin[num_massBins] = {tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep,tMin+tStep};
	double pi0FitMax[num_massBins] = {tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep,tMax-2*tStep};

	std::vector<TF1> expPlusLinFits_eta;
	std::vector<TF1> expPlusLinFits_pi0;
	std::vector<double> tSlopes_eta;
	std::vector<double> tSlopes_pi0;
	std::vector<double> tSlopesError_eta;
	std::vector<double> tSlopesError_pi0;


	// *********************************************
	// *********** FIT THE t-distributions AND DRAW ***************
	// *********************************************

	branchIdx=-1;
	int skipLastN;
	for ( string branchName: branchNames ) {
		++branchIdx;
		std::vector<string> massBinTitles; 
		for (Int_t i=0; i<num_massBins; ++i){
			massBinTitles.push_back( std::to_string( mMin+i*mStep )+" < M(#pi^{0}#eta) < " + std::to_string( mMin+(i+1)*mStep ) );
		}
		TH1F* massHist;
		TH1F* unscaledMassHist;
		TF1* expFit_t;
		gStyle->SetErrorX(0.000001); // remove the x-error bars


		if ( selectionString == "tLT05" ){
			skipLastN = 1;
		}	
		else {
			skipLastN = 0;
		}	
		
		for (Int_t iMass=0; iMass<massBinTitles.size(); ++iMass){
			cout << "Plotting histogram for mass bin: " << iMass << endl;
			massHist = new TH1F(("massHist_"+std::to_string(iMass)+"_"+std::to_string(branchIdx)).c_str(),(massBinTitles[iMass]+";t_{#eta}/t_{#pi^{0}} GeV^{2}").c_str(),num_tBins,tMin,tMax);	
			unscaledMassHist = new TH1F(("unscaledMassHist_"+std::to_string(iMass)+"_"+std::to_string(branchIdx)).c_str(),(massBinTitles[iMass]+";t_{#eta}/t_{#pi^{0}} GeV^{2}").c_str(),num_tBins,tMin,tMax);	
			// I think we have to zero suppress our graphs to make the scaling on the second axis to make sense
			massHist->SetAxisRange(1,maxYield*1.1,"Y");

			gStyle->SetTitleSize(0.1,"t");
			massHist->GetXaxis()->SetLabelSize(0.06);
			massHist->GetYaxis()->SetLabelSize(0.06);
			unscaledMassHist->SetAxisRange(1,maxUnscaledYield*1.1,"Y");
			unscaledMassHist->SetTitleSize(1,"t");
			unscaledMassHist->GetXaxis()->SetLabelSize(0.06);
			unscaledMassHist->GetYaxis()->SetLabelSize(0.06);
			if ( branchName == "mandelstam_teta_meas") { 
				expFit_t = new TF1("expFit_t","expo",etaFitMin[iMass],etaFitMax[iMass]);//+pol0(2)) ,tMin+tStep,tMax);
			}
			else if ( branchName == "mandelstam_tpi0_meas") { 
				expFit_t = new TF1("expFit_t","expo",pi0FitMin[iMass],pi0FitMax[iMass]);//+pol0(2)) ,tMin+tStep,tMax);
			}
			expFit_t->SetLineStyle(2);
			for (Int_t i=0; i<num_tBins; ++i){
				massHist->SetBinContent( i+1, yields[branchIdx][i+num_tBins*iMass]  );
				massHist->SetBinError( i+1, yieldErrors[branchIdx][i+num_tBins*iMass]  );
				unscaledMassHist->SetBinContent( i+1, unscaledYields[0][branchIdx][i+num_tBins*iMass]  );
				unscaledMassHist->SetBinError( i+1, unscaledYieldErrors[0][branchIdx][i+num_tBins*iMass]  );
				cout << i << ") efficiency corrected yield in bin: " << i+num_tBins*iMass << " is " << yields[branchIdx][i+num_tBins*iMass] << endl;
			}
			cout << "Fitting mass hist" << endl;
			massHist->Fit("expFit_t","RBNQ");

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
				massHist->SetAxisRange(30,maxYield*10,"Y");
				massHist->GetXaxis()->SetTitleSize(0.08);
				massHist->GetYaxis()->SetTitleSize(0.08);
				allCanvases_tSlope->Update();
				expFit_t->SetLineColor(kBlue);
				expFit_t->Draw("SAME");

				allCanvases_yields->cd( iMass+1 );	
				gStyle->SetOptStat(0);
				massHist->SetMarkerColor(kBlue);
				massHist->Draw("E1 PMC");
				massHist->SetAxisRange(1,maxYield*1.1,"Y");
				massHist->GetXaxis()->SetTitleSize(0.08);
				massHist->GetYaxis()->SetTitleSize(0.08);
				allCanvases_yields->Update();
				scaleAxis = gPad->GetUymax()/(maxEfficiency*1.1); // this should be the same for all pads since I set the axisRange above
				hist_efficiencies_eta[iMass]->Scale(scaleAxis);  
				hist_efficiencies_eta[iMass]->SetLineColor(kBlue);
				hist_efficiencies_eta[iMass]->Draw("E SAME");
				hist_efficiencies_eta[iMass]->Draw("SAME Lhist");

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
				massHist->SetAxisRange(30,maxYield*10,"Y");
				massHist->SetMarkerColor(kRed);
				massHist->Draw("E1 PMC SAME");
				allCanvases_tSlope->Update();
				expFit_t->SetLineColor(kRed);
				expFit_t->Draw("SAME");

				allCanvases_yields->cd( iMass+1 );	
				massHist->SetAxisRange(1,maxYield*1.1,"Y");
				massHist->SetMarkerColor(kRed);
				massHist->Draw("E1 PMC SAME");
				hist_efficiencies_pi0[iMass]->Scale(scaleAxis);  
				hist_efficiencies_pi0[iMass]->SetLineColor(kRed);
				hist_efficiencies_pi0[iMass]->Draw("E SAME");
				hist_efficiencies_pi0[iMass]->Draw("SAME Lhist");

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

		cout << "ALL THE (YIELDS, EFFICIENCIES) FOR " << branchName << ": " <<endl;
		for ( int iHist=0; iHist<numHists; ++iHist ) {
			cout << "(" << yields[branchIdx][iHist] << ", " << efficiencies[branchIdx][iHist] << "), ";
		}
		cout << endl;
	logFile << "\n\n\n======================= END OF THIS ITERATION ===============================\n" << endl;
    	}

	// PLOTTING THE YIELD HISTOGRAMS AND THE t-Slope FITS
	 TLegend *leg;
	 for (int iMass=0; iMass<num_massBins-skipLastN; ++iMass){
	 	allCanvases_tSlope->cd(iMass+1);
		if ( selectionString == "tLT05" ) {
	 		leg = new TLegend(0.65,0.11,0.9,0.35);
		}
		else{
	 		leg = new TLegend(0.1,0.11,0.35,0.35);
		}
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


	 // Setting up the histograms for the t-slopes Vs Mpi0eta before I draw them
	 TFile *tSlopesFile = new TFile(("deckPlots/"+selectionString+"/tSlopesVsMpi0eta.root").c_str(),"RECREATE");
	 allCanvases->Clear();
	 allCanvases->cd();
	 TH1F* hist_tSlopes_eta = new TH1F("hist_tSlope_eta","; M(#pi^{0}#eta) (GeV); t-slope",num_massBins,mMin,mMax);
	 TH1F* hist_tSlopes_pi0 = new TH1F("hist_tSlope_pi0","; M(#pi^{0}#eta) (GeV); t-slope",num_massBins,mMin,mMax);
	 TH1F* hist_tSlopes_eta_scaledMpi0eta = new TH1F("hist_tSlope_eta_scaledMpi0eta","; M(#pi^{0}#eta) (GeV); t-slope",num_massBins,mMin,mMax);
	 TH1F* hist_tSlopes_pi0_scaledMpi0eta = new TH1F("hist_tSlope_pi0_scaledMpi0eta","; M(#pi^{0}#eta) (GeV); t-slope",num_massBins,mMin,mMax);
	 double logMass;
	 for ( int iMass=0; iMass < num_massBins-skipLastN; ++iMass ){
		cout << "Mass at center of massBin " << to_string(iMass) << " is: " << mMin+(iMass+0.5)*mStep << endl;
		logMass = TMath::Log(mMin+(iMass+0.5)*mStep); 
		cout << "   log of above mass is: " << logMass << endl;
		cout << "tSlope_eta for massBin " << to_string(iMass) << " is: " << tSlopes_eta[iMass] << endl;
		cout << "tSlope_pi0 for massBin " << to_string(iMass) << " is: " << tSlopes_pi0[iMass] << endl;
		cout << "tSlope_eta_scaledMpi0eta for massBin " << to_string(iMass) << " is: " << tSlopes_eta[iMass]/logMass << endl;
		cout << "tSlope_pi0_scaledMpi0eta for massBin " << to_string(iMass) << " is: " << tSlopes_pi0[iMass]/logMass << endl;
	 	hist_tSlopes_eta->SetBinContent( iMass+1, tSlopes_eta[iMass]);
	 	hist_tSlopes_pi0->SetBinContent( iMass+1, tSlopes_pi0[iMass]);
	 	hist_tSlopes_eta->SetBinError( iMass+1, tSlopesError_eta[iMass]);
	 	hist_tSlopes_pi0->SetBinError( iMass+1, tSlopesError_pi0[iMass]);

		hist_tSlopes_eta_scaledMpi0eta->SetBinContent( iMass+1, tSlopes_eta[iMass]/logMass );
		hist_tSlopes_pi0_scaledMpi0eta->SetBinContent( iMass+1, tSlopes_pi0[iMass]/logMass );
	 	hist_tSlopes_eta_scaledMpi0eta->SetBinError( iMass+1, tSlopesError_eta[iMass]/logMass );
	 	hist_tSlopes_pi0_scaledMpi0eta->SetBinError( iMass+1, tSlopesError_pi0[iMass]/logMass );
	 }

	 // Plotting the unscaled t-slopes first
	 leg = new TLegend(0.65,0.11,0.9,0.35);
	 leg->AddEntry(hist_tSlopes_eta,"t_{#eta}");
	 hist_tSlopes_eta->SetMarkerStyle(kFullCircle);
	 hist_tSlopes_eta->SetMarkerSize(1.5);
	 hist_tSlopes_eta->SetMarkerColor(kBlue);
	 hist_tSlopes_eta->GetXaxis()->SetLabelSize(0.05);
	 hist_tSlopes_eta->GetYaxis()->SetLabelSize(0.05);
	 hist_tSlopes_eta->SetMinimum(0);
	 hist_tSlopes_eta->Draw("E1 PMC");

	 hist_tSlopes_pi0->SetMarkerStyle(kFullCircle);
	 leg->AddEntry(hist_tSlopes_pi0,"t_{#pi^{0}}");
	 hist_tSlopes_pi0->SetMarkerSize(1.5);
	 hist_tSlopes_pi0->SetMarkerColor(kRed);
	 hist_tSlopes_pi0->SetMinimum(0);
	 hist_tSlopes_pi0->Draw("E1 PMC SAME");
	 leg->Draw();
	 allCanvases->SaveAs(("deckPlots/"+selectionString+"/tSlopesVsMpi0eta.png").c_str());
	 allCanvases->SaveAs(("deckPlots/"+selectionString+"/tSlopesVsMpi0eta.C").c_str());
	 allCanvases->Write();

	 // now we finally plot the scaled t-slopes. Scaled by log(Mpi0eta)
	 allCanvases->Clear();
	 leg = new TLegend(0.65,0.11,0.9,0.35);
	 leg->AddEntry(hist_tSlopes_eta_scaledMpi0eta,"t_{#eta}/log(M(#pi^{0}#eta))");
	 hist_tSlopes_eta_scaledMpi0eta->SetMarkerStyle(kFullCircle);
	 hist_tSlopes_eta_scaledMpi0eta->SetMarkerSize(1.5);
	 hist_tSlopes_eta_scaledMpi0eta->SetMarkerColor(kBlue);
	 hist_tSlopes_eta_scaledMpi0eta->GetXaxis()->SetLabelSize(0.05);
	 hist_tSlopes_eta_scaledMpi0eta->GetYaxis()->SetLabelSize(0.05);
	 hist_tSlopes_eta_scaledMpi0eta->SetMinimum(0);
	 hist_tSlopes_eta_scaledMpi0eta->Draw("E1 PMC");

	 hist_tSlopes_pi0_scaledMpi0eta->SetMarkerStyle(kFullCircle);
	 leg->AddEntry(hist_tSlopes_pi0_scaledMpi0eta,"t_{#pi^{0}}/log(M(#pi^{0}#eta))");
	 hist_tSlopes_pi0_scaledMpi0eta->SetMarkerSize(1.5);
	 hist_tSlopes_pi0_scaledMpi0eta->SetMarkerColor(kRed);
	 hist_tSlopes_pi0_scaledMpi0eta->SetMinimum(0);
	 hist_tSlopes_pi0_scaledMpi0eta->Draw("E1 PMC SAME");
	 leg->Draw();
	 allCanvases->SaveAs(("deckPlots/"+selectionString+"/tSlopesScaledVsMpi0eta.png").c_str());
	 allCanvases->SaveAs(("deckPlots/"+selectionString+"/tSlopesScaledVsMpi0eta.C").c_str());
	 allCanvases->Write();

}


void makeDeckPlots(){
	gSystem->Exec("rm -rf deckPlots");
	gSystem->Exec("mkdir -p deckPlots/tAll/mandelstam_teta_meas/data");
	gSystem->Exec("mkdir -p deckPlots/tAll/mandelstam_tpi0_meas/data");
	gSystem->Exec("mkdir -p deckPlots/tAll/tSlope/data/reco");
	gSystem->Exec("mkdir -p deckPlots/tAll/mandelstam_teta_meas/reco");
	gSystem->Exec("mkdir -p deckPlots/tAll/mandelstam_tpi0_meas/reco");
	gSystem->Exec("mkdir -p deckPlots/tAll/tSlope/reco");
	makeDeckPlot("tAll");

	gSystem->Exec("mkdir -p deckPlots/tGT1/mandelstam_teta_meas/data");
	gSystem->Exec("mkdir -p deckPlots/tGT1/mandelstam_tpi0_meas/data");
	gSystem->Exec("mkdir -p deckPlots/tGT1/tSlope/data");
	gSystem->Exec("mkdir -p deckPlots/tGT1/mandelstam_teta_meas/reco");
	gSystem->Exec("mkdir -p deckPlots/tGT1/mandelstam_tpi0_meas/reco");
	gSystem->Exec("mkdir -p deckPlots/tGT1/tSlope/reco");
	makeDeckPlot("tGT1");

	gSystem->Exec("mkdir -p deckPlots/tLT05/mandelstam_teta_meas/data");
	gSystem->Exec("mkdir -p deckPlots/tLT05/mandelstam_tpi0_meas/data");
	gSystem->Exec("mkdir -p deckPlots/tLT05/tSlope/data");
	gSystem->Exec("mkdir -p deckPlots/tLT05/mandelstam_teta_meas/reco");
	gSystem->Exec("mkdir -p deckPlots/tLT05/mandelstam_tpi0_meas/reco");
	gSystem->Exec("mkdir -p deckPlots/tLT05/tSlope/reco");
	makeDeckPlot("tLT05");

	gSystem->Exec("mkdir -p deckPlots/tGT05LT1/mandelstam_teta_meas/data");
	gSystem->Exec("mkdir -p deckPlots/tGT05LT1/mandelstam_tpi0_meas/data");
	gSystem->Exec("mkdir -p deckPlots/tGT05LT1/tSlope/data");
	gSystem->Exec("mkdir -p deckPlots/tGT05LT1/mandelstam_teta_meas/reco");
	gSystem->Exec("mkdir -p deckPlots/tGT05LT1/mandelstam_tpi0_meas/reco");
	gSystem->Exec("mkdir -p deckPlots/tGT05LT1/tSlope/reco");
	makeDeckPlot("tGT05LT1");
}
