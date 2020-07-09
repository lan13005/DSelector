//Looking at teta/tpi and bins of Metap and Mpi0p which is s23 in vincent/colin language


double degToRad=TMath::Pi()/180;
// par[3] is used to shift phase by the para or perp orientation, either 0 for para or 90 for perp. 0/-45 is para and 45/90 is perp. 
int numDOFsig_sc = 2;
Double_t shiftedCos000(Double_t *x, Double_t *par){
	return par[0]*(1.0 - par[1]*TMath::Cos(2*degToRad*x[0]));
}
Double_t shiftedCos045(Double_t *x, Double_t *par){
	return par[0]*(1.0 + par[1]*TMath::Cos(2*degToRad*(x[0]-45-90)));
}
Double_t shiftedCos090(Double_t *x, Double_t *par){
	return par[0]*(1.0 - par[1]*TMath::Cos(2*degToRad*(x[0]-90)));
}
Double_t shiftedCos135(Double_t *x, Double_t *par){
	return par[0]*(1.0 - par[1]*TMath::Cos(2*degToRad*(x[0]-(-45))));
}
int numDOFsig_scAMO=3;
Double_t shiftedCosAMO(Double_t *x, Double_t *par){
	return par[0]*(1.0 - par[1]*TMath::Cos(2*degToRad*(x[0]-par[2])));
}

int numDOFsig_flat = 1;
Double_t flat(Double_t *x, Double_t *par){
	return par[0];
}

int numDOFsig_asym=4;
Double_t asymmetry(Double_t *x, Double_t *par){
	return ((par[0]+par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3]))/(2+(par[0]-par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3]))));
}

void fitAsymmetryPlots_baryon(){
        // Do some cleaning
        gSystem->Exec("rm -rf asymmetryPlots/SigVst1_binS23");
        gSystem->Exec("mkdir asymmetryPlots/SigVst1_binS23");
        //////////////////

	static const int nSetsBS=1; // (number of sets to bootstrap-1). 1 is for not bootstrapping at all
	int maxPrintBS=1; // Only show up to this value when going through nSetsBS. This includes the "full" set in the count

	// seems like this this would somehow instantiate gMinuit maybe? If I dont do this I get some errors
	//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // If I dont get this I get the error "Error in <TMinuitMinimizer::Scan>:  Error executing command SCAN" 
	static const int num_tBins=5;
	static const int num_mBins=5;
	string phis_orientation[5] =  {"phi000", "phi045", "phi090", "phi135", "phiAMO"};

	gStyle->SetOptFit();
	gStyle->SetStatY(1);
	gStyle->SetStatX(1);
	gStyle->SetStatW(0.16);
	gStyle->SetStatH(0.16);

	TCanvas *allCanvases = new TCanvas("","",1440,900);


	// *****************************
	// Define flux ratios for 2017, 2018_1, 2018_2 
	// *****************************
	static const int nDataSets = 3;
	double fluxRatios_90_0[nDataSets] = {  4.346818e+12/4.188001e+12, 0.965429, 0.918503 };
	double fluxRatios_45_135[nDataSets] = {  4.076065e+12/4.095013e+12, 1.02261, 1.03254 };
	string dataSetTag[nDataSets] = { "data_2017", "data_2018_1", "data_2018_8" };
	string dataFolders[nDataSets] = {"deg000_data_2017","deg000_data_2018_1", "deg000_data_2018_8"};

	static const int nTagEta = 1;
	string tagEta[nTagEta] = {""};
	string tagPi0[nTagEta] = {""};
	string tag[nTagEta] = {""};// same as above but a single name to group them both.
	std::vector<int> nEventsPhiEta[nTagEta]; // second one is for the vanHove selected
	std::vector<int> nEventsPhiPi0[nTagEta];

	// *****************************
	// Loading histograms for  2017, 2018_1, 2018_2, then scaling yield by the flux ratio. Then we can finally add the data sets together
	// *****************************
	//
	//  ************* these are not really used, only as a place holder before scaling them to get the total *************
	
	// There are 10 histograms as shown below
	static const int nHists=10;
	cout << "Shape of the tensor of histograms: ["<<nTagEta<<"]["<<num_tBins<<"]["<<num_mBins<<"]["<<nDataSets<<"]["<<nSetsBS<<"]"<< endl;
	TH1F *phi000_eta_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS]; // 2 for the tag and 3 for the datasets
	TH1F *phi045_eta_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS]; 
	TH1F *phi090_eta_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS];
	TH1F *phi135_eta_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS];
	TH1F *phiAMO_eta_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS];
	TH1F *phi000_pi0_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS];
	TH1F *phi045_pi0_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS];
	TH1F *phi090_pi0_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS];
	TH1F *phi135_pi0_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS];
	TH1F *phiAMO_pi0_unscaled[nTagEta][num_mBins][num_tBins][nDataSets][nSetsBS];
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// ************************************************************ LOAD THE DATA WHATEVER WAY YOU WANT **********************************************************************
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//
	// *************************************
	// Loading the data in from the histogram 
	// *************************************
	if (nSetsBS == 1){
		for (int iSet=0; iSet<nSetsBS; ++iSet){ // we add 1 since using the full dataset takes the first element of the array
			for (int iData=0; iData <nDataSets; ++iData){
				string dataFileName = "/d/grid15/ln16/pi0eta/092419/degALL_"+dataSetTag[iData]+"_BA_hists_DSelector.root";
				//string dataFileName = "/d/grid15/ln16/pi0eta/092419/newGraphs_histValues/rootFiles/deg000_data_"+dataSetTag[iData]+"_hists_DSelector.root";
				TFile *dataFile = new TFile(dataFileName.c_str());
				cout << "LOADING ROOT FILE: " << dataFileName << endl; 
				for (int iTag=0; iTag < sizeof(tagEta)/sizeof(tagEta[0]); ++iTag){
					for (int iMass=0; iMass<num_mBins; ++iMass){
						for (int it=0; it<num_tBins; ++it){
							string affixEta = "Bin"+to_string(iMass)+"_tetaBin"+to_string(it);
							string affixPi0 = "Bin"+to_string(iMass)+"_tpi0Bin"+to_string(it);
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_045_Mpi0p"+affixEta).c_str(), phi045_eta_unscaled[iTag][iMass][it][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_045_Metap"+affixPi0).c_str(), phi045_pi0_unscaled[iTag][iMass][it][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_090_Mpi0p"+affixEta).c_str(), phi090_eta_unscaled[iTag][iMass][it][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_090_Metap"+affixPi0).c_str(), phi090_pi0_unscaled[iTag][iMass][it][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_AMO_Mpi0p"+affixEta).c_str(), phiAMO_eta_unscaled[iTag][iMass][it][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_AMO_Metap"+affixPi0).c_str(), phiAMO_pi0_unscaled[iTag][iMass][it][iData][iSet]);
							// Need to scale these para yields by flux ratio.
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_000_Mpi0p"+affixEta).c_str(), phi000_eta_unscaled[iTag][iMass][it][iData][iSet]);
							cout << "phi000_eta_unscaled["<<iTag<<"]["<<iMass<<"]["<<it<<"]["<<iData<<"]["<<iSet<<"]" << endl;
							cout << (" -- prodPlanePSphi"+tagEta[iTag]+"_000_Mpi0p"+affixEta).c_str() << endl;
							cout << " -- nentries=" << phi000_eta_unscaled[iTag][iMass][it][iData][iSet]->GetEntries() << endl;
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_000_Metap"+affixPi0).c_str(), phi000_pi0_unscaled[iTag][iMass][it][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_135_Mpi0p"+affixEta).c_str(), phi135_eta_unscaled[iTag][iMass][it][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_135_Metap"+affixPi0).c_str(), phi135_pi0_unscaled[iTag][iMass][it][iData][iSet]);
						}
					}
				}
			}
		}
	}
	
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// ************************************************** NOW THAT WE LOADED THE DATA WE CAN SCALE IT AND FIT IT *************************************************************
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------------------------------------------------------------------------------------------
	double asymmetries_000_eta[nTagEta][nSetsBS][num_mBins][num_tBins];
	double asymmetries_045_eta[nTagEta][nSetsBS][num_mBins][num_tBins];
	double asymmetries_000_eta_err[nTagEta][nSetsBS][num_mBins][num_tBins];
	double asymmetries_045_eta_err[nTagEta][nSetsBS][num_mBins][num_tBins];
	double asymmetries_000_pi0[nTagEta][nSetsBS][num_mBins][num_tBins];
	double asymmetries_045_pi0[nTagEta][nSetsBS][num_mBins][num_tBins];
	double asymmetries_000_pi0_err[nTagEta][nSetsBS][num_mBins][num_tBins];
	double asymmetries_045_pi0_err[nTagEta][nSetsBS][num_mBins][num_tBins];
	cout << "--------------------------------------------\nLOADING AND SCALING HISTOGRAMS\n--------------------------------------------" << endl;
	cout << "As a simple test we can GetMaximum before and after the scaling to see if it actually worked" << endl;
	cout << "\tWe will follow phi000_eta_total throughout the process" << endl;
	for (int iSet=0; iSet <nSetsBS; ++iSet) {
		// These are the weighted summed histograms. weighted according to the Flux ratio
		TH1F *phi000_eta_total[nTagEta][num_mBins][num_tBins];
		TH1F *phi045_eta_total[nTagEta][num_mBins][num_tBins];
		TH1F *phi090_eta_total[nTagEta][num_mBins][num_tBins];
		TH1F *phi135_eta_total[nTagEta][num_mBins][num_tBins];
		TH1F *phiAMO_eta_total[nTagEta][num_mBins][num_tBins];
		TH1F *phi000_pi0_total[nTagEta][num_mBins][num_tBins];
		TH1F *phi045_pi0_total[nTagEta][num_mBins][num_tBins];
		TH1F *phi090_pi0_total[nTagEta][num_mBins][num_tBins];
		TH1F *phi135_pi0_total[nTagEta][num_mBins][num_tBins];
		TH1F *phiAMO_pi0_total[nTagEta][num_mBins][num_tBins];
		// Now we sum the scaled histograms over the 3 datasets
		double maximum_before;
		double maximum_after;
		int totalEntriesInPhi000_eta_total;
		for (int iTag=0; iTag < nTagEta; ++iTag) {
			for (int iMass=0; iMass<num_mBins; ++iMass){
				for (int it=0; it<num_tBins; ++it){
					totalEntriesInPhi000_eta_total=0;
					cout << "--------------\n" << tagEta[0] << "\n--------------" << endl;
					// We first have to clone the first dataTag (2017 run) and scale by flux ratio if in para config
					// we will scale the total in this to track the number of maxima and entires before and after scaling
					cout << "phi000_eta_unscaled["<<iTag<<"]["<<iMass<<"]["<<it<<"][0]["<<iSet<<"]" << endl;
					phi000_eta_total[iTag][iMass][it] = (TH1F *)phi000_eta_unscaled[iTag][iMass][it][0][iSet]->Clone();
					phi000_pi0_total[iTag][iMass][it] = (TH1F *)phi000_pi0_unscaled[iTag][iMass][it][0][iSet]->Clone();
					phi135_eta_total[iTag][iMass][it] = (TH1F *)phi135_eta_unscaled[iTag][iMass][it][0][iSet]->Clone();
					phi135_pi0_total[iTag][iMass][it] = (TH1F *)phi135_pi0_unscaled[iTag][iMass][it][0][iSet]->Clone();
					maximum_before = phi000_eta_total[iTag][iMass][it]->GetMaximum();
					cout << "-Maximum before scaling of phi000_eta_total[" << iTag << "][" << it << "]: " << maximum_before << endl;
					cout << "-NEntries before scaling of phi000_eta_total[" << iTag << "]" << it << "]: " << phi000_eta_total[iTag][iMass][it]->GetEntries() << endl;
					phi000_eta_total[iTag][iMass][it]->Scale(fluxRatios_90_0[0]);
					phi000_pi0_total[iTag][iMass][it]->Scale(fluxRatios_90_0[0]);
					phi135_eta_total[iTag][iMass][it]->Scale(fluxRatios_45_135[0]);
					phi135_pi0_total[iTag][iMass][it]->Scale(fluxRatios_45_135[0]);
					maximum_after = phi000_eta_total[iTag][iMass][it]->GetMaximum();
					cout << "-Maximum after scaling of phi000_eta_total[" << iTag << "]" << it << "]: " << maximum_after << endl;
					cout << "-NEntries after scaling of phi000_eta_total[" << iTag << "]" << it << "]: " << phi000_eta_total[iTag][iMass][it]->GetEntries() << endl; 
					cout << "-Scale factor vs Expected: " << maximum_after/maximum_before << " vs " << fluxRatios_90_0[0] << endl; 
					phi045_eta_total[iTag][iMass][it] = (TH1F *)phi045_eta_unscaled[iTag][iMass][it][0][iSet]->Clone();
					phi045_pi0_total[iTag][iMass][it] = (TH1F *)phi045_pi0_unscaled[iTag][iMass][it][0][iSet]->Clone();
					phi090_eta_total[iTag][iMass][it] = (TH1F *)phi090_eta_unscaled[iTag][iMass][it][0][iSet]->Clone();
					phi090_pi0_total[iTag][iMass][it] = (TH1F *)phi090_pi0_unscaled[iTag][iMass][it][0][iSet]->Clone();
					phiAMO_eta_total[iTag][iMass][it] = (TH1F *)phiAMO_eta_unscaled[iTag][iMass][it][0][iSet]->Clone();
					phiAMO_pi0_total[iTag][iMass][it] = (TH1F *)phiAMO_pi0_unscaled[iTag][iMass][it][0][iSet]->Clone();

					totalEntriesInPhi000_eta_total += phi000_eta_total[iTag][iMass][it]->GetEntries();
					cout << "-Cumulative entries after adding dataSet=0: " << totalEntriesInPhi000_eta_total << endl;
					// for the last 2 runs we simply add with a weight if in para config and add with weight=1 if in perp config
					for (int iData=1; iData < nDataSets; ++iData){
						// keeping a count of the total entries to check if we are doing the Add right
						totalEntriesInPhi000_eta_total += phi000_eta_unscaled[iTag][iMass][it][iData][iSet]->GetEntries();
						cout << "-Cumulative entries after adding dataSet=" << iData << ": " << totalEntriesInPhi000_eta_total << endl;

						// scale the rest of the data sets before adding them to the scaled+cloned 2017 dataset
						// we will never use unscaled histograms anymore so we can just scale it directly
						phi000_eta_unscaled[iTag][iMass][it][iData][iSet]->Scale(fluxRatios_90_0[iData]);
						phi000_pi0_unscaled[iTag][iMass][it][iData][iSet]->Scale(fluxRatios_90_0[iData]);
						phi135_eta_unscaled[iTag][iMass][it][iData][iSet]->Scale(fluxRatios_45_135[iData]);
						phi135_pi0_unscaled[iTag][iMass][it][iData][iSet]->Scale(fluxRatios_45_135[iData]);
						phi000_eta_total[iTag][iMass][it]->Add(phi000_eta_unscaled[iTag][iMass][it][iData][iSet]);
        	        			phi000_pi0_total[iTag][iMass][it]->Add(phi000_pi0_unscaled[iTag][iMass][it][iData][iSet]);
        	        			phi135_eta_total[iTag][iMass][it]->Add(phi135_eta_unscaled[iTag][iMass][it][iData][iSet]);
        	        			phi135_pi0_total[iTag][iMass][it]->Add(phi135_pi0_unscaled[iTag][iMass][it][iData][iSet]);

						// Dont have to scale these guys so just add them
        	        			phi045_eta_total[iTag][iMass][it]->Add(phi045_eta_unscaled[iTag][iMass][it][iData][iSet]);
        	        			phi045_pi0_total[iTag][iMass][it]->Add(phi045_pi0_unscaled[iTag][iMass][it][iData][iSet]);
        	        			phi090_eta_total[iTag][iMass][it]->Add(phi090_eta_unscaled[iTag][iMass][it][iData][iSet]);
        	        			phi090_pi0_total[iTag][iMass][it]->Add(phi090_pi0_unscaled[iTag][iMass][it][iData][iSet]);
        	        			phiAMO_eta_total[iTag][iMass][it]->Add(phiAMO_eta_unscaled[iTag][iMass][it][iData][iSet]);
        	        			phiAMO_pi0_total[iTag][iMass][it]->Add(phiAMO_pi0_unscaled[iTag][iMass][it][iData][iSet]);
                                                cout << "Total entries in phi000_eta["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phi000_eta_total[iTag][iMass][it]->GetEntries() << endl;
                                                cout << "Total entries in phi045_eta["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phi045_eta_total[iTag][iMass][it]->GetEntries() << endl;
                                                cout << "Total entries in phi090_eta["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phi090_eta_total[iTag][iMass][it]->GetEntries() << endl;
                                                cout << "Total entries in phi135_eta["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phi135_eta_total[iTag][iMass][it]->GetEntries() << endl;
                                                cout << "Total entries in phiAMO_eta["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phiAMO_eta_total[iTag][iMass][it]->GetEntries() << endl;
                                                cout << "Total entries in phi000_pi0["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phi000_pi0_total[iTag][iMass][it]->GetEntries() << endl;
                                                cout << "Total entries in phi045_pi0["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phi045_pi0_total[iTag][iMass][it]->GetEntries() << endl;
                                                cout << "Total entries in phi090_pi0["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phi090_pi0_total[iTag][iMass][it]->GetEntries() << endl;
                                                cout << "Total entries in phi135_pi0["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phi135_pi0_total[iTag][iMass][it]->GetEntries() << endl;
                                                cout << "Total entries in phiAMO_pi0["<<iTag<<"]["<<iMass<<"]["<<it<<"]: " << phiAMO_pi0_total[iTag][iMass][it]->GetEntries() << endl;
					}
				}
			}
		}	


		cout << "--------------------------------------------\nSTARTING FITTING\n--------------------------------------------" << endl;


	        double minPi0P[num_mBins] = {1.15, 1.4, 1.6, 1.75, 1.95};
	        double maxPi0P[num_mBins] = {1.4, 1.6, 1.75, 1.95, 2.4};
	        double minEtaP[num_mBins] = {1.5, 1.65, 1.9, 2.2, 2.4};
	        double maxEtaP[num_mBins] = {1.65, 1.9, 2.2, 2.4, 2.8};
	        //double minPi0P[3] = {1.1, 1.5, 1.9};
	        //double maxPi0P[3] = {1.5, 1.9, 2.3};
	        //double minEtaP[3] = {1.5, 1.9, 2.3};
	        //double maxEtaP[3] = {1.9, 2.3, 2.7};
		for (int iMass=0; iMass<num_mBins; ++iMass){
			// *****************************
			// Defining array to hold Asymmetry values from fiting asymmetry histograms
			// *****************************
			static const int numPolarizations = 4;
			string names[4] = {"phi000","phi045","phi090","phi135"};
			//double perpOrPara[4] = {0, 90, 90, 0}; //p2
			double orientation[4] = {0, 45, 0, -45}; //p1
			double Phi0_0_90 = 0+3.1;
			double Phi0_45_135 = -45+3.2;

			// *****************************
			// Defining array to hold P*Sigma values from fitting prodPlanePhi
			// *****************************
			vector<double> phi_000_pi0; // will have a value for each tBin
			vector<double> phi_045_pi0;
			vector<double> phi_090_pi0;
			vector<double> phi_135_pi0;
			vector<double> phi_000_eta;
			vector<double> phi_045_eta;
			vector<double> phi_090_eta;
			vector<double> phi_135_eta;
			vector< vector<double > > phis_pi0_PSig; // holds all the diamond orientations
       			phis_pi0_PSig.push_back(phi_000_pi0);
			phis_pi0_PSig.push_back(phi_045_pi0);
			phis_pi0_PSig.push_back(phi_090_pi0);
			phis_pi0_PSig.push_back(phi_135_pi0);
			vector< vector<double > > phis_eta_PSig;
			phis_eta_PSig.push_back(phi_000_eta);
			phis_eta_PSig.push_back(phi_045_eta);
			phis_eta_PSig.push_back(phi_090_eta);
			phis_eta_PSig.push_back(phi_135_eta);
			vector<double> phi_000_pi0_err;
			vector<double> phi_045_pi0_err;
			vector<double> phi_090_pi0_err;
			vector<double> phi_135_pi0_err;
			vector<double> phi_000_eta_err;
			vector<double> phi_045_eta_err;
			vector<double> phi_090_eta_err;
			vector<double> phi_135_eta_err;
			vector< vector<double > > phis_pi0_PSig_err;
			phis_pi0_PSig_err.push_back(phi_000_pi0_err);
			phis_pi0_PSig_err.push_back(phi_045_pi0_err);
			phis_pi0_PSig_err.push_back(phi_090_pi0_err);
			phis_pi0_PSig_err.push_back(phi_135_pi0_err);
			vector< vector<double > > phis_eta_PSig_err;
       			phis_eta_PSig_err.push_back(phi_000_eta_err);
			phis_eta_PSig_err.push_back(phi_045_eta_err);
			phis_eta_PSig_err.push_back(phi_090_eta_err);
			phis_eta_PSig_err.push_back(phi_135_eta_err);

		        double tBins[num_tBins];
		        double tBins_err[num_tBins];
			Int_t fitStatus;
			for (int iTag=0; iTag<nTagEta; ++iTag){
				cout << "Starting iTag=" << iTag << endl;
		                for (int it=0; it<num_tBins; ++it){
			                double lowerT = it*0.2;
			                double upperT = (it+1)*0.2;
                                        double midT = (it+0.5)*0.2; 
					// Need to define the polarization with the systematic shift
					// the shift for 0/90 is 3.1 and shifting 45/135 is 3.2. Since 0 and 135/-45 is para the beginning orientation is 0 and -45 whih then we shift by 3.1, 3.2 respectively
					string fitOption = "E S";
					TH1F *phi000_eta = phi000_eta_total[iTag][iMass][it];
					TH1F *phi045_eta = phi045_eta_total[iTag][iMass][it];
					TH1F *phi090_eta = phi090_eta_total[iTag][iMass][it];
					TH1F *phi135_eta = phi135_eta_total[iTag][iMass][it];
					TH1F *phiAMO_eta = phiAMO_eta_total[iTag][iMass][it];
					TH1F *phi000_pi0 = phi000_pi0_total[iTag][iMass][it];
					TH1F *phi045_pi0 = phi045_pi0_total[iTag][iMass][it];
					TH1F *phi090_pi0 = phi090_pi0_total[iTag][iMass][it];
					TH1F *phi135_pi0 = phi135_pi0_total[iTag][iMass][it];
					TH1F *phiAMO_pi0 = phiAMO_pi0_total[iTag][iMass][it];
					TH1F *phis_eta[4] =  {phi000_eta, phi045_eta, phi090_eta, phi135_eta};
					TH1F *phis_pi0[4] =  {phi000_pi0, phi045_pi0, phi090_pi0, phi135_pi0};

					// *****************************
					// Defining the xBins which will be shared
					// *****************************
					//tBins[iMass] = 1.5+(iMass+0.5)*xBinSize;
					//tBins_err[iMass] = xBinSize/2;
					//tBins[iMass] = 1.4+(iMass+0.5)*xBinSize;
					//tBins_err[iMass] = xBinSize/2;

			 	        tBins[it] = midT;
			 	        tBins_err[it] = 0.1;


					// *****************************
					// GETTING ASYMMETRIES
					// *****************************
					// **** IF WE EVER DECIDE TO CALCULATE ASYMMETRY OURSELVES THIS SHOULD WORK I THINK **
					//TH1F* asymmetry000_090_eta = (TH1F *)phi090_eta->Clone();
					//TH1F* asymmetry000_090_eta_denom = (TH1F *)phi090_eta->Clone();
					//TH1F* asymmetry045_135_eta = (TH1F *)phi045_eta->Clone();
					//TH1F* asymmetry045_135_eta_denom = (TH1F *)phi045_eta->Clone();
					//asymmetry000_090_eta->Add(phi000_eta,-1);
					//asymmetry045_135_eta->Add(phi135_eta,-1);
					//asymmetry000_090_eta_denom->Add(phi000_eta,1);
					//asymmetry045_135_eta_denom->Add(phi135_eta,1);
					//asymmetry000_090_eta->Divide(asymmetry000_090_eta_denom);
					//asymmetry045_135_eta->Divide(asymmetry045_135_eta_denom);
					// ------------------------------------------------------------------------------------
					TH1* asymmetry000_090_eta = phi090_eta->GetAsymmetry(phi000_eta);
					TH1* asymmetry045_135_eta = phi045_eta->GetAsymmetry(phi135_eta);
					asymmetry000_090_eta->SetTitle("0/90 Asymmetry");
					asymmetry045_135_eta->SetTitle("45/135 Asymmetry");
					TH1* asymmetry000_090_pi0 = phi090_pi0->GetAsymmetry(phi000_pi0);
					TH1* asymmetry045_135_pi0 = phi045_pi0->GetAsymmetry(phi135_pi0);
					asymmetry000_090_pi0->SetTitle(("0/90 Asymmetry"+tagEta[iTag]).c_str());
					asymmetry045_135_pi0->SetTitle(("45/135 Asymmetry"+tagPi0[iTag]).c_str());

					// *****************************
					// Fitting asymmetry for eta
					// *****************************
					allCanvases->Clear();
					allCanvases->Divide(2,1);
					allCanvases->cd(1);
					// We set P_perp = P_para = 0.35 which is close to the expected. We initialize asymmetry to be 0 since it can vary from [-1,1]. 
					TF1 * fit_asym = new TF1("fit_asym",asymmetry,-180,180,numDOFsig_asym); 
					fit_asym->SetParameters(0.35,0.35,0.5,Phi0_0_90);
					fit_asym->FixParameter(0,0.35);
					fit_asym->FixParameter(1,0.35);
					fit_asym->SetParLimits(2,-1,1);
					fit_asym->FixParameter(3,Phi0_0_90);

					TFitResultPtr fitPointer = asymmetry000_090_eta->Fit(fit_asym,fitOption.c_str());
					asymmetries_000_eta[iTag][iSet][iMass][it] = fit_asym->GetParameter(2);
					asymmetries_000_eta_err[iTag][iSet][iMass][it] = fit_asym->GetParError(2);
					asymmetry000_090_eta->Draw("SAME");
					asymmetry000_090_eta->SetAxisRange(-0.5,0.6,"Y");
					TGraph* likelihoodFit_000_090_eta = new TGraph(500); 
					//fitPointer->Scan(2,likelihoodFit_000_090_eta,0,1);

					allCanvases->cd(2);
					fit_asym->SetParameters(0.35,0.35,0.5,Phi0_45_135);
					fit_asym->FixParameter(0,0.35);
					fit_asym->FixParameter(1,0.35);
					fit_asym->SetParLimits(2,-1,1);
					fit_asym->FixParameter(3,Phi0_45_135);

					fitPointer = asymmetry045_135_eta->Fit(fit_asym,fitOption.c_str());
					asymmetries_045_eta[iTag][iSet][iMass][it] = fit_asym->GetParameter(2);
					asymmetries_045_eta_err[iTag][iSet][iMass][it] = fit_asym->GetParError(2);
					asymmetry045_135_eta->Draw("SAME");
					asymmetry045_135_eta->SetAxisRange(-0.5,0.6,"Y");
					TGraph* likelihoodFit_045_135_eta = new TGraph(500); 
					//fitPointer->Scan(2,likelihoodFit_045_135_eta,0,1);

					if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVst1_binS23/asymmetry"+tagEta[iTag]+"_Mpi0pBin"+to_string(iMass)+"_tetaBin"+to_string(it)+"_iSet"+to_string(iSet)+".png").c_str()); }
                                        cout << "Fitted asym" << endl;
					
					//allCanvases->Clear();
					//allCanvases->Divide(2,1);
					//allCanvases->cd(1);
					//likelihoodFit_000_090_eta->Draw("ALP");
					//likelihoodFit_000_090_eta->SetTitle(("0/90 - tBin"+to_string(it)).c_str());
					//allCanvases->cd(2);
					//likelihoodFit_045_135_eta->Draw("ALP");
					//likelihoodFit_045_135_eta->SetTitle(("045_135 - tBin"+to_string(it)).c_str());
					//if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVst1_binS23/likelihood"+tagEta[iTag]+"_Mpi0etaBin"+to_string(iMass)+"_tetaBin"+to_string(it)+"_iSet"+to_string(iSet)+".root").c_str()); }

					// *****************************
					// Fitting asymmetry for pi0
					// *****************************
					allCanvases->Clear();
					allCanvases->Divide(2,1);
					allCanvases->cd(1);
					// We set P_perp = P_para = 0.35 which is close to the expected. We initialize asymmetry to be 0 since it can vary from [-1,1]. 
					fit_asym = new TF1("fit_asym",asymmetry,-180,180,numDOFsig_asym); 
					fit_asym->SetParameters(0.35,0.35,0.5,Phi0_0_90);
					fit_asym->FixParameter(0,0.35);
					fit_asym->FixParameter(1,0.35);
					fit_asym->SetParLimits(2,-1,1);
					fit_asym->FixParameter(3,Phi0_0_90);
					fitPointer = asymmetry000_090_pi0->Fit(fit_asym,fitOption.c_str());
					asymmetries_000_pi0[iTag][iSet][iMass][it] = fit_asym->GetParameter(2);
					asymmetries_000_pi0_err[iTag][iSet][iMass][it] = fit_asym->GetParError(2);
					asymmetry000_090_pi0->Draw("SAME");
					asymmetry000_090_pi0->SetAxisRange(-0.5,0.6,"Y");
					TGraph* likelihoodFit_000_090_pi0 = new TGraph(50); 
					//fitPointer->Scan(2,likelihoodFit_000_090_pi0,-2,2);
			
					allCanvases->cd(2);
					fit_asym->SetParameters(0.35,0.35,0.5,Phi0_45_135);
					fit_asym->FixParameter(0,0.35);
					fit_asym->FixParameter(1,0.35);
					fit_asym->SetParLimits(2,-1,1);
					fit_asym->FixParameter(3,Phi0_45_135);
					fitPointer = asymmetry045_135_pi0->Fit(fit_asym,fitOption.c_str());
					asymmetry045_135_pi0->Draw("SAME");
					asymmetry045_135_pi0->SetAxisRange(-0.5,0.6,"Y");
					asymmetries_045_pi0[iTag][iSet][iMass][it] = fit_asym->GetParameter(2);
					asymmetries_045_pi0_err[iTag][iSet][iMass][it] = fit_asym->GetParError(2);
					TGraph* likelihoodFit_045_135_pi0 = new TGraph(50); 
					//fitPointer->Scan(2,likelihoodFit_045_135_pi0,-2,2);

					if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVst1_binS23/asymmetry"+tagPi0[iTag]+"_MetapBin"+to_string(iMass)+"_tpi0Bin"+to_string(it)+"_iSet"+to_string(iSet)+".png").c_str()); }
                                        cout << "Fitted asym" << endl;

					//allCanvases->Clear();
					//allCanvases->Divide(2,1);
					//allCanvases->cd(1);
					//likelihoodFit_000_090_eta->Draw("ALP");
					//likelihoodFit_000_090_eta->SetTitle(("0/90 - tBin"+to_string(it)).c_str());
					//allCanvases->cd(2);
					//likelihoodFit_045_135_eta->Draw("ALP");
					//likelihoodFit_045_135_eta->SetTitle(("045/135 - tBin"+to_string(it)).c_str());
					//if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVst1_binS23/likelihood"+tagEta[iTag]+"_Mpi0etaBin"+to_string(iMass)+"_tetaBin"+to_string(it)+"_iSet"+to_string(iSet)+".root").c_str()); }
			
					// *****************************
					// Fitting prodPlanePhi for eta
					// *****************************
					cout << "Starting phi eta fit" << endl;
					allCanvases->Clear();
					allCanvases->Divide(3,2);
					TF1 * fit_sc000 = new TF1("fit_sc000",shiftedCos000,-180,180,numDOFsig_sc); 
					TF1 * fit_sc045 = new TF1("fit_sc045",shiftedCos045,-180,180,numDOFsig_sc); 
					TF1 * fit_sc090 = new TF1("fit_sc090",shiftedCos090,-180,180,numDOFsig_sc); 
					TF1 * fit_sc135 = new TF1("fit_sc135",shiftedCos135,-180,180,numDOFsig_sc); 
					std::vector<TF1*> fit_scs = {fit_sc000, fit_sc045, fit_sc090, fit_sc135 };
					double p0;
					int counter=0;
					for (auto phi: phis_eta){
						allCanvases->cd(counter+1);
						phi->SetTitle((names[counter]).c_str());
						nEventsPhiEta[iTag].push_back(phi->GetEntries());
						p0 = phi->GetEntries()/phi->GetNbinsX();
						fit_scs[counter]->SetParameters(p0,0.1);
						cout << "\n\nFixing orientation for phi fit for fast eta to: " << orientation[counter] << endl;
						fitStatus = phi->Fit(fit_scs[counter],"E S");
						phis_eta_PSig[counter].push_back(fit_scs[counter]->GetParameter(1));
						phis_eta_PSig_err[counter].push_back(fit_scs[counter]->GetParError(1));
						phi->Draw("SAME");
						++counter;
					}
					TF1 * fit_flat = new TF1("fit_flat",shiftedCosAMO,-180,180,numDOFsig_scAMO); 
					allCanvases->cd(5);
					phiAMO_eta->SetTitle("phiAMO");
					p0 = phiAMO_eta->GetEntries()/phiAMO_eta->GetNbinsX();
					nEventsPhiEta[iTag].push_back(phiAMO_eta->GetEntries());
					cout << "Doing flat fit to AMO for eta has entries=" << phiAMO_eta->GetEntries()  << endl;
					fitStatus = phiAMO_eta->Fit(fit_flat,"RQE");
					phiAMO_eta->Draw("SAME");
					if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVst1_binS23/phiYieldFits"+tagEta[iTag]+"_Mpi0pBin"+to_string(iMass)+"_tetaBin"+to_string(it)+"_iSet"+to_string(iSet)+".png").c_str()); }
					
			
					// *****************************
					// Fitting prodPlanePhi for pi0
					// *****************************
					cout << "Starting phi pi0 fit" << endl;
					allCanvases->Clear();
					allCanvases->Divide(3,2);
					fit_sc000 = new TF1("fit_sc000",shiftedCos000,-180,180,numDOFsig_sc); 
					fit_sc045 = new TF1("fit_sc045",shiftedCos045,-180,180,numDOFsig_sc); 
					fit_sc090 = new TF1("fit_sc090",shiftedCos090,-180,180,numDOFsig_sc); 
					fit_sc135 = new TF1("fit_sc135",shiftedCos135,-180,180,numDOFsig_sc); 
					fit_scs = {fit_sc000, fit_sc045, fit_sc090, fit_sc135 };
					counter=0;
					for (auto phi: phis_pi0){
						allCanvases->cd(counter+1);
						phi->SetTitle((names[counter]).c_str());
						p0 = phi->GetEntries()/phi->GetNbinsX();
						fit_scs[counter]->SetParameters(p0,0.1);
						cout << "Fixing orientation for phi fit for fast pi0 to: " << orientation[counter] << endl;
						fitStatus = phi->Fit(fit_scs[counter],"E S");
						phis_pi0_PSig[counter].push_back(fit_scs[counter]->GetParameter(1));
						phis_pi0_PSig_err[counter].push_back(fit_scs[counter]->GetParError(1));
						phi->Draw("SAME");
						nEventsPhiPi0[iTag].push_back(phi->GetEntries());
						cout << "(iTag=" << iTag << ")Entries in nEventsPhiPi0 if different orientations: " << phi->GetEntries() << endl;
						++counter;
					}
					fit_flat = new TF1("fit_flat",shiftedCosAMO,-180,180,numDOFsig_scAMO); 
					allCanvases->cd(5);
					phiAMO_pi0->SetTitle("phiAMO");
					nEventsPhiPi0[iTag].push_back(phiAMO_pi0->GetEntries());
					p0 = phiAMO_pi0->GetEntries()/phiAMO_pi0->GetNbinsX();
					cout << "Doing flat fit to AMO for pi0 has entries=" << phiAMO_pi0->GetEntries() << endl;
					fitStatus = phiAMO_pi0->Fit(fit_flat,"E S");
					phiAMO_pi0->Draw("SAME");
					cout << "(iTag=" << iTag << ")Entries in nEventsPhiPi0 if different orientations: " << phiAMO_pi0->GetEntries() << endl;
					if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVst1_binS23/phiYieldFits"+tagPi0[iTag]+"_MetapBin"+to_string(iMass)+"_tpi0Bin"+to_string(it)+"_iSet"+to_string(iSet)+".png").c_str()); }

				}

				// *****************************
				// Now we overlay the pi0 and eta beam asymmetries 
				// *****************************
				cout << "Starting asymmetry fit for eta" << endl;
				auto leg1 = new TLegend(0.7,0.15,0.9,0.35);
	                        allCanvases = new TCanvas("","",1440,900);
                                TPaveLabel* title = new TPaveLabel(0.1,0.91,0.9,0.97,"Beam Asymmetry 0/90 AND 45/135 Orientation");
                                title->SetTextSize(0.8);
                                title->Draw();
				allCanvases->Divide(2,1,0,0);
				allCanvases->cd(1);
				gPad->SetBottomMargin(0.15);
				auto gr_000 = new TGraphErrors(num_tBins,tBins,asymmetries_000_eta[iTag][iSet][iMass],tBins_err,asymmetries_000_eta_err[iTag][iSet][iMass]);
				leg1->AddEntry(gr_000,"0/90","lep");
				gr_000->SetMarkerColor(4);
				gr_000->SetMarkerStyle(21);
				gr_000->SetLineColor(4);
				gr_000->GetXaxis()->SetTitle("t_{#eta} (GeV^{2})");
				gr_000->GetXaxis()->SetTitleSize(0.05);
				gr_000->Draw("AP");
				gr_000->SetTitle((to_string(minPi0P[iMass])+" < M(#pi^{0}p) < "+to_string(maxPi0P[iMass])).c_str());
				gr_000->GetHistogram()->SetMaximum(1.2);
				gr_000->GetHistogram()->SetMinimum(-1);
				gr_000 = new TGraphErrors(num_tBins,tBins,asymmetries_045_eta[iTag][iSet][iMass],tBins_err,asymmetries_045_eta_err[iTag][iSet][iMass]);
				leg1->AddEntry(gr_000,"45/135","lep");
				gr_000->SetLineColor(2);
				gr_000->SetMarkerColor(2);
				gr_000->SetMarkerStyle(20);
				gr_000->Draw("P SAME");
				leg1->Draw();
		
				allCanvases->cd(2);
				gPad->SetBottomMargin(0.15);
				cout << "Starting asymmetry fit for pi0" << endl;
				auto gr_045 = new TGraphErrors(num_tBins,tBins,asymmetries_000_pi0[iTag][iSet][iMass],tBins_err,asymmetries_000_pi0_err[iTag][iSet][iMass]);
				gr_045->SetMarkerColor(4);
				gr_045->SetLineColor(4);
				gr_045->SetMarkerStyle(21);
				gr_045->GetXaxis()->SetTitleSize(0.06);
				gr_045->GetXaxis()->SetTitle("t_{#pi} (GeV^{2})");
				gr_045->Draw("AP");
				gr_045->SetTitle((to_string(minEtaP[iMass])+" < M(#eta p) < "+to_string(maxEtaP[iMass])).c_str());
				gr_045->GetHistogram()->SetMaximum(1.2);
				gr_045->GetHistogram()->SetMinimum(-1);
				gr_045 = new TGraphErrors(num_tBins,tBins,asymmetries_045_pi0[iTag][iSet][iMass],tBins_err,asymmetries_045_pi0_err[iTag][iSet][iMass]);
				gr_045->SetMarkerColor(2);
				gr_045->SetLineColor(2);
				gr_045->SetMarkerStyle(20);
				gr_045->Draw("P SAME");
				if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVst1_binS23/asymVsMbaryon"+tag[iTag]+"_iteta"+to_string(iMass)+"_iSet"+to_string(iSet)+".png").c_str()); }
			}
		}
	}
}
















