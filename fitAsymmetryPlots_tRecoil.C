// This one is for Asymmetry vs t_recoil which is equal to u3 in vincent/colin language.
#include "/d/grid13/gluex/gluex_top/gluex_style.C"

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



string fileType="C";

void fitAsymmetryPlots_tRecoil(){
        gluex_style();
        // Do some cleaning
        gSystem->Exec("rm -rf asymmetryPlots/SigVsU3");
        gSystem->Exec("mkdir asymmetryPlots/SigVsU3");
        //////////////////

	static const int nSetsBS=1; // (number of sets to bootstrap-1). 1 is for not bootstrapping at all
	int maxPrintBS=1; // Only show up to this value when going through nSetsBS. This includes the "full" set in the count

	// seems like this this would somehow instantiate gMinuit maybe? If I dont do this I get some errors
	//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // If I dont get this I get the error "Error in <TMinuitMinimizer::Scan>:  Error executing command SCAN" 
	static const int num_tBins=5;
	static const int num_trecoilBins=4;
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
	double fluxRatios_45_135[nDataSets] = {  4.076065e+12/4.095013e+12 , 1.02261, 1.03254 };
	string dataSetTag[nDataSets] = { "2017" , "2018_1", "2018_8" };
	string dataFolders[nDataSets] = {"deg000_data_2017" , "deg000_data_2018_1", "deg000_data_2018_8"};

	static const int nTagEta = 1;
	string tagEta[nTagEta] = {""};
	string tagPi0[nTagEta] = {""};
	string tag[nTagEta] = {""}; // same as above but a single name to group them both.
	std::vector<int> nEventsPhiEta[nTagEta]; // second one is for the vanHove selected
	std::vector<int> nEventsPhiPi0[nTagEta];

	// *****************************
	// Loading histograms for  2017, 2018_1, 2018_2, then scaling yield by the flux ratio. Then we can finally add the data sets together
	// *****************************
	//
	//  ************* these are not really used, only as a place holder before scaling them to get the total *************
	
	// There are 10 histograms as shown below
	//static const int nHists=10;
	cout << "Shape of the tensor of histograms: ["<<nTagEta<<"]["<<num_tBins<<"]["<<num_trecoilBins<<"]["<<nDataSets<<"]["<<nSetsBS<<"]"<< endl;
	TH1F *phi000_eta_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS]; // 2 for the tag and 3 for the datasets
	TH1F *phi045_eta_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS]; 
	TH1F *phi090_eta_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS];
	TH1F *phi135_eta_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS];
	TH1F *phiAMO_eta_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS];
	TH1F *phi000_pi0_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS];
	TH1F *phi045_pi0_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS];
	TH1F *phi090_pi0_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS];
	TH1F *phi135_pi0_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS];
	TH1F *phiAMO_pi0_unscaled[nTagEta][num_tBins][num_trecoilBins][nDataSets][nSetsBS];
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
				string dataFileName = "/d/grid15/ln16/pi0eta/092419/degALL_data_"+dataSetTag[iData]+"_BA_hists_DSelector.root";
				//string dataFileName = "/d/grid15/ln16/pi0eta/092419/newGraphs_histValues/rootFiles/deg000_data_"+dataSetTag[iData]+"_hists_DSelector.root";
				TFile *dataFile = new TFile(dataFileName.c_str());
				cout << "LOADING ROOT FILE: " << dataFileName << endl; 
				for (int iTag=0; iTag < nTagEta; ++iTag){
					for (int it=0; it<num_tBins; ++it){
						for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){
							string affixEta = "trecoilBin"+to_string(itrecoilBin)+"_tetaBin"+to_string(it);
							string affixPi0 = "trecoilBin"+to_string(itrecoilBin)+"_tpi0Bin"+to_string(it);
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_045_"+affixEta).c_str(), phi045_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_045_"+affixPi0).c_str(), phi045_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_090_"+affixEta).c_str(), phi090_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_090_"+affixPi0).c_str(), phi090_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_AMO_"+affixEta).c_str(), phiAMO_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_AMO_"+affixPi0).c_str(), phiAMO_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);
							// Need to scale these para yields by flux ratio.
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_000_"+affixEta).c_str(), phi000_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
							cout << "phi000_eta_unscaled["<<iTag<<"]["<<it<<"]["<<itrecoilBin<<"]["<<iData<<"]["<<iSet<<"]" << endl;
							cout << (" -- prodPlanePSphi"+tagEta[iTag]+"_000_"+affixEta).c_str() << endl;
							cout << " -- nentries=" << phi000_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]->GetEntries() << endl;
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_000_"+affixPi0).c_str(), phi000_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_135_"+affixEta).c_str(), phi135_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
							dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_135_"+affixPi0).c_str(), phi135_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);
						}
					}
				}
			}
		}
	}
	//else { 
	//	// *************************************
	//	// Loading the data in from the textFile
	//	// ONLY DO BS FOR LOWEST trecoil THRESHOLD
	//	// *************************************
	//	string baseFileLoc = "/d/grid15/ln16/pi0eta/092419/newGraphs_histValues/";
	//	std::vector<double> values[nTagEta][num_trecoilBins][nDataSets][nHists];
	//	std::vector<double> weights[nTagEta][num_trecoilBins][nDataSets][nHists];
	//	// These two vectors will hold all the BA and their errors respectively. The array of two will be for the two orientations
	//	std::vector<double> sigmas[2];
	//	std::vector<double> sigma_errors[2];
	//	
	//	double value;
	//	double weight;

	//	// This histogram was used to check atleast the shape and counts are correct. We dont do FR scaling yet but w.e.
	//	//TH1F* prodPlane = new TH1F("","",100,-180,180);
	//	for ( int iData=0; iData<nDataSets; ++iData){ // the dataset we use, 2017 or the two 2018 sets
	//		for (int iTag=0; iTag < nTagEta; ++iTag){ // nothing or backwardPi0/Etap
	//			for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){ // t bins
	//				string name0 = "prodPlanePSphi"+tagEta[iTag]+"_045_trecoil_fastEtaBin"+to_string(itrecoilBin);
	//				string name1 = "prodPlanePSphi"+tagPi0[iTag]+"_045_trecoil_fastPi0Bin"+to_string(itrecoilBin);
	//				string name2 = "prodPlanePSphi"+tagEta[iTag]+"_090_trecoil_fastEtaBin"+to_string(itrecoilBin);
	//				string name3 = "prodPlanePSphi"+tagPi0[iTag]+"_090_trecoil_fastPi0Bin"+to_string(itrecoilBin);
	//				string name4 = "prodPlanePSphi"+tagEta[iTag]+"_AMO_trecoil_fastEtaBin"+to_string(itrecoilBin);
	//				string name5 = "prodPlanePSphi"+tagPi0[iTag]+"_AMO_trecoil_fastPi0Bin"+to_string(itrecoilBin);
	//				string name6 = "prodPlanePSphi"+tagEta[iTag]+"_000_trecoil_fastEtaBin"+to_string(itrecoilBin);
	//				string name7 = "prodPlanePSphi"+tagPi0[iTag]+"_000_trecoil_fastPi0Bin"+to_string(itrecoilBin);
	//				string name8 = "prodPlanePSphi"+tagEta[iTag]+"_135_trecoil_fastEtaBin"+to_string(itrecoilBin);
	//				string name9 = "prodPlanePSphi"+tagPi0[iTag]+"_135_trecoil_fastPi0Bin"+to_string(itrecoilBin);
	//				string nameTag = "_iData"+to_string(iData)+"_iSet0"; 

        //				phi045_eta_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name0+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //				phi045_pi0_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name1+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //				phi090_eta_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name2+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //				phi090_pi0_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name3+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //				phiAMO_eta_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name4+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //				phiAMO_pi0_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name5+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //				phi000_eta_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name6+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //				phi000_pi0_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name7+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //				phi135_eta_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name8+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //				phi135_pi0_unscaled[iTag][it][itrecoilBin][iData][0] = new TH1F((name9+nameTag).c_str(), "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
	//				
	//				ifstream* inFile0 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name0+".txt").c_str());
	//				cout << "First file located at: " << baseFileLoc+dataFolders[iData]+"/"+name0 << endl;
	//				ifstream* inFile1 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name1+".txt").c_str());
	//				ifstream* inFile2 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name2+".txt").c_str());
	//				ifstream* inFile3 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name3+".txt").c_str());
	//				ifstream* inFile4 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name4+".txt").c_str());
	//				ifstream* inFile5 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name5+".txt").c_str());
	//				// Have to scale the para yields by flux ratio.
	//				ifstream* inFile6 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name6+".txt").c_str());
	//				ifstream* inFile7 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name7+".txt").c_str());
	//				ifstream* inFile8 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name8+".txt").c_str());
	//				ifstream* inFile9 = new ifstream((baseFileLoc+dataFolders[iData]+"/"+name9+".txt").c_str());
	//				// Arrange everything in a vector so we can easily loop through them
	//				std::vector<ifstream*> inFiles={inFile0,inFile1,inFile2,inFile3,inFile4,inFile5,inFile6,inFile7,inFile8,inFile9};
	//				// since these are pointers I think after I fill the histograms here I could use the oringal histogram pointer names and it should work still
	//				std::vector<TH1F*> hists = {phi045_eta_unscaled[iTag][itrecoilBin][iData][0],phi045_pi0_unscaled[iTag][itrecoilBin][iData][0]
	//					,phi090_eta_unscaled[iTag][itrecoilBin][iData][0]
	//					,phi090_pi0_unscaled[iTag][itrecoilBin][iData][0],phiAMO_eta_unscaled[iTag][itrecoilBin][iData][0],phiAMO_pi0_unscaled[iTag][itrecoilBin][iData][0]
	//					,phi000_eta_unscaled[iTag][itrecoilBin][iData][0],phi000_pi0_unscaled[iTag][itrecoilBin][iData][0],phi135_eta_unscaled[iTag][itrecoilBin][iData][0]
	//					,phi135_pi0_unscaled[iTag][itrecoilBin][iData][0]};
	//				int nevents[nHists]; // since we split up the loop for the full dataset and the bootstrapped samples we have to keep track of how much events each hist has
	//				for ( int iHist=0; iHist<nHists; ++iHist){
	//					while (*(inFiles[iHist]) >> value >> weight){
	//						//cout << value << endl;
	//						values[iTag][itrecoilBin][iData][iHist].push_back(value);
	//						weights[iTag][itrecoilBin][iData][iHist].push_back(weight);
	//						hists[iHist]->Fill(value,weight);
	//					}
	//					hists[iHist]->Print();
	//					nevents[iHist] = values[iTag][itrecoilBin][iData][iHist].size();	
	//					cout << "(iData=" << iData << ")nevents in histogram: " << nevents[iHist] << endl;
	//				}

	//				// We will now bootstrap some datasets and calculate std of sigma
	//				for (int iSet=1; iSet<nSetsBS; ++iSet){ // we add 1 since using the full dataset takes the first element of the array
	//					string nameTag = "_iData"+to_string(iData)+"_iSet"+to_string(iSet);
	//					string cutTag = "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees";
        //					phi045_eta_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name0+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
        //					phi045_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name1+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
        //					phi090_eta_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name2+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
        //					phi090_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name3+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
        //					phiAMO_eta_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name4+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
        //					phiAMO_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name5+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
        //					phi000_eta_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name6+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
        //					phi000_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name7+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
        //					phi135_eta_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name8+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
        //					phi135_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet] = new TH1F((name9+nameTag).c_str(), cutTag.c_str(), 40, -180, 180);
	//					hists = {phi045_eta_unscaled[iTag][it][itrecoilBin][iData][iSet],phi045_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet],phi090_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]
	//						,phi090_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet],phiAMO_eta_unscaled[iTag][it][itrecoilBin][iData][iSet],phiAMO_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]
	//						,phi000_eta_unscaled[iTag][it][itrecoilBin][iData][iSet],phi000_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet],phi135_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]
	//						,phi135_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]};
	//					int randIdx;
	//					for ( int iHist=0; iHist<nHists; ++iHist){
	//						for (int ievent=0; ievent<nevents[iHist]; ++ievent){
	//							randIdx = rand() % nevents[iHist];		
	//							value = values[iTag][itrecoilBin][iData][iHist][randIdx];
	//							weight = weights[iTag][itrecoilBin][iData][iHist][randIdx];
	//							hists[iHist]->Fill(value, weight);
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	
	
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// ************************************************** NOW THAT WE LOADED THE DATA WE CAN SCALE IT AND FIT IT *************************************************************
	// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------------------------------------------------------------------------------------------
	double asymmetries_000_eta[nTagEta][nSetsBS][num_trecoilBins][num_tBins];
	double asymmetries_045_eta[nTagEta][nSetsBS][num_trecoilBins][num_tBins];
	double asymmetries_000_eta_err[nTagEta][nSetsBS][num_trecoilBins][num_tBins];
	double asymmetries_045_eta_err[nTagEta][nSetsBS][num_trecoilBins][num_tBins];
	double asymmetries_000_pi0[nTagEta][nSetsBS][num_trecoilBins][num_tBins];
	double asymmetries_045_pi0[nTagEta][nSetsBS][num_trecoilBins][num_tBins];
	double asymmetries_000_pi0_err[nTagEta][nSetsBS][num_trecoilBins][num_tBins];
	double asymmetries_045_pi0_err[nTagEta][nSetsBS][num_trecoilBins][num_tBins];
	cout << "--------------------------------------------\nLOADING AND SCALING HISTOGRAMS\n--------------------------------------------" << endl;
	cout << "As a simple test we can GetMaximum before and after the scaling to see if it actually worked" << endl;
	cout << "\tWe will follow phi000_eta_total throughout the process" << endl;
	for (int iSet=0; iSet <nSetsBS; ++iSet) {
		// These are the weighted summed histograms. weighted according to the Flux ratio
		TH1F *phi000_eta_total[nTagEta][num_tBins][num_trecoilBins];
		TH1F *phi045_eta_total[nTagEta][num_tBins][num_trecoilBins];
		TH1F *phi090_eta_total[nTagEta][num_tBins][num_trecoilBins];
		TH1F *phi135_eta_total[nTagEta][num_tBins][num_trecoilBins];
		TH1F *phiAMO_eta_total[nTagEta][num_tBins][num_trecoilBins];
		TH1F *phi000_pi0_total[nTagEta][num_tBins][num_trecoilBins];
		TH1F *phi045_pi0_total[nTagEta][num_tBins][num_trecoilBins];
		TH1F *phi090_pi0_total[nTagEta][num_tBins][num_trecoilBins];
		TH1F *phi135_pi0_total[nTagEta][num_tBins][num_trecoilBins];
		TH1F *phiAMO_pi0_total[nTagEta][num_tBins][num_trecoilBins];
		// Now we sum the scaled histograms over the 3 datasets
		double maximum_before;
		double maximum_after;
		int totalEntriesInPhi000_eta_total;
		for (int iTag=0; iTag < nTagEta; ++iTag) {
			for (int it=0; it<num_tBins; ++it){
				for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){
					totalEntriesInPhi000_eta_total=0;
					cout << "--------------\n" << tagEta[0] << "\n--------------" << endl;
					// We first have to clone the first dataTag (2017 run) and scale by flux ratio if in para config
					// we will scale the total in this to track the number of maxima and entires before and after scaling
					phi000_eta_total[iTag][it][itrecoilBin] = (TH1F *)phi000_eta_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();
					phi000_pi0_total[iTag][it][itrecoilBin] = (TH1F *)phi000_pi0_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();
					phi135_eta_total[iTag][it][itrecoilBin] = (TH1F *)phi135_eta_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();
					phi135_pi0_total[iTag][it][itrecoilBin] = (TH1F *)phi135_pi0_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();
					maximum_before = phi000_eta_total[iTag][it][itrecoilBin]->GetMaximum();
					cout << "-Maximum before scaling of phi000_eta_total[" << iTag << "][" << it << "][" << itrecoilBin << "]: " << maximum_before << endl;
					cout << "-NEntries before scaling of phi000_eta_total[" << iTag << "][" << it << "]["<< itrecoilBin << "]: " << phi000_eta_total[iTag][it][itrecoilBin]->GetEntries() << endl;
					phi000_eta_total[iTag][it][itrecoilBin]->Scale(fluxRatios_90_0[0]);
					phi000_pi0_total[iTag][it][itrecoilBin]->Scale(fluxRatios_90_0[0]);
					phi135_eta_total[iTag][it][itrecoilBin]->Scale(fluxRatios_45_135[0]);
					phi135_pi0_total[iTag][it][itrecoilBin]->Scale(fluxRatios_45_135[0]);
					maximum_after = phi000_eta_total[iTag][it][itrecoilBin]->GetMaximum();
					cout << "-Maximum after scaling of phi000_eta_total[" << iTag << "]" << it << "]["<< itrecoilBin << "]: " << maximum_after << endl;
					cout << "-NEntries after scaling of phi000_eta_total[" << iTag << "][" << it << "]["<< itrecoilBin << "]: " << phi000_eta_total[iTag][it][itrecoilBin]->GetEntries() << endl; 
					cout << "-Scale factor vs Expected: " << maximum_after/maximum_before << " vs " << fluxRatios_90_0[0] << endl; 
					phi045_eta_total[iTag][it][itrecoilBin] = (TH1F *)phi045_eta_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();
					phi045_pi0_total[iTag][it][itrecoilBin] = (TH1F *)phi045_pi0_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();
					phi090_eta_total[iTag][it][itrecoilBin] = (TH1F *)phi090_eta_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();
					phi090_pi0_total[iTag][it][itrecoilBin] = (TH1F *)phi090_pi0_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();
					phiAMO_eta_total[iTag][it][itrecoilBin] = (TH1F *)phiAMO_eta_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();
					phiAMO_pi0_total[iTag][it][itrecoilBin] = (TH1F *)phiAMO_pi0_unscaled[iTag][it][itrecoilBin][0][iSet]->Clone();

					totalEntriesInPhi000_eta_total += phi000_eta_total[iTag][it][itrecoilBin]->GetEntries();
					cout << "-Cumulative entries after adding dataSet=0: " << totalEntriesInPhi000_eta_total << endl;
					// for the last 2 runs we simply add with a weight if in para config and add with weight=1 if in perp config
					for (int iData=1; iData < nDataSets; ++iData){
						// keeping a count of the total entries to check if we are doing the Add right
						totalEntriesInPhi000_eta_total += phi000_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]->GetEntries();
						cout << "-Cumulative entries after adding dataSet=" << iData << ": " << totalEntriesInPhi000_eta_total << endl;

						// scale the rest of the data sets before adding them to the scaled+cloned 2017 dataset
						// we will never use unscaled histograms anymore so we can just scale it directly
						phi000_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]->Scale(fluxRatios_90_0[iData]);
						phi000_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]->Scale(fluxRatios_90_0[iData]);
						phi135_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]->Scale(fluxRatios_45_135[iData]);
						phi135_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]->Scale(fluxRatios_45_135[iData]);
						phi000_eta_total[iTag][it][itrecoilBin]->Add(phi000_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
        	        			phi000_pi0_total[iTag][it][itrecoilBin]->Add(phi000_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);
        	        			phi135_eta_total[iTag][it][itrecoilBin]->Add(phi135_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
        	        			phi135_pi0_total[iTag][it][itrecoilBin]->Add(phi135_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);

						// Dont have to scale these guys so just add them
        	        			phi045_eta_total[iTag][it][itrecoilBin]->Add(phi045_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
        	        			phi045_pi0_total[iTag][it][itrecoilBin]->Add(phi045_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);
        	        			phi090_eta_total[iTag][it][itrecoilBin]->Add(phi090_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
        	        			phi090_pi0_total[iTag][it][itrecoilBin]->Add(phi090_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);
        	        			phiAMO_eta_total[iTag][it][itrecoilBin]->Add(phiAMO_eta_unscaled[iTag][it][itrecoilBin][iData][iSet]);
        	        			phiAMO_pi0_total[iTag][it][itrecoilBin]->Add(phiAMO_pi0_unscaled[iTag][it][itrecoilBin][iData][iSet]);
					}
				}
			}
		}	

		cout << "--------------------------------------------\nSTARTING FITTING\n--------------------------------------------" << endl;


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

		//double mintrecoil = 1.6;
		//double maxtrecoil = 2.8;
		//double mBinSize = (maxtrecoil-mintrecoil)/5;

		Int_t fitStatus;
		for (int iTag=0; iTag<nTagEta; ++iTag){
			for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){ //num_trecoilBins
			        for (int it=0; it<num_tBins; ++it){
				        double lowerT = it*0.2;
				        double upperT = (it+1)*0.2;
                                        double midT = (it+0.5)*0.2; 
					// Need to define the polarization with the systematic shift
					// the shift for 0/90 is 3.1 and shifting 45/135 is 3.2. Since 0 and 135/-45 is para the beginning orientation is 0 and -45 whih then we shift by 3.1, 3.2 respectively
					string fitOption = "E S";
					TH1F *phi000_eta = phi000_eta_total[iTag][it][itrecoilBin];
					TH1F *phi045_eta = phi045_eta_total[iTag][it][itrecoilBin];
					TH1F *phi090_eta = phi090_eta_total[iTag][it][itrecoilBin];
					TH1F *phi135_eta = phi135_eta_total[iTag][it][itrecoilBin];
					TH1F *phiAMO_eta = phiAMO_eta_total[iTag][it][itrecoilBin];
					TH1F *phi000_pi0 = phi000_pi0_total[iTag][it][itrecoilBin];
					TH1F *phi045_pi0 = phi045_pi0_total[iTag][it][itrecoilBin];
					TH1F *phi090_pi0 = phi090_pi0_total[iTag][it][itrecoilBin];
					TH1F *phi135_pi0 = phi135_pi0_total[iTag][it][itrecoilBin];
					TH1F *phiAMO_pi0 = phiAMO_pi0_total[iTag][it][itrecoilBin];
					TH1F *phis_eta[4] =  {phi000_eta, phi045_eta, phi090_eta, phi135_eta};
					TH1F *phis_pi0[4] =  {phi000_pi0, phi045_pi0, phi090_pi0, phi135_pi0};

					// *****************************
					// Defining the trecoilBins which will be shared
					// *****************************
					//trecoilBins[itrecoilBin] = mintrecoil+(itrecoilBin+0.5)*mBinSize;
					//trecoilBins_err[itrecoilBin] = mBinSize/2;
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
					fit_asym->FixParameter(3,Phi0_0_90);

					TFitResultPtr fitPointer = asymmetry000_090_eta->Fit(fit_asym,fitOption.c_str());
					asymmetries_000_eta[iTag][iSet][itrecoilBin][it] = fit_asym->GetParameter(2);
					asymmetries_000_eta_err[iTag][iSet][itrecoilBin][it] = fit_asym->GetParError(2);
					asymmetry000_090_eta->Draw("SAME");
					asymmetry000_090_eta->SetAxisRange(-0.5,0.6,"Y");
					TGraph* likelihoodFit_000_090_eta = new TGraph(500); 
					//fitPointer->Scan(2,likelihoodFit_000_090_eta,0,1);

					allCanvases->cd(2);
					fit_asym->SetParameters(0.35,0.35,0.5,Phi0_45_135);
					fit_asym->FixParameter(0,0.35);
					fit_asym->FixParameter(1,0.35);
					fit_asym->FixParameter(3,Phi0_45_135);

					fitPointer = asymmetry045_135_eta->Fit(fit_asym,fitOption.c_str());
					asymmetries_045_eta[iTag][iSet][itrecoilBin][it] = fit_asym->GetParameter(2);
					asymmetries_045_eta_err[iTag][iSet][itrecoilBin][it] = fit_asym->GetParError(2);
					asymmetry045_135_eta->Draw("SAME");
					asymmetry045_135_eta->SetAxisRange(-0.5,0.6,"Y");
					TGraph* likelihoodFit_045_135_eta = new TGraph(500); 
					//fitPointer->Scan(2,likelihoodFit_045_135_eta,0,1);

					if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVsU3/asymmetry"+tagEta[iTag]+"_trecoilBin"+to_string(itrecoilBin)+"_tetaBin"+to_string(it)+"_iSet"+to_string(iSet)+"."+fileType).c_str()); }
					
					//allCanvases->Clear();
					//allCanvases->Divide(2,1);
					//allCanvases->cd(1);
					//likelihoodFit_000_090_eta->Draw("ALP");
					//likelihoodFit_000_090_eta->SetTitle(("0/90 - tBin"+to_string(itrecoilBin)).c_str());
					//allCanvases->cd(2);
					//likelihoodFit_045_135_eta->Draw("ALP");
					//likelihoodFit_045_135_eta->SetTitle(("045_135 - tBin"+to_string(itrecoilBin)).c_str());
					//if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVsU3/likelihood"+tagEta[iTag]+"_trecoilBin"+to_string(iMass)+"_trecoil_fastEtaBin"+to_string(itrecoilBin)+"_iSet"+to_string(iSet)+".root").c_str()); }

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
					fit_asym->FixParameter(3,Phi0_0_90);
					fitPointer = asymmetry000_090_pi0->Fit(fit_asym,fitOption.c_str());
					asymmetries_000_pi0[iTag][iSet][itrecoilBin][it] = fit_asym->GetParameter(2);
					asymmetries_000_pi0_err[iTag][iSet][itrecoilBin][it] = fit_asym->GetParError(2);
					asymmetry000_090_pi0->Draw("SAME");
					asymmetry000_090_pi0->SetAxisRange(-0.5,0.6,"Y");
					TGraph* likelihoodFit_000_090_pi0 = new TGraph(50); 
					//fitPointer->Scan(2,likelihoodFit_000_090_pi0,-2,2);
		
					allCanvases->cd(2);
					fit_asym->SetParameters(0.35,0.35,0.5,Phi0_45_135);
					fit_asym->FixParameter(0,0.35);
					fit_asym->FixParameter(1,0.35);
					fit_asym->FixParameter(3,Phi0_45_135);
					fitPointer = asymmetry045_135_pi0->Fit(fit_asym,fitOption.c_str());
					asymmetry045_135_pi0->Draw("SAME");
					asymmetry045_135_pi0->SetAxisRange(-0.5,0.6,"Y");
					asymmetries_045_pi0[iTag][iSet][itrecoilBin][it] = fit_asym->GetParameter(2);
					asymmetries_045_pi0_err[iTag][iSet][itrecoilBin][it] = fit_asym->GetParError(2);
					TGraph* likelihoodFit_045_135_pi0 = new TGraph(50); 
					//fitPointer->Scan(2,likelihoodFit_045_135_pi0,-2,2);

					if ( iSet < maxPrintBS ){ 
						allCanvases->SaveAs(("asymmetryPlots/SigVsU3/asymmetry"+tagPi0[iTag]+"_trecoilBin"+to_string(itrecoilBin)+"_tpi0Bin"+to_string(it)+"_iSet"+to_string(iSet)+"."+fileType).c_str()); 
					}

					//allCanvases->Clear();
					//allCanvases->Divide(2,1);
					//allCanvases->cd(1);
					//likelihoodFit_000_090_eta->Draw("ALP");
					//likelihoodFit_000_090_eta->SetTitle(("0/90 - tBin"+to_string(itrecoilBin)).c_str());
					//allCanvases->cd(2);
					//likelihoodFit_045_135_eta->Draw("ALP");
					//likelihoodFit_045_135_eta->SetTitle(("045/135 - tBin"+to_string(itrecoilBin)).c_str());
					//if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVsU3/likelihood"+tagEta[iTag]+"_trecoilBin"+to_string(iMass)+"_trecoil_fastEtaBin"+to_string(itrecoilBin)+"_iSet"+to_string(iSet)+".root").c_str()); }
		
					// *****************************
					// Fitting prodPlanePhi for eta
					// *****************************
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
					cout << "Doing flat fit to AMO for eta" << endl;
					fitStatus = phiAMO_eta->Fit(fit_flat,"RQE");
					phiAMO_eta->Draw("SAME");
					if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVsU3/phiYieldFits"+tagEta[iTag]+"_trecoilBin"+to_string(itrecoilBin)+"_tetaBin"+to_string(it)+"_iSet"+to_string(iSet)+"."+fileType).c_str()); }
					
					// *****************************
					// Fitting prodPlanePhi for pi0
					// *****************************
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
					cout << "Doing flat fit to AMO for pi0" << endl;
					fitStatus = phiAMO_pi0->Fit(fit_flat,"E S");
					phiAMO_pi0->Draw("SAME");
					cout << "(iTag=" << iTag << ")Entries in nEventsPhiPi0 if different orientations: " << phiAMO_pi0->GetEntries() << endl;
					if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVsU3/phiYieldFits"+tagPi0[iTag]+"_trecoilBin"+to_string(itrecoilBin)+"_tpi0Bin"+to_string(it)+"_iSet"+to_string(iSet)+"."+fileType).c_str()); }

				}

				// *****************************
				// Now we overlay the pi0 and eta beam asymmetries 
				// *****************************
				auto leg1 = new TLegend(0.7,0.15,0.9,0.35);
	                        allCanvases = new TCanvas("","",1440,900);
                                TPaveLabel* title;
                                if (itrecoilBin==0){ 
                                    title = new TPaveLabel(0.1,0.91,0.9,0.96,"0 < t_{proton} < 0.5");
                                }
                                else if (itrecoilBin==1){ 
                                    title = new TPaveLabel(0.1,0.91,0.9,0.96,"0.5 < t_{proton} < 1");
                                }
                                else if (itrecoilBin==2){ 
                                    title = new TPaveLabel(0.1,0.91,0.9,0.96,"1 < t_{proton}");
                                }
                                else{ 
                                    title = new TPaveLabel(0.1,0.91,0.9,0.96,"all t_{proton}");
                                }
                                title->SetTextSize(0.8);
                                title->Draw();
				allCanvases->Divide(2,1,0,0);
				allCanvases->cd(1);
				gPad->SetBottomMargin(0.15);
				auto gr_000 = new TGraphErrors(num_tBins,tBins,asymmetries_000_eta[iTag][iSet][itrecoilBin],tBins_err,asymmetries_000_eta_err[iTag][iSet][itrecoilBin]);
				leg1->AddEntry(gr_000,"t_{#eta}","lep");
				gr_000->SetTitle("0/90 Orientation");
				gr_000->SetMarkerColor(4);
				gr_000->SetMarkerStyle(21);
				gr_000->SetLineColor(4);
				gr_000->Draw("AP");
				gr_000->GetXaxis()->SetTitle("t_{#eta}/t_{#pi})");
                                gr_000->GetYaxis()->SetTitle("#Sigma");
				//gr_000->GetXaxis()->SetTitleSize(0.06);
                                gr_000->GetXaxis()->SetTitleSize(0.08);
                                gr_000->GetYaxis()->SetTitleSize(0.08);
				gr_000->GetHistogram()->SetMaximum(1.2);
				gr_000->GetHistogram()->SetMinimum(-1);
				gr_000 = new TGraphErrors(num_tBins,tBins,asymmetries_000_pi0[iTag][iSet][itrecoilBin],tBins_err,asymmetries_000_pi0_err[iTag][iSet][itrecoilBin]);
				leg1->AddEntry(gr_000,"t_{#pi^{0}}","lep");
				gr_000->SetTitle("0/90 Orientation");
				gr_000->SetLineColor(2);
				gr_000->SetMarkerColor(2);
				gr_000->SetMarkerStyle(20);
				gr_000->Draw("P SAME");
				leg1->Draw();
		
				allCanvases->cd(2);
				gPad->SetBottomMargin(0.15);
				auto gr_045 = new TGraphErrors(num_tBins,tBins,asymmetries_045_eta[iTag][iSet][itrecoilBin],tBins_err,asymmetries_045_eta_err[iTag][iSet][itrecoilBin]);
				gr_045->SetTitle("45/135 Orientation");
				gr_045->SetMarkerColor(4);
				gr_045->SetLineColor(4);
				gr_045->SetMarkerStyle(21);
				gr_045->GetXaxis()->SetTitle("t_{#eta}/t_{#pi}");
                                gr_045->GetYaxis()->SetTitle("#Sigma");
				//gr_000->GetXaxis()->SetTitleSize(0.06);
				gr_045->Draw("AP");
                                gr_045->GetXaxis()->SetTitleSize(0.08);
                                gr_045->GetYaxis()->SetTitleSize(0.08);
				gr_045->GetHistogram()->SetMaximum(1.2);
				gr_045->GetHistogram()->SetMinimum(-1);
				gr_045 = new TGraphErrors(num_tBins,tBins,asymmetries_045_pi0[iTag][iSet][itrecoilBin],tBins_err,asymmetries_045_pi0_err[iTag][iSet][itrecoilBin]);
				gr_045->SetTitle("45/135 Orientation");
				gr_045->SetMarkerColor(2);
				gr_045->SetLineColor(2);
				gr_045->SetMarkerStyle(20);
				gr_045->Draw("P SAME");
				if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVsU3/asymVst_trecoilBin"+to_string(itrecoilBin)+"_iSet"+to_string(iSet)+"."+fileType).c_str()); }

	                        allCanvases = new TCanvas("","",1440,900);
                                int nBins = sizeof(asymmetries_000_eta[iTag][iSet][itrecoilBin])/sizeof(double);
                                vector<double> asymmetries_000_045_eta;
                                vector<double> asymmetries_000_045_pi0;
                                vector<double> asymmetries_000_045_eta_err;
                                vector<double> asymmetries_000_045_pi0_err;
                                for (int iBin=0; iBin<nBins; ++iBin){
                                   asymmetries_000_045_eta.push_back(0.5*(asymmetries_000_eta[iTag][iSet][itrecoilBin][iBin]+asymmetries_045_eta[iTag][iSet][itrecoilBin][iBin])); 
                                   asymmetries_000_045_pi0.push_back(0.5*(asymmetries_000_pi0[iTag][iSet][itrecoilBin][iBin]+asymmetries_045_pi0[iTag][iSet][itrecoilBin][iBin])); 
                                   double asym_eta_err = 0.5*TMath::Sqrt(asymmetries_000_eta_err[iTag][iSet][itrecoilBin][iBin]*asymmetries_000_eta_err[iTag][iSet][itrecoilBin][iBin]+asymmetries_045_eta_err[iTag][iSet][itrecoilBin][iBin]*asymmetries_045_eta_err[iTag][iSet][itrecoilBin][iBin]);
                                   double asym_pi0_err = 0.5*TMath::Sqrt(asymmetries_000_pi0_err[iTag][iSet][itrecoilBin][iBin]*asymmetries_000_pi0_err[iTag][iSet][itrecoilBin][iBin]+asymmetries_045_pi0_err[iTag][iSet][itrecoilBin][iBin]*asymmetries_045_pi0_err[iTag][iSet][itrecoilBin][iBin]);
                                   asymmetries_000_045_eta_err.push_back(asym_eta_err);
                                   asymmetries_000_045_pi0_err.push_back(asym_pi0_err);
                                }
                                auto gr_avg_eta = new TGraphErrors(num_tBins,tBins,&asymmetries_000_045_eta[0],tBins_err,&asymmetries_000_045_eta_err[0]);
                                auto gr_avg_pi0 = new TGraphErrors(num_tBins,tBins,&asymmetries_000_045_pi0[0],tBins_err,&asymmetries_000_045_pi0_err[0]);
				gr_avg_eta->SetTitle("Beam Asymmetry");
				gr_avg_eta->SetMarkerColor(4);
				gr_avg_eta->SetLineColor(4);
				gr_avg_eta->SetMarkerStyle(21);
				gr_avg_eta->GetXaxis()->SetTitle("t_{#eta}/t_{#pi}");
                                gr_avg_eta->GetYaxis()->SetTitle("#Sigma");
				//gr_000->GetXaxis()->SetTitleSize(0.06);
				gr_avg_eta->Draw("AP");
                                gr_avg_eta->GetXaxis()->SetTitleSize(0.08);
                                gr_avg_eta->GetYaxis()->SetTitleSize(0.08);
				gr_avg_eta->GetHistogram()->SetMaximum(1.2);
				gr_avg_eta->GetHistogram()->SetMinimum(-1);
				gr_avg_pi0->SetMarkerColor(2);
				gr_avg_pi0->SetLineColor(2);
				gr_avg_pi0->SetMarkerStyle(20);
				gr_avg_pi0->Draw("P SAME");
				if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVsU3/avg_asymVst_trecoilBin"+to_string(itrecoilBin)+"_iSet"+to_string(iSet)+"."+fileType).c_str()); }



				allCanvases->Clear();
				allCanvases->Divide(2,1);
				int counter = 0;
   				auto legend_eta = new TLegend(0.6,0.7,0.9,0.9);
   				legend_eta->SetHeader("","C"); // option "C" allows to center the header
   				auto legend_pi0 = new TLegend(0.6,0.7,0.9,0.9);
   				legend_pi0->SetHeader("","C"); // option "C" allows to center the header
				TGraphErrors* gr_etas[numPolarizations];
				TGraphErrors* gr_pi0s[numPolarizations];
				for (int iphi=0; iphi<numPolarizations; ++iphi){
					gr_etas[counter] = new TGraphErrors(num_tBins,tBins,&(phis_eta_PSig[counter][itrecoilBin]),tBins_err,&(phis_eta_PSig_err[counter][itrecoilBin]));
					gr_etas[counter]->SetMarkerStyle(20+counter);
					gr_etas[counter]->SetLineColor(kBlue+counter);
					gr_etas[counter]->SetMarkerColor(kBlue+counter);
					gr_pi0s[counter] = new TGraphErrors(num_tBins,tBins,&(phis_pi0_PSig[counter][itrecoilBin]),tBins_err,&(phis_pi0_PSig_err[counter][itrecoilBin]));
					gr_pi0s[counter]->SetMarkerStyle(20+counter);
					gr_pi0s[counter]->SetLineColor(kRed+counter);
					gr_pi0s[counter]->SetMarkerColor(kRed+counter);

					gr_etas[counter]->SetName(("eta_trecoilBin"+to_string(itrecoilBin)+"_"+names[counter]).c_str());
					gr_pi0s[counter]->SetName(("pi0_trecoilBin"+to_string(itrecoilBin)+"_"+names[counter]).c_str());
					if (counter==0){
						allCanvases->cd(1);
                                                if (itrecoilBin==0){
						    gr_etas[counter]->SetTitle("PSigma t_recoil < 0.5");
						    gr_pi0s[counter]->SetTitle("PSigma t_recoil < 0.5");
                                                }
                                                else if (itrecoilBin==1){
						    gr_etas[counter]->SetTitle("PSigma 0.5 < t_recoil < 1");
						    gr_pi0s[counter]->SetTitle("PSigma 0.5 < t_recoil < 1");
                                                }
                                                else { 
						    gr_etas[counter]->SetTitle("PSigma t_recoil > 1");
						    gr_pi0s[counter]->SetTitle("PSigma t_recoil > 1");
                                                }


						gr_etas[counter]->Draw("AP");
						// apparently we have to AddEntry Before we cd? Or after we draw?
   						legend_eta->AddEntry(("eta_"+names[counter]).c_str(),(names[counter]).c_str(),"lep");
						allCanvases->cd(2);
						gr_pi0s[counter]->Draw("AP");
   						legend_pi0->AddEntry(("pi0_"+names[counter]).c_str(),(names[counter]).c_str(),"lep");
					}
					else {
						allCanvases->cd(1);
						gr_etas[counter]->Draw("P SAME");
						// apparently we have to AddEntry Before we cd? Or after we draw?
   						legend_eta->AddEntry(("eta_"+names[counter]).c_str(),(names[counter]).c_str(),"lep");
						allCanvases->cd(2);
						gr_pi0s[counter]->Draw("P SAME");
   						legend_pi0->AddEntry(("pi0_"+names[counter]).c_str(),(names[counter]).c_str(),"lep");
					}
					gr_etas[counter]->GetHistogram()->SetMaximum(0.3);
					gr_etas[counter]->GetHistogram()->SetMinimum(-0.3);
					gr_pi0s[counter]->GetHistogram()->SetMaximum(0.3);
					gr_pi0s[counter]->GetHistogram()->SetMinimum(-0.3);
					++counter;
				}
				allCanvases->cd(1);
   				legend_eta->Draw();
				allCanvases->cd(2);
   				legend_pi0->Draw();
				if ( iSet < maxPrintBS ){ allCanvases->SaveAs(("asymmetryPlots/SigVsU3/PSigma_trecoilBin"+to_string(itrecoilBin)+"_iSet"+to_string(iSet)+"."+fileType).c_str()); }
			}
		}
		
	}

	// **************************************************
	// Calculate BS errors normalized to the "full" set errors.
	// **************************************************
	// We will keep the same structure but now for each iSet we will calculate the std up to that point. So we can see how the std curve looks like
	// we will also neglect the 0th nSetsBS since it is the "full"
	//if (nSetsBS != 1) {
	//	cout << "\n\n----------------------------------Begin calculating boostrapped errors\n-----------------------------------"<<endl;
	//	static const int nSetsJustBS=nSetsBS-1; // just used to initialize arrays and as a size count for TGraphs.
	//	double asymmetries_000_eta_err_BS[nTagEta][num_trecoilBins][nSetsJustBS]; // nSetsJustBS is the index with values = the running population std
	//	double asymmetries_045_eta_err_BS[nTagEta][num_trecoilBins][nSetsJustBS];
	//	double asymmetries_000_pi0_err_BS[nTagEta][num_trecoilBins][nSetsJustBS];
	//	double asymmetries_045_pi0_err_BS[nTagEta][num_trecoilBins][nSetsJustBS];
	//	double nBoostdstrapped[nSetsJustBS];
	//	double nBoostdstrappedErrs[nSetsJustBS];

	//	for (int iTag=0; iTag < 1; ++iTag) { //nTagEta
	//		// all the std should be zero if there is only 1 element...
	//		for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){
	//			asymmetries_000_eta_err_BS[iTag][itrecoilBin][0] = 0;
	//			asymmetries_045_eta_err_BS[iTag][itrecoilBin][0] = 0;
	//			asymmetries_000_pi0_err_BS[iTag][itrecoilBin][0] = 0;
	//			asymmetries_045_pi0_err_BS[iTag][itrecoilBin][0] = 0;
	//			nBoostdstrapped[0] = 1;
	//			nBoostdstrappedErrs[0] = 0;
	//		}
	//		// we dont begin with iSetJustBS=1 since we would have no elements for the inner loop: for ( int iSet=1; iSet < iSetJustBS; ++iSet )
	//		// dont begin with 2 since we have already filled up that one in the code above since we know it is zero
	//		// so we begin with 3 and end with nSetsBS+2. So if like nSetsBS=6 then: for 3 to 8 gives 3,4,5,6,7 which contains 5 iterations or 5 histograms as expected
	//		ofstream allAsymsBSFile;
	//		allAsymsBSFile.open(("asymmetryPlots/SigVsU3/allAsymsBSFile_iTag"+to_string(iTag)+".txt").c_str(),std::ios_base::trunc);
	//		allAsymsBSFile << "asymmetries_000_eta[iTag][iSet][itrecoilBin]" << endl;
	//		for(int iSetJustBS=3; iSetJustBS<nSetsBS+2; ++iSetJustBS){
	//			nBoostdstrapped[iSetJustBS-2] = iSetJustBS-1;  // so following the logic of nBoostdstrapped[0] = 1 we have to shift iSetJustBS by 2 in the index and by 1 in the rvalue 
	//			nBoostdstrappedErrs[iSetJustBS-2] = 0;
	//			for(int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){
	//				double mean_000_eta = 0;
	//				double mean_045_eta = 0;
	//				double mean_000_pi0 = 0;
	//				double mean_045_pi0 = 0;
	//				for ( int iSet=1; iSet < iSetJustBS; ++iSet ){
	//					mean_000_eta += asymmetries_000_eta[iTag][iSet][itrecoilBin]; // these asymmetries begin with the "full" fits. Only do BS for lowest trecoil threshold
	//					mean_045_eta += asymmetries_045_eta[iTag][iSet][itrecoilBin];
	//					mean_000_pi0 += asymmetries_000_pi0[iTag][iSet][itrecoilBin];
	//					mean_045_pi0 += asymmetries_045_pi0[iTag][iSet][itrecoilBin];
	//				}
	//				mean_000_eta /= (iSetJustBS-1);// need to subtract 1 here since if iSetJustBS=3 then the above for loop would only have 2 elements 
	//				mean_045_eta /= (iSetJustBS-1);
	//				mean_000_pi0 /= (iSetJustBS-1);
        //				mean_045_pi0 /= (iSetJustBS-1); 
	//				double summedSqDiff_000_eta = 0;
	//				double summedSqDiff_045_eta = 0;
	//				double summedSqDiff_000_pi0 = 0;
	//				double summedSqDiff_045_pi0 = 0;
	//				for ( int iSet=1; iSet<iSetJustBS; ++iSet ){ //again we start at 1 since "full" uses the 0 position
	//					summedSqDiff_000_eta += TMath::Power(asymmetries_000_eta[iTag][iSet][itrecoilBin]-mean_000_eta,2);
	//					summedSqDiff_045_eta += TMath::Power(asymmetries_045_eta[iTag][iSet][itrecoilBin]-mean_045_eta,2);
	//					summedSqDiff_000_pi0 += TMath::Power(asymmetries_000_pi0[iTag][iSet][itrecoilBin]-mean_000_pi0,2);
	//					summedSqDiff_045_pi0 += TMath::Power(asymmetries_045_pi0[iTag][iSet][itrecoilBin]-mean_045_pi0,2);
	//				}
	//				summedSqDiff_000_eta /= (iSetJustBS-2); // subtract by 2 here since 1 is from the "full" shift and one is for hte population std
	//				summedSqDiff_000_eta = sqrt(summedSqDiff_000_eta);
	//				// normalized to the full fit's error, the second zero since we do not bootstrap for diff trecoil thresholds
	//				summedSqDiff_000_eta /= asymmetries_000_eta_err[iTag][0][itrecoilBin];
	//				summedSqDiff_045_eta /= (iSetJustBS-2);
	//				summedSqDiff_045_eta = sqrt(summedSqDiff_045_eta);
	//				summedSqDiff_045_eta /= asymmetries_045_eta_err[iTag][0][itrecoilBin]; 
	//				summedSqDiff_000_pi0 /= (iSetJustBS-2);
	//				summedSqDiff_000_pi0 = sqrt(summedSqDiff_000_pi0);
	//				summedSqDiff_000_pi0 /= asymmetries_000_pi0_err[iTag][0][itrecoilBin]; 
	//				summedSqDiff_045_pi0 /= (iSetJustBS-2);
	//				summedSqDiff_045_pi0 = sqrt(summedSqDiff_045_pi0);
	//				summedSqDiff_045_pi0 /= asymmetries_045_pi0_err[iTag][0][itrecoilBin]; 

	//				asymmetries_000_eta_err_BS[iTag][itrecoilBin][iSetJustBS-2] = summedSqDiff_000_eta;
	//				asymmetries_045_eta_err_BS[iTag][itrecoilBin][iSetJustBS-2] = summedSqDiff_045_eta;
	//				asymmetries_000_pi0_err_BS[iTag][itrecoilBin][iSetJustBS-2] = summedSqDiff_000_pi0;
	//				asymmetries_045_pi0_err_BS[iTag][itrecoilBin][iSetJustBS-2] = summedSqDiff_045_pi0;
	//				cout << "(itrecoilBin, iSetJustBS)=(" << itrecoilBin << ", " << iSetJustBS << ") currentSTD:" << summedSqDiff_000_eta << endl;
	//				allAsymsBSFile << "asymmetries_000_eta[" << iTag << "]["<< iSetJustBS-2 << "][" << itrecoilBin << "] = " << asymmetries_000_eta[iTag][iSetJustBS-2][itrecoilBin] << endl;
	//			}
	//		}
	//		allAsymsBSFile.close();

	//		allCanvases->Clear();
	//		allCanvases->Divide(2,2);
	//		allCanvases->cd(1);
	//		TGraph* gr_stds;
	//		for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){
	//			gr_stds = new TGraph(nSetsJustBS, nBoostdstrapped, asymmetries_000_eta_err_BS[iTag][itrecoilBin]);
	//			gr_stds->SetTitle("BS errors 0/90 eta");
	//			//gr_stds->SetMarkerColor(itrecoilBin+1);
	//			//gr_stds->SetMarkerStyle(21);
	//			gr_stds->SetLineColor(itrecoilBin+1);
	//			gr_stds->SetLineWidth(2);
	//			if(itrecoilBin==0){
	//				gr_stds->Draw("AL");
	//				gr_stds->GetHistogram()->SetMaximum(2);
	//				gr_stds->GetHistogram()->SetMinimum(0.0);
	//			}
	//			else{ 
	//				gr_stds->Draw("SAME L");
	//			}
	//		}
	//		allCanvases->cd(2);
	//		for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){
	//			gr_stds = new TGraph(nSetsJustBS, nBoostdstrapped, asymmetries_045_eta_err_BS[iTag][itrecoilBin]);
	//			gr_stds->SetTitle("BS errors 45/135 eta");
	//			//gr_stds->SetMarkerColor(itrecoilBin+1);
	//			//gr_stds->SetMarkerStyle(21);
	//			gr_stds->SetLineColor(itrecoilBin+1);
	//			gr_stds->SetLineWidth(2);
	//			if(itrecoilBin==0){
	//				gr_stds->Draw("AL");
	//				gr_stds->GetHistogram()->SetMaximum(2);
	//				gr_stds->GetHistogram()->SetMinimum(0.0);
	//			}
	//			else{ 
	//				gr_stds->Draw("SAME L");
	//			}
	//		}
	//		allCanvases->cd(3);
	//		for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){
	//			gr_stds = new TGraph(nSetsJustBS, nBoostdstrapped, asymmetries_000_pi0_err_BS[iTag][itrecoilBin]);
	//			gr_stds->SetTitle("BS errors 0/90 pi0");
	//			//gr_stds->SetMarkerColor(itrecoilBin+1);
	//			//gr_stds->SetMarkerStyle(21);
	//			gr_stds->SetLineColor(itrecoilBin+1);
	//			gr_stds->SetLineWidth(2);
	//			if(itrecoilBin==0){
	//				gr_stds->Draw("AL");
	//				gr_stds->GetHistogram()->SetMaximum(2);
	//				gr_stds->GetHistogram()->SetMinimum(0.0);
	//			}
	//			else{ 
	//				gr_stds->Draw("SAME L");
	//			}
	//		}
	//		allCanvases->cd(4);
	//		for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){
	//			gr_stds = new TGraph(nSetsJustBS, nBoostdstrapped, asymmetries_045_pi0_err_BS[iTag][itrecoilBin]);
	//			gr_stds->SetTitle("BS errors 45/135 pi0");
	//			//gr_stds->SetMarkerColor(itrecoilBin+1);
	//			//gr_stds->SetMarkerStyle(21);
	//			gr_stds->SetLineColor(itrecoilBin+1);
	//			gr_stds->SetLineWidth(2);
	//			if(itrecoilBin==0){
	//				gr_stds->Draw("AL");
	//				gr_stds->GetHistogram()->SetMaximum(2);
	//				gr_stds->GetHistogram()->SetMinimum(0.0);
	//			}
	//			else{ 
	//				gr_stds->Draw("SAME L");
	//			}
	//		}
	//		allCanvases->SaveAs(("asymmetryPlots/SigVsU3/bootstrappedErrors_iTag"+to_string(iTag)+"."+fileType).c_str());
	//	}
	//}

	cout << "\n" << endl;
	cout << "Number entries in nEventsPhiEta: " << nEventsPhiEta[0].size() << endl;
	cout << "Number entries in nEventsPhiPi0: " << nEventsPhiPi0[0].size() << endl;
	int currentIndex;
	int sizeOfOrientationList = sizeof(phis_orientation)/sizeof(phis_orientation[0]);
	for (int itrecoilBin=0; itrecoilBin<num_trecoilBins; ++itrecoilBin){
		cout << "--------------------------------------" << endl;
		for (int phiIndex=0; phiIndex<sizeOfOrientationList; ++phiIndex) { 
			currentIndex = itrecoilBin*sizeOfOrientationList + phiIndex;
			cout << "CURRENT INDEX: " << currentIndex << endl;
			cout << "(" << phis_orientation[phiIndex] << ")" << "Number events in teta Bin " << itrecoilBin << ": " << (double)nEventsPhiEta[0][currentIndex] << endl;
			cout << "(" << phis_orientation[phiIndex] << ")" << "Number events in tpi0 Bin " << itrecoilBin << ": " << (double)nEventsPhiPi0[0][currentIndex] << endl;
		}
	}
	cout << "\n\n" << endl;
}
















