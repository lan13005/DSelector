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

void asymmetryMC(){
	int maxPrintBS=1; // Only show up to this value when going through nSetsBS. This includes the "full" set in the count

	// seems like this this would somehow instantiate gMinuit maybe? If I dont do this I get some errors
	//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // If I dont get this I get the error "Error in <TMinuitMinimizer::Scan>:  Error executing command SCAN" 
	static const int num_tBins=5;
	static const int num_mBins=1;
	string phis_orientation[5] =  {"phi000", "phi045", "phi090", "phi135", "phiAMO"};

	//gStyle->SetOptFit();
	//gStyle->SetStatY(1);
	//gStyle->SetStatX(1);
	//gStyle->SetStatW(0.16);
	//gStyle->SetStatH(0.16);
	gStyle->SetOptStat(0);

	TCanvas *allCanvases = new TCanvas("","",1440,900);


	// *****************************
	// Define flux ratios for 2017, 2018_1, 2018_2 
	// *****************************
	static const int nDataSets = 3;
	double fluxRatios_90_0[nDataSets] = {  4.346818e+12/4.188001e+12, 0.965429, 0.918503 };
	double fluxRatios_45_135[nDataSets] = {  4.076065e+12/4.095013e+12, 1.02261, 1.03254 };
	string dataSetTag[nDataSets] = { "data_2017", "data_2018_1", "data_2018_8" };
	string dataFolders[nDataSets] = {"deg000_data_2017", "deg000_data_2018_1", "deg000_data_2018_8"};

	TH1F *phi000_eta_unscaled[num_mBins][num_tBins][nDataSets]; // 2 for the tag and 3 for the datasets
	TH1F *phi045_eta_unscaled[num_mBins][num_tBins][nDataSets]; 
	TH1F *phi090_eta_unscaled[num_mBins][num_tBins][nDataSets];
	TH1F *phi135_eta_unscaled[num_mBins][num_tBins][nDataSets];
	TH1F *phiAMO_eta_unscaled[num_mBins][num_tBins][nDataSets];
	TH1F *phi000_pi0_unscaled[num_mBins][num_tBins][nDataSets];
	TH1F *phi045_pi0_unscaled[num_mBins][num_tBins][nDataSets];
	TH1F *phi090_pi0_unscaled[num_mBins][num_tBins][nDataSets];
	TH1F *phi135_pi0_unscaled[num_mBins][num_tBins][nDataSets];
	TH1F *phiAMO_pi0_unscaled[num_mBins][num_tBins][nDataSets];
	for (int iData=0; iData <nDataSets; ++iData){
		string dataFileName = "/d/grid15/ln16/pi0eta/092419/degALL_"+dataSetTag[iData]+"_hists_DSelector.root";
		//string dataFileName = "/d/grid15/ln16/pi0eta/092419/newGraphs_histValues/rootFiles/deg000_data_"+dataSetTag[iData]+"_hists_DSelector.root";
		TFile *dataFile = new TFile(dataFileName.c_str());
		cout << "LOADING ROOT FILE: " << dataFileName << endl; 
		for (int iMass=0; iMass<num_mBins; ++iMass){
			for (int iteta=0; iteta<num_tBins; ++iteta){
				string affix = "Bin"+to_string(iteta)+"_Mpi0etaBin"+to_string(iMass);
				dataFile->GetObject(("prodPlanePSphi_045_teta"+affix).c_str(), phi045_eta_unscaled[iMass][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi_045_tpi0"+affix).c_str(), phi045_pi0_unscaled[iMass][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi_090_teta"+affix).c_str(), phi090_eta_unscaled[iMass][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi_090_tpi0"+affix).c_str(), phi090_pi0_unscaled[iMass][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi_AMO_teta"+affix).c_str(), phiAMO_eta_unscaled[iMass][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi_AMO_tpi0"+affix).c_str(), phiAMO_pi0_unscaled[iMass][iteta][iData]);
				// Need to scale these para yields by flux ratio.
				dataFile->GetObject(("prodPlanePSphi_000_teta"+affix).c_str(), phi000_eta_unscaled[iMass][iteta][iData]);
				cout << "phi000_eta_unscaled["<<iMass<<"]["<<iteta<<"]["<<iData<<"]" << endl;
				cout << (" -- prodPlanePSphi_000_teta"+affix).c_str() << endl;
				cout << " -- nentries=" << phi000_eta_unscaled[iMass][iteta][iData]->GetEntries() << endl;
				dataFile->GetObject(("prodPlanePSphi_000_tpi0"+affix).c_str(), phi000_pi0_unscaled[iMass][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi_135_teta"+affix).c_str(), phi135_eta_unscaled[iMass][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi_135_tpi0"+affix).c_str(), phi135_pi0_unscaled[iMass][iteta][iData]);
			}
		}
		
	}
	

	cout << "--------------------------------------------\nLOADING AND SCALING HISTOGRAMS\n--------------------------------------------" << endl;
	cout << "As a simple test we can GetMaximum before and after the scaling to see if it actually worked" << endl;
	cout << "\tWe will follow phi000_eta_total throughout the process" << endl;
	// These are the weighted summed histograms. weighted according to the Flux ratio
	TH1F *phi000_eta_total[num_mBins][num_tBins];
	TH1F *phi045_eta_total[num_mBins][num_tBins];
	TH1F *phi090_eta_total[num_mBins][num_tBins];
	TH1F *phi135_eta_total[num_mBins][num_tBins];
	TH1F *phiAMO_eta_total[num_mBins][num_tBins];
	TH1F *phi000_pi0_total[num_mBins][num_tBins];
	TH1F *phi045_pi0_total[num_mBins][num_tBins];
	TH1F *phi090_pi0_total[num_mBins][num_tBins];
	TH1F *phi135_pi0_total[num_mBins][num_tBins];
	TH1F *phiAMO_pi0_total[num_mBins][num_tBins];
	// Now we sum the scaled histograms over the 3 datasets
	double maximum_before;
	double maximum_after;
	int totalEntriesInPhi000_eta_total;
	for (int iMass=0; iMass<num_mBins; ++iMass){
		for (int iteta=0; iteta<num_tBins; ++iteta){
			totalEntriesInPhi000_eta_total=0;
			// We first have to clone the first dataTag (2017 run) and scale by flux ratio if in para config
			// we will scale the total in this to track the number of maxima and entires before and after scaling
			cout << "phi000_eta_unscaled["<<iMass<<"]["<<iteta<<"][0]" << endl;
			phi000_eta_total[iMass][iteta] = (TH1F *)phi000_eta_unscaled[iMass][iteta][0]->Clone();
			phi000_pi0_total[iMass][iteta] = (TH1F *)phi000_pi0_unscaled[iMass][iteta][0]->Clone();
			phi135_eta_total[iMass][iteta] = (TH1F *)phi135_eta_unscaled[iMass][iteta][0]->Clone();
			phi135_pi0_total[iMass][iteta] = (TH1F *)phi135_pi0_unscaled[iMass][iteta][0]->Clone();
			maximum_before = phi000_eta_total[iMass][iteta]->GetMaximum();
			cout << "-Maximum before scaling of phi000_eta_total[" << iteta << "]: " << maximum_before << endl;
			cout << "-NEntries before scaling of phi000_eta_total[" << iteta << "]: " << phi000_eta_total[iMass][iteta]->GetEntries() << endl;
			phi000_eta_total[iMass][iteta]->Scale(fluxRatios_90_0[0]);
			phi000_pi0_total[iMass][iteta]->Scale(fluxRatios_90_0[0]);
			phi135_eta_total[iMass][iteta]->Scale(fluxRatios_45_135[0]);
			phi135_pi0_total[iMass][iteta]->Scale(fluxRatios_45_135[0]);
			maximum_after = phi000_eta_total[iMass][iteta]->GetMaximum();
			cout << "-Maximum after scaling of phi000_eta_total[" << iteta << "]: " << maximum_after << endl;
			cout << "-NEntries after scaling of phi000_eta_total[" << iteta << "]: " << phi000_eta_total[iMass][iteta]->GetEntries() << endl; 
			cout << "-Scale factor vs Expected: " << maximum_after/maximum_before << " vs " << fluxRatios_90_0[0] << endl; 
			phi045_eta_total[iMass][iteta] = (TH1F *)phi045_eta_unscaled[iMass][iteta][0]->Clone();
			phi045_pi0_total[iMass][iteta] = (TH1F *)phi045_pi0_unscaled[iMass][iteta][0]->Clone();
			phi090_eta_total[iMass][iteta] = (TH1F *)phi090_eta_unscaled[iMass][iteta][0]->Clone();
			phi090_pi0_total[iMass][iteta] = (TH1F *)phi090_pi0_unscaled[iMass][iteta][0]->Clone();
			phiAMO_eta_total[iMass][iteta] = (TH1F *)phiAMO_eta_unscaled[iMass][iteta][0]->Clone();
			phiAMO_pi0_total[iMass][iteta] = (TH1F *)phiAMO_pi0_unscaled[iMass][iteta][0]->Clone();

			totalEntriesInPhi000_eta_total += phi000_eta_total[iMass][iteta]->GetEntries();
			cout << "-Cumulative entries after adding dataSet=0: " << totalEntriesInPhi000_eta_total << endl;
			// for the last 2 runs we simply add with a weight if in para config and add with weight=1 if in perp config
			for (int iData=1; iData < nDataSets; ++iData){
				// keeping a count of the total entries to check if we are doing the Add right
				totalEntriesInPhi000_eta_total += phi000_eta_unscaled[iMass][iteta][iData]->GetEntries();
				cout << "-Cumulative entries after adding dataSet=" << iData << ": " << totalEntriesInPhi000_eta_total << endl;

				// scale the rest of the data sets before adding them to the scaled+cloned 2017 dataset
				// we will never use unscaled histograms anymore so we can just scale it directly
				phi000_eta_unscaled[iMass][iteta][iData]->Scale(fluxRatios_90_0[iData]);
				phi000_pi0_unscaled[iMass][iteta][iData]->Scale(fluxRatios_90_0[iData]);
				phi135_eta_unscaled[iMass][iteta][iData]->Scale(fluxRatios_45_135[iData]);
				phi135_pi0_unscaled[iMass][iteta][iData]->Scale(fluxRatios_45_135[iData]);
				phi000_eta_total[iMass][iteta]->Add(phi000_eta_unscaled[iMass][iteta][iData]);
        			phi000_pi0_total[iMass][iteta]->Add(phi000_pi0_unscaled[iMass][iteta][iData]);
        			phi135_eta_total[iMass][iteta]->Add(phi135_eta_unscaled[iMass][iteta][iData]);
        			phi135_pi0_total[iMass][iteta]->Add(phi135_pi0_unscaled[iMass][iteta][iData]);

				// Dont have to scale these guys so just add them
        			phi045_eta_total[iMass][iteta]->Add(phi045_eta_unscaled[iMass][iteta][iData]);
        			phi045_pi0_total[iMass][iteta]->Add(phi045_pi0_unscaled[iMass][iteta][iData]);
        			phi090_eta_total[iMass][iteta]->Add(phi090_eta_unscaled[iMass][iteta][iData]);
        			phi090_pi0_total[iMass][iteta]->Add(phi090_pi0_unscaled[iMass][iteta][iData]);
        			phiAMO_eta_total[iMass][iteta]->Add(phiAMO_eta_unscaled[iMass][iteta][iData]);
        			phiAMO_pi0_total[iMass][iteta]->Add(phiAMO_pi0_unscaled[iMass][iteta][iData]);
			}
		}
	}
	// Just look at the lowest M(pi0eta) threshold and lowest t-bin for this MC test
	//

	for (int iBin=0; iBin<5; ++iBin){
		int numEntries000 = phi000_eta_total[0][iBin]->GetEntries();
		int numEntries090 = phi090_eta_total[0][iBin]->GetEntries();

		allCanvases->Clear();
		allCanvases->Divide(1,3);
		allCanvases->cd(1);
		phi000_eta_total[0][iBin]->Draw();
		allCanvases->cd(2);
		phi090_eta_total[0][iBin]->Draw();
		TH1* asymmetry000_090_eta_true = phi090_eta_total[0][iBin]->GetAsymmetry(phi000_eta_total[0][iBin]);
		allCanvases->cd(3);
		TF1 * fit_asym_true = new TF1(("fit_asym_true_tBin"+to_string(iBin)).c_str(),asymmetry,-180,180,numDOFsig_asym); 
		fit_asym_true->SetParameters(0.35,0.35,0.5,0);
		fit_asym_true->FixParameter(0,0.35);
		fit_asym_true->FixParameter(1,0.35);
		fit_asym_true->FixParameter(3,0);
		string fitOption = "E S";
		TFitResultPtr fitPointer = asymmetry000_090_eta_true->Fit(fit_asym_true,fitOption.c_str());
		asymmetry000_090_eta_true->Draw();
		double trueAsym = fit_asym_true->GetParameter(2);
		double trueAsymErr = fit_asym_true->GetParError(2);
		allCanvases->SaveAs(("asymmetryMC/true_asymmetry_tBin"+to_string(iBin)+".png").c_str());
		
		cout << numEntries000 << " entries are in the 000 orientation for eta" << endl;
		cout << numEntries090 << " entries are in the 090 orientation for eta" << endl;

		TF1* fit000 = new TF1("fit000",("1.0-0.35*"+to_string(trueAsym)+"*cos(2*0.01745*x)").c_str(),-180,180);
		TF1* fit090 = new TF1("fit000",("1.0-0.35*"+to_string(trueAsym)+"*cos(2*0.01745*(x-90))").c_str(),-180,180);
		std::vector<double> simulated000;
		std::vector<double> simulated090;
		double randValue000;
		double randValue090;

		std::vector<double> asymmetries;
		std::vector<double> asymmetryErrors;
		//TH1F* dHist_asymmetries = new TH1F(("dHist_asymmetries_tBin"+to_string(iBin)).c_str(),"dHist_asymmetries;#Sigma",100,0.16,0.36);
		TH1F* dHist_asymmetries = new TH1F(("dHist_asymmetries_tBin"+to_string(iBin)).c_str(),"dHist_asymmetries;#Sigma",50,0,0);
		TH1F* dHist_asymmetryErrors = new TH1F(("dHist_asymmetryErrors_tBin"+to_string(iBin)).c_str(),"dHist_asymmetryErrors;#sigma_{#Sigma}/#Sigma",50,0,0);
		TH1F* dHist_pulls = new TH1F(("dHist_pulls_tBin"+to_string(iBin)).c_str(), "dHist_pulls; #Pull", 50, 0, 0);
		double pull;
		dHist_asymmetries->GetXaxis()->SetTitleSize(0.06);
		dHist_asymmetryErrors->GetXaxis()->SetTitleSize(0.06);
		int maxPrint = 2;
		// We will just make double sure that rgen is doing the right thing for each orientation so we will just make 2
		cout << "Creating random number generators" << endl;
		TRandom3 rgen000;
		TRandom3 rgen090;
		rgen000.SetSeed(1992);
		rgen090.SetSeed(1992);
		double toyAsym;
		double toyAsymErr;
		for (int iter=0; iter<1000; ++iter){
			cout << "Sampling number of events for 000 with mean " << numEntries000 << endl;
			cout << "Sampling number of events for 090 with mean " << numEntries090 << endl;
			Int_t rNumEntries000 = rgen000.Poisson((double)numEntries000);	
			Int_t rNumEntries090 = rgen090.Poisson((double)numEntries090);	
			cout << "Using " << rNumEntries000 << " for 000" << endl;
			cout << "Using " << rNumEntries090 << " for 090" << endl;

			TH1F* dHist_simulated000 = new TH1F(("dHist_simulated000_tBin"+to_string(iBin)+"_"+to_string(iter)).c_str(),"dHist_simulated000; #phi degrees", 40,-180,180);
			TH1F* dHist_simulated090 = new TH1F(("dHist_simulated090_tBin"+to_string(iBin)+"_"+to_string(iter)).c_str(),"dHist_simulated090; #phi degrees", 40,-180,180);
			for (int ientry=0; ientry < rNumEntries000; ++ientry){
				randValue000 = fit000->GetRandom();
				simulated000.push_back(randValue000);
				dHist_simulated000->Fill(randValue000);
			}
			for (int ientry=0; ientry < rNumEntries090; ++ientry){
				randValue090 = fit090->GetRandom();
				simulated090.push_back(randValue090);
				dHist_simulated090->Fill(randValue090);
			}
			TH1* asymmetry000_090_eta = dHist_simulated090->GetAsymmetry(dHist_simulated000);
		
			TF1 * fit_sc000 = new TF1("fit_sc000",shiftedCos000,-180,180,numDOFsig_sc); 
			dHist_simulated000->Fit(fit_sc000);

			TF1 * fit_sc090 = new TF1("fit_sc090",shiftedCos090,-180,180,numDOFsig_sc); 
			dHist_simulated090->Fit(fit_sc090);

			TF1 * fit_asym = new TF1("fit_asym",asymmetry,-180,180,numDOFsig_asym); 
			fit_asym->SetParameters(0.35,0.35,0.5,0);
			fit_asym->FixParameter(0,0.35);
			fit_asym->FixParameter(1,0.35);
			fit_asym->FixParameter(3,0);
			string fitOption = "E S";
			TFitResultPtr fitPointer = asymmetry000_090_eta->Fit(fit_asym,fitOption.c_str());

			toyAsym = fit_asym->GetParameter(2);
			toyAsymErr = abs(fit_asym->GetParError(2)); 

			asymmetries.push_back(toyAsym);
			asymmetryErrors.push_back(toyAsymErr);
			dHist_asymmetries->Fill(toyAsym);
			dHist_asymmetryErrors->Fill(toyAsymErr/toyAsym);

			pull = (trueAsym-toyAsym)/sqrt(trueAsymErr*trueAsymErr-toyAsymErr*toyAsymErr);
			dHist_pulls->Fill(pull);

			if (iter<maxPrint){
				allCanvases->Clear();
				allCanvases->Divide(2,2);
				allCanvases->cd(1);
				dHist_simulated000->Draw("E");
				dHist_simulated000->GetXaxis()->SetTitleSize(0.06);
				allCanvases->cd(2);
				dHist_simulated090->Draw("E");
				dHist_simulated090->GetXaxis()->SetTitleSize(0.06);
				allCanvases->cd(3);
				asymmetry000_090_eta->Draw("E");
				asymmetry000_090_eta->GetXaxis()->SetTitleSize(0.06);
				allCanvases->cd(4);
				allCanvases->SaveAs(("asymmetryMC/simulatedAsymmetry_tBin"+to_string(iBin)+"_set"+to_string(iter)+".png").c_str());
			}
		}

		allCanvases->Clear();
		allCanvases->Divide(2,2);
		allCanvases->cd(1);
		TLine* newLine = new TLine();
		newLine->SetLineColor(kRed);
		dHist_asymmetries->Draw();
		newLine->DrawLine(trueAsym,0,trueAsym,dHist_asymmetries->GetMaximum());
		cout << "Drawing a line at " << trueAsym << " for the Asymmetry" << endl;
		allCanvases->cd(2);
		dHist_asymmetryErrors->Draw();
		newLine->DrawLine(trueAsymErr/trueAsym,0,trueAsymErr/trueAsym,dHist_asymmetries->GetMaximum());
		cout << "Drawing a line at " << trueAsymErr/trueAsym << " for the Asymmetry fraction unceratinty" << endl;
		dHist_asymmetryErrors->SetTitle(("True Fractional Uncertainty = "+to_string(trueAsymErr/trueAsym)).c_str());
		allCanvases->cd(3);
		TF1* pullFit = new TF1("pullFit","gausn",-5,5); 
		dHist_pulls->Fit(pullFit,"S");
		double pullMean = pullFit->GetParameter(1);
		double pullMeanSig = pullFit->GetParError(1);
		double pullSig = pullFit->GetParameter(2);
		double pullSigSig = pullFit->GetParError(2);
		dHist_pulls->Draw();
		dHist_pulls->SetTitle(("#splitline{Mean, MeanSig, Sigma, SigmaSig = "+to_string(pullMean)+", "+to_string(pullMeanSig)+"}{"+to_string(pullSig)+", "+to_string(pullSigSig)+"}").c_str());
		dHist_pulls->GetXaxis()->SetTitleSize(0.06);

		allCanvases->SaveAs(("asymmetryMC/simulatedSigma_tBin"+to_string(iBin)+"_SigErr.png").c_str());
		dHist_asymmetries->SaveAs(("asymmetryMC/Asymmetries_tBin"+to_string(iBin)+".C").c_str());
		dHist_asymmetryErrors->SaveAs(("asymmetryMC/AsymmetryErrors_tBin"+to_string(iBin)+".C").c_str());
	}
}

