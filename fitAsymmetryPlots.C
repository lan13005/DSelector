double degToRad=TMath::Pi()/180;
int numDOFsig_sc = 3;
Double_t shiftedCos(Double_t *x, Double_t *par){
	return par[0]*(1.0 + par[1]*TMath::Cos(2*degToRad*(x[0]-par[2])));
}

int numDOFsig_flat = 1;
Double_t flat(Double_t *x, Double_t *par){
	return par[0];
}

int numDOFsig_asym=4;
Double_t asymmetry(Double_t *x, Double_t *par){
	return ((par[0]+par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3])))/(2+(par[0]-par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3])));
}

void fitAsymmetryPlots(){
	// seems like this this would somehow instantiate gMinuit maybe? If I dont do this I get some errors
	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); // If I dont get this I get the error "Error in <TMinuitMinimizer::Scan>:  Error executing command SCAN" 
	static const int num_tBins=5;

	gStyle->SetOptFit();
	gStyle->SetStatY(1);
	gStyle->SetStatX(1);
	gStyle->SetStatW(0.16);
	gStyle->SetStatH(0.16);

	TCanvas *allCanvases = new TCanvas("","",1440,900);


	// *****************************
	// Define flux ratios for 2017, 2018_1, 2018_2 
	// *****************************
	double fluxRatios_90_0[3] = {  4.346818e+12/4.188001e+12, 0.965429, 0.918503 };
	double fluxRatios_45_135[3] = {  4.076065e+12/4.095013e+12, 1.02261, 1.03254 };
	string dataSetTag[3] = { "2017", "2018_1", "2018_8" };

	string tagEta[2] = {"","_backwardPi0P"};
	string tagPi0[2] = {"","_backwardEtaP"};
	string tag[2] = {"","vanHoveSelected"};
	std::vector<int> nEventsPhiEta[2]; // second one is for the vanHove selected
	std::vector<int> nEventsPhiPi0[2];

	// *****************************
	// Loading histograms for  2017, 2018_1, 2018_2, then scaling yield by the flux ratio. Then we can finally add the data sets together
	// *****************************
	TH1F *phi000_eta_unscaled[2][num_tBins][3]; // 2 for the tag and 3 for the datasets
	TH1F *phi045_eta_unscaled[2][num_tBins][3];
	TH1F *phi090_eta_unscaled[2][num_tBins][3];
	TH1F *phi135_eta_unscaled[2][num_tBins][3];
	TH1F *phiAMO_eta_unscaled[2][num_tBins][3];
	TH1F *phi000_pi0_unscaled[2][num_tBins][3];
	TH1F *phi045_pi0_unscaled[2][num_tBins][3];
	TH1F *phi090_pi0_unscaled[2][num_tBins][3];
	TH1F *phi135_pi0_unscaled[2][num_tBins][3];
	TH1F *phiAMO_pi0_unscaled[2][num_tBins][3];
	TH1F *phi000_eta_total[2][num_tBins];
	TH1F *phi045_eta_total[2][num_tBins];
	TH1F *phi090_eta_total[2][num_tBins];
	TH1F *phi135_eta_total[2][num_tBins];
	TH1F *phiAMO_eta_total[2][num_tBins];
	TH1F *phi000_pi0_total[2][num_tBins];
	TH1F *phi045_pi0_total[2][num_tBins];
	TH1F *phi090_pi0_total[2][num_tBins];
	TH1F *phi135_pi0_total[2][num_tBins];
	TH1F *phiAMO_pi0_total[2][num_tBins];
	for (int iData=0; iData < sizeof(fluxRatios_90_0)/sizeof(fluxRatios_90_0[0]); ++iData){
		string dataFileName = "deg000_data_"+dataSetTag[iData]+"_hists_DSelector.root";
		TFile *dataFile = new TFile(dataFileName.c_str());
		cout << "LOADING ROOT FILE: " << dataFileName << endl; 
		for (int iTag=0; iTag < sizeof(tagEta)/sizeof(tagEta[0]); ++iTag){
			for (int iteta=0; iteta<num_tBins; ++iteta){
				dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_045_tetaBin"+to_string(iteta)).c_str(), phi045_eta_unscaled[iTag][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_045_tpi0Bin"+to_string(iteta)).c_str(), phi045_pi0_unscaled[iTag][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_090_tetaBin"+to_string(iteta)).c_str(), phi090_eta_unscaled[iTag][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_090_tpi0Bin"+to_string(iteta)).c_str(), phi090_pi0_unscaled[iTag][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_AMO_tetaBin"+to_string(iteta)).c_str(), phiAMO_eta_unscaled[iTag][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_AMO_tpi0Bin"+to_string(iteta)).c_str(), phiAMO_pi0_unscaled[iTag][iteta][iData]);
				// Need to scale these para yields by flux ratio.
				dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_000_tetaBin"+to_string(iteta)).c_str(), phi000_eta_unscaled[iTag][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_000_tpi0Bin"+to_string(iteta)).c_str(), phi000_pi0_unscaled[iTag][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi"+tagEta[iTag]+"_135_tetaBin"+to_string(iteta)).c_str(), phi135_eta_unscaled[iTag][iteta][iData]);
				dataFile->GetObject(("prodPlanePSphi"+tagPi0[iTag]+"_135_tpi0Bin"+to_string(iteta)).c_str(), phi135_pi0_unscaled[iTag][iteta][iData]);
			}
		}
	}
	// Now we sum the scaled histograms over the 3 datasets
	double maximum_before;
	double maximum_after;
	int totalEntriesInPhi000_eta_total;
	cout << "--------------------------------------------\nLOADING AND SCALING HISTOGRAMS\n--------------------------------------------" << endl;
	cout << "As a simple test we can GetMaximum before and after the scaling to see if it actually worked" << endl;
	cout << "\tWe will follow phi000_eta_total throughout the process" << endl;
	for (int iTag=0; iTag < sizeof(tagEta)/sizeof(tagEta[0]); ++iTag){
		for (int iteta=0; iteta<num_tBins; ++iteta){
			totalEntriesInPhi000_eta_total=0;
			cout << "--------------\n" << tagEta[0] << "\n--------------" << endl;
			// We first have to clone the first dataTag (2017 run) and scale by flux ratio if in para config
			phi000_eta_total[iTag][iteta] = (TH1F *)phi000_eta_unscaled[iTag][iteta][0]->Clone();
			phi000_pi0_total[iTag][iteta] = (TH1F *)phi000_pi0_unscaled[iTag][iteta][0]->Clone();
			phi135_eta_total[iTag][iteta] = (TH1F *)phi135_eta_unscaled[iTag][iteta][0]->Clone();
			phi135_pi0_total[iTag][iteta] = (TH1F *)phi135_pi0_unscaled[iTag][iteta][0]->Clone();
			maximum_before = phi000_eta_total[iTag][iteta]->GetMaximum();
			cout << "-Maximum before scaling of phi000_eta_total[" << iTag << "][" << iteta << "]: " << maximum_before << endl;
			cout << "-NEntries before scaling of phi000_eta_total[" << iTag << "]" << iteta << "]: " << phi000_eta_total[iTag][iteta]->GetEntries() << endl;
			phi000_eta_total[iTag][iteta]->Scale(fluxRatios_90_0[0]);
			phi000_pi0_total[iTag][iteta]->Scale(fluxRatios_90_0[0]);
			phi135_eta_total[iTag][iteta]->Scale(fluxRatios_45_135[0]);
			phi135_pi0_total[iTag][iteta]->Scale(fluxRatios_45_135[0]);
			maximum_after = phi000_eta_total[iTag][iteta]->GetMaximum();
			cout << "-Maximum after scaling of phi000_eta_total[" << iTag << "]" << iteta << "]: " << maximum_after << endl;
			cout << "-NEntries after scaling of phi000_eta_total[" << iTag << "]" << iteta << "]: " << phi000_eta_total[iTag][iteta]->GetEntries() << endl; 
			cout << "-Scale factor vs Expected: " << maximum_after/maximum_before << " vs " << fluxRatios_90_0[0] << endl; 
			phi045_eta_total[iTag][iteta] = (TH1F *)phi045_eta_unscaled[iTag][iteta][0]->Clone();
			phi045_pi0_total[iTag][iteta] = (TH1F *)phi045_pi0_unscaled[iTag][iteta][0]->Clone();
			phi090_eta_total[iTag][iteta] = (TH1F *)phi090_eta_unscaled[iTag][iteta][0]->Clone();
			phi090_pi0_total[iTag][iteta] = (TH1F *)phi090_pi0_unscaled[iTag][iteta][0]->Clone();
			phiAMO_eta_total[iTag][iteta] = (TH1F *)phiAMO_eta_unscaled[iTag][iteta][0]->Clone();
			phiAMO_pi0_total[iTag][iteta] = (TH1F *)phiAMO_pi0_unscaled[iTag][iteta][0]->Clone();

			totalEntriesInPhi000_eta_total += phi000_eta_total[iTag][iteta]->GetEntries();
			cout << "-Cumulative entries after adding dataSet=0: " << totalEntriesInPhi000_eta_total << endl;
			// for the last 2 runs we simply add with a weight if in para config and add with weight=1 if in perp config
			for (int iData=1; iData < sizeof(fluxRatios_90_0)/sizeof(fluxRatios_90_0[0]); ++iData){
				phi000_eta_total[iTag][iteta]->Add(phi000_eta_unscaled[iTag][iteta][iData]);
				totalEntriesInPhi000_eta_total += phi000_eta_unscaled[iTag][iteta][iData]->GetEntries();
				cout << "-Cumulative entries after adding dataSet=" << iData << ": " << totalEntriesInPhi000_eta_total << endl;
                		phi000_pi0_total[iTag][iteta]->Add(phi000_pi0_unscaled[iTag][iteta][iData]);
                		phi045_eta_total[iTag][iteta]->Add(phi045_eta_unscaled[iTag][iteta][iData]);
                		phi045_pi0_total[iTag][iteta]->Add(phi045_pi0_unscaled[iTag][iteta][iData]);
                		phi090_eta_total[iTag][iteta]->Add(phi090_eta_unscaled[iTag][iteta][iData]);
                		phi090_pi0_total[iTag][iteta]->Add(phi090_pi0_unscaled[iTag][iteta][iData]);
                		phi135_eta_total[iTag][iteta]->Add(phi135_eta_unscaled[iTag][iteta][iData]);
                		phi135_pi0_total[iTag][iteta]->Add(phi135_pi0_unscaled[iTag][iteta][iData]);
                		phiAMO_eta_total[iTag][iteta]->Add(phiAMO_eta_unscaled[iTag][iteta][iData]);
                		phiAMO_pi0_total[iTag][iteta]->Add(phiAMO_pi0_unscaled[iTag][iteta][iData]);
			}
		}
	}	
	

	cout << "--------------------------------------------\nSTARTING FITTING\n--------------------------------------------" << endl;

	// *****************************
	// Defining some variables that we will use
	// *****************************
	string names[4] = {"phi000","phi045","phi090","phi135"};
	double asymmetries_000_eta[num_tBins];
	double asymmetries_045_eta[num_tBins];
	double asymmetries_000_eta_err[num_tBins];
	double asymmetries_045_eta_err[num_tBins];
	double asymmetries_000_pi0[num_tBins];
	double asymmetries_045_pi0[num_tBins];
	double asymmetries_000_pi0_err[num_tBins];
	double asymmetries_045_pi0_err[num_tBins];
	double tBins[num_tBins];
	double tBins_err[num_tBins];
	double tBinSize=0.2;
	double orientation[4] = {0,45,90,135};


	string phis_orientation[5] =  {"phi000", "phi045", "phi090", "phi135", "phiAMO"};
	for (int iTag=0; iTag < sizeof(tagEta)/sizeof(tagEta[0]); ++iTag){
		for (int iteta=0; iteta<num_tBins; ++iteta){ //num_tBins
			// Need to define the polarization with the systematic shift
			double Phi0_0_90 = 3.1;
			double Phi0_45_135 = 3.2+45;
			string fitOption = "E S";
			TH1F *phi000_eta = phi000_eta_total[iTag][iteta];
			TH1F *phi045_eta = phi045_eta_total[iTag][iteta];
			TH1F *phi090_eta = phi090_eta_total[iTag][iteta];
			TH1F *phi135_eta = phi135_eta_total[iTag][iteta];
			TH1F *phiAMO_eta = phiAMO_eta_total[iTag][iteta];
			TH1F *phi000_pi0 = phi000_pi0_total[iTag][iteta];
			TH1F *phi045_pi0 = phi045_pi0_total[iTag][iteta];
			TH1F *phi090_pi0 = phi090_pi0_total[iTag][iteta];
			TH1F *phi135_pi0 = phi135_pi0_total[iTag][iteta];
			TH1F *phiAMO_pi0 = phiAMO_pi0_total[iTag][iteta];

			// *****************************
			// Some definitions for Eta
			// *****************************
			TH1F *phis_eta[4] =  {phi000_eta, phi045_eta, phi090_eta, phi135_eta};
			TH1* asymmetry000_090_eta = phi000_eta->GetAsymmetry(phi090_eta);
			TH1* asymmetry045_135_eta = phi045_eta->GetAsymmetry(phi135_eta);
			asymmetry000_090_eta->SetTitle("0/90 Asymmetry");
			asymmetry045_135_eta->SetTitle("45/135 Asymmetry");
			// *****************************
			// Some definitions for Pi0
			// *****************************
			TH1F *phis_pi0[4] =  {phi000_pi0, phi045_pi0, phi090_pi0, phi135_pi0};
			TH1* asymmetry000_090_pi0 = phi000_pi0->GetAsymmetry(phi090_pi0);
			TH1* asymmetry045_135_pi0 = phi045_pi0->GetAsymmetry(phi135_pi0);
			asymmetry000_090_pi0->SetTitle(("0/90 Asymmetry"+tagEta[iTag]).c_str());
			asymmetry045_135_pi0->SetTitle(("45/135 Asymmetry"+tagPi0[iTag]).c_str());

			// *****************************
			// Fitting asymmetry for eta
			// *****************************
			allCanvases->Clear();
			allCanvases->Divide(2,1);
			allCanvases->cd(1);
			Int_t fitStatus;
			// We set P_perp = P_para = 0.35 which is close to the expected. We initialize asymmetry to be 0 since it can vary from [-1,1]. 
			TF1 * fit_asym = new TF1("fit_asym",asymmetry,-180,180,numDOFsig_asym); 
			fit_asym->SetParameters(0.35,0.35,0.5,0);
			fit_asym->FixParameter(0,0.35);
			fit_asym->FixParameter(1,0.35);
			fit_asym->FixParameter(3,Phi0_0_90);
			allCanvases->Clear();

			TFitResultPtr fitPointer = asymmetry000_090_eta->Fit(fit_asym,fitOption.c_str());
			asymmetries_000_eta[iteta] = fit_asym->GetParameter(2);
			asymmetries_000_eta_err[iteta] = fit_asym->GetParError(2);
			asymmetry000_090_eta->Draw("SAME");
			asymmetry000_090_eta->SetAxisRange(-0.5,0.6,"Y");
			TGraph* likelihoodFit_000_090_eta = new TGraph(50); 
			fitPointer->Scan(2,likelihoodFit_000_090_eta,-2,2);

			allCanvases->cd(2);
			fit_asym->SetParameters(0.35,0.35,0.5,Phi0_45_135);
			fit_asym->FixParameter(0,0.35);
			fit_asym->FixParameter(1,0.35);
			fit_asym->FixParameter(3,Phi0_45_135);

			fitPointer = asymmetry045_135_eta->Fit(fit_asym,fitOption.c_str());
			asymmetries_045_eta[iteta] = fit_asym->GetParameter(2);
			asymmetries_045_eta_err[iteta] = fit_asym->GetParError(2);
			asymmetry045_135_eta->Draw("SAME");
			asymmetry045_135_eta->SetAxisRange(-0.5,0.6,"Y");
			TGraph* likelihoodFit_045_135_eta = new TGraph(50); 
			fitPointer->Scan(2,likelihoodFit_045_135_eta,-2,2);

			allCanvases->SaveAs(("asymmetryPlots/asymmetry"+tagEta[iTag]+"_tetaBin"+to_string(iteta)+".png").c_str());
			
			allCanvases->Clear();
			allCanvases->Divide(2,1);
			allCanvases->cd(1);
			likelihoodFit_000_090_eta->Draw("ALP");
			likelihoodFit_000_090_eta->SetTitle(("0/90 - tBin"+to_string(iteta)).c_str());
			allCanvases->cd(2);
			likelihoodFit_045_135_eta->Draw("ALP");
			likelihoodFit_045_135_eta->SetTitle(("045_135 - tBin"+to_string(iteta)).c_str());
			allCanvases->SaveAs(("asymmetryPlots/likelihood"+tagEta[iTag]+"_tetaBin"+to_string(iteta)+".png").c_str());

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
			asymmetries_000_pi0[iteta] = fit_asym->GetParameter(2);
			asymmetries_000_pi0_err[iteta] = fit_asym->GetParError(2);
			asymmetry000_090_pi0->Draw("SAME");
			asymmetry000_090_pi0->SetAxisRange(-0.5,0.6,"Y");
			TGraph* likelihoodFit_000_090_pi0 = new TGraph(50); 
			fitPointer->Scan(2,likelihoodFit_000_090_pi0,-2,2);
	
			allCanvases->cd(2);
			fit_asym->SetParameters(0.35,0.35,0.5,Phi0_45_135);
			fit_asym->FixParameter(0,0.35);
			fit_asym->FixParameter(1,0.35);
			fit_asym->FixParameter(3,Phi0_45_135);
			fitPointer = asymmetry045_135_pi0->Fit(fit_asym,fitOption.c_str());
			asymmetry045_135_pi0->Draw("SAME");
			asymmetry045_135_pi0->SetAxisRange(-0.5,0.6,"Y");
			asymmetries_045_pi0[iteta] = fit_asym->GetParameter(2);
			asymmetries_045_pi0_err[iteta] = fit_asym->GetParError(2);
			TGraph* likelihoodFit_045_135_pi0 = new TGraph(50); 
			fitPointer->Scan(2,likelihoodFit_045_135_pi0,-2,2);

			allCanvases->SaveAs(("asymmetryPlots/asymmetry"+tagPi0[iTag]+"_tpi0Bin"+to_string(iteta)+".png").c_str());

			allCanvases->Clear();
			allCanvases->Divide(2,1);
			allCanvases->cd(1);
			likelihoodFit_000_090_eta->Draw("ALP");
			likelihoodFit_000_090_eta->SetTitle(("0/90 - tBin"+to_string(iteta)).c_str());
			allCanvases->cd(2);
			likelihoodFit_045_135_eta->Draw("ALP");
			likelihoodFit_045_135_eta->SetTitle(("045/135 - tBin"+to_string(iteta)).c_str());
			allCanvases->SaveAs(("asymmetryPlots/likelihood"+tagEta[iTag]+"_tetaBin"+to_string(iteta)+".png").c_str());


			// *****************************
			// Defining the tBins which will be shared
			// *****************************
			tBins[iteta] = iteta*tBinSize;
			tBins_err[iteta] = tBinSize/2;
	
			// *****************************
			// Fitting prodPlanePhi for eta
			// *****************************
			allCanvases->Clear();
			allCanvases->Divide(3,2);
			TF1 * fit_sc = new TF1("fit_sc",shiftedCos,-180,180,numDOFsig_sc); 
			double p0;
			int counter=0;
			for (auto phi: phis_eta){
				allCanvases->cd(counter+1);
				phi->SetTitle((names[counter]).c_str());
				nEventsPhiEta[iTag].push_back(phi->GetEntries());
				p0 = phi->GetEntries()/phi->GetNbinsX();
				fit_sc->SetParameters(p0,0.1,orientation[counter]);
				fit_sc->FixParameter(2, orientation[counter]);
				fitStatus = phi->Fit(fit_sc,"RLQE");
				phi->Draw("SAME");
				++counter;
			}
			TF1 * fit_flat = new TF1("fit_flat",flat,-180,180,numDOFsig_flat); 
			allCanvases->cd(5);
			phiAMO_eta->SetTitle("phiAMO");
			p0 = phiAMO_eta->GetEntries()/phiAMO_eta->GetNbinsX();
			nEventsPhiEta[iTag].push_back(phiAMO_eta->GetEntries());
			fit_sc->SetParameters(p0,p0/2,orientation[counter]);
			fit_sc->FixParameter(2, orientation[counter]);
			fitStatus = phiAMO_eta->Fit(fit_flat,"RLQE");
			phiAMO_eta->Draw("SAME");
			allCanvases->SaveAs(("asymmetryPlots/phiYieldFits"+tagEta[iTag]+"_tetaBin"+to_string(iteta)+".png").c_str());
			
	
			// *****************************
			// Fitting prodPlanePhi for pi0
			// *****************************
			allCanvases->Clear();
			allCanvases->Divide(3,2);
			fit_sc = new TF1("fit_sc",shiftedCos,-180,180,numDOFsig_sc); 
			counter=0;
			for (auto phi: phis_pi0){
				allCanvases->cd(counter+1);
				phi->SetTitle((names[counter]).c_str());
				p0 = phi->GetEntries()/phi->GetNbinsX();
				fit_sc->SetParameters(p0,0.1,orientation[counter]);
				fit_sc->FixParameter(2, orientation[counter]);
				fitStatus = phi->Fit(fit_sc,"RLQE");
				phi->Draw("SAME");
				nEventsPhiPi0[iTag].push_back(phi->GetEntries());
				cout << "(iTag=" << iTag << ")Entries in nEventsPhiPi0 if different orientations: " << phi->GetEntries() << endl;
				++counter;
			}
			fit_flat = new TF1("fit_flat",flat,-180,180,numDOFsig_flat); 
			allCanvases->cd(5);
			phiAMO_pi0->SetTitle("phiAMO");
			nEventsPhiPi0[iTag].push_back(phiAMO_pi0->GetEntries());
			p0 = phiAMO_pi0->GetEntries()/phiAMO_pi0->GetNbinsX();
			fit_sc->SetParameters(p0,0.1,orientation[counter]);
			fit_sc->FixParameter(2, orientation[counter]);
			fitStatus = phiAMO_pi0->Fit(fit_flat,"RLQE");
			phiAMO_pi0->Draw("SAME");
			cout << "(iTag=" << iTag << ")Entries in nEventsPhiPi0 if different orientations: " << phiAMO_pi0->GetEntries() << endl;
			allCanvases->SaveAs(("asymmetryPlots/phiYieldFits"+tagPi0[iTag]+"_tpi0Bin"+to_string(iteta)+".png").c_str());
		}
	
		// *****************************
		// Now we overlay the pi0 and eta beam asymmetries 
		// *****************************
		allCanvases->Clear();
		allCanvases->Divide(2,1);
		allCanvases->cd(1);
		auto gr_000 = new TGraphErrors(num_tBins,tBins,asymmetries_000_eta,tBins_err,asymmetries_000_eta_err);
		gr_000->SetTitle("Beam Asymmetry 0/90 Orientation");
		gr_000->SetMarkerColor(4);
		gr_000->SetMarkerStyle(21);
		gr_000->SetLineColor(4);
		gr_000->Draw("AP");
		gr_000->GetHistogram()->SetMaximum(1.2);
		gr_000->GetHistogram()->SetMinimum(-1);
		gr_000 = new TGraphErrors(num_tBins,tBins,asymmetries_000_pi0,tBins_err,asymmetries_000_pi0_err);
		gr_000->SetTitle("Beam Asymmetry 0/90 Orientation");
		gr_000->SetLineColor(2);
		gr_000->SetMarkerColor(2);
		gr_000->SetMarkerStyle(20);
		gr_000->Draw("P SAME");
	
		allCanvases->cd(2);
		auto gr_045 = new TGraphErrors(num_tBins,tBins,asymmetries_045_eta,tBins_err,asymmetries_045_eta_err);
		gr_045->SetTitle("Beam Asymmetry 45/135 Orientation");
		gr_045->SetMarkerColor(4);
		gr_045->SetLineColor(4);
		gr_045->SetMarkerStyle(21);
		gr_045->Draw("AP");
		gr_045->GetHistogram()->SetMaximum(1.2);
		gr_045->GetHistogram()->SetMinimum(-1);
		gr_045 = new TGraphErrors(num_tBins,tBins,asymmetries_045_pi0,tBins_err,asymmetries_045_pi0_err);
		gr_045->SetTitle("Beam Asymmetry 45/135 Orientation");
		gr_045->SetMarkerColor(2);
		gr_045->SetLineColor(2);
		gr_045->SetMarkerStyle(20);
		gr_045->Draw("P SAME");
		allCanvases->SaveAs(("asymmetryPlots/asymVst"+tag[iTag]+".png").c_str());
	}



	cout << "\n" << endl;
	cout << "Number entries in nEventsPhiEta: " << nEventsPhiEta[0].size() << endl;
	cout << "Number entries in nEventsPhiEta VH selected: " << nEventsPhiEta[1].size() << endl;
	cout << "Number entries in nEventsPhiPi0: " << nEventsPhiPi0[0].size() << endl;
	cout << "Number entries in nEventsPhiPi0 VH selected: " << nEventsPhiPi0[1].size() << endl;
	int currentIndex;
	int sizeOfOrientationList = sizeof(phis_orientation)/sizeof(phis_orientation[0]);
	for (int iteta=0; iteta<num_tBins; ++iteta){
		cout << "--------------------------------------" << endl;
		for (int phiIndex=0; phiIndex<sizeOfOrientationList; ++phiIndex) { 
			currentIndex = iteta*sizeOfOrientationList + phiIndex;
			cout << "CURRENT INDEX: " << currentIndex << endl;
			cout << "(" << phis_orientation[phiIndex] << ")" << "Number events in teta Bin " << iteta << ": " << (double)nEventsPhiEta[0][currentIndex] << endl;
			cout << "(" << phis_orientation[phiIndex] << ")" << "Number VH events in teta Bin " << iteta << ": " << (double)nEventsPhiEta[1][currentIndex] << endl;
			cout << "(" << phis_orientation[phiIndex] << ")" << "VH selection eff in teta Bin " << iteta << ": " << (double)nEventsPhiEta[1][currentIndex]/nEventsPhiEta[0][currentIndex] << endl;
			cout << "(" << phis_orientation[phiIndex] << ")" << "Number events in tpi0 Bin " << iteta << ": " << (double)nEventsPhiPi0[0][currentIndex] << endl;
			cout << "(" << phis_orientation[phiIndex] << ")" << "Number VH events in tpi0 Bin " << iteta << ": " << (double)nEventsPhiPi0[1][currentIndex] << endl;
			cout << "(" << phis_orientation[phiIndex] << ")" << "VH selection eff in tpi0 Bin " << iteta << ": " << (double)nEventsPhiPi0[1][currentIndex]/nEventsPhiPi0[0][currentIndex] << endl;
		}
	}
	cout << "\n\n" << endl;
}
















