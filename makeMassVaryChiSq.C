
// 1d double gaussian
int numDOFsig=5;
Double_t signal(Double_t *x, Double_t *par){
	return par[0]/par[2]/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]))
	     + par[3]*par[0]/(par[4]*par[2])/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/(par[4]*par[2]))*((x[0]-par[1])/(par[4]*par[2])));

}

int numDOFbkg = 2;
Double_t background(Double_t *x, Double_t *par){
	return par[0]+par[1]*x[0];
}

Double_t fitFunc(Double_t *x, Double_t *par){
	return background(x,par)+signal(x,&par[numDOFbkg]);
}

void makeMassVaryChiSq(){
	TFile* dataFile = TFile::Open("pi0eta_test_hists_DSelector.root");
    	//ofstream logFile;
    	//logFile.open("newGraphs/sigEffs/purityResults.txt");

	gStyle->SetOptStat(0);

	TCanvas *allCanvases_puritys = new TCanvas("allCanvases_puritys","",1440,900);
	allCanvases_puritys->Divide(2,1);

	TH1F *anyHist;
	TH1F *pi0Hist;
	string particles[2] = {"eta","pi0"};
	string baseName = "Mass_Kin_mEllipsePre_ChiMpi0etaBin";
	gStyle->SetTitleSize(0.08,"t");

	std::vector<double> fitRangeEta;
	std::vector<double> fitRangePi0;
	fitRangeEta={0.44,0.64};
	fitRangePi0={0.11,0.165};

	double maxY = DBL_MIN;
	double minY = DBL_MAX;
	TF1* fit;
	TF1* bkgFit;
	TF1* sigFit;
	Double_t par[numDOFbkg+numDOFsig];

	double minMpi0eta = 0.6;	
	double maxMpi0eta = 2.2;
	int numMpi0etaBins = 8;
	double binStepMpi0eta = (maxMpi0eta-minMpi0eta)/numMpi0etaBins;
	static const int numRegions_ChiSq=7;
	double integralBKGs[2][numRegions_ChiSq*numMpi0etaBins];
	double integralSIGs[2][numRegions_ChiSq*numMpi0etaBins];
	double puritys[2][numRegions_ChiSq*numMpi0etaBins];
	double chiSqs[2][numRegions_ChiSq*numMpi0etaBins]; // from the individual fits
	double integralBKG;
	double integralSIG;
	double purity;
	int cumulativeBinNum;
	double chiSq;

	double chiSqThresholds[numRegions_ChiSq] = {1,3,5,7,9,11,13}; // to use for getting the names of the histograms
	cout << "Starting to loop through all the histograms" << endl;
	TMultiGraph *mg_eta = new TMultiGraph("mg_eta","Purities In Eta");
	TMultiGraph *mg_pi0 = new TMultiGraph("mg_pi0","Purities In Pi0");
	auto legend = new TLegend(0.1,0.7,0.48,0.9);
	std::vector<TGraph*> gr_eta;
	std::vector<TGraph*> gr_pi0;
	for (int bin_Mpi0eta=0; bin_Mpi0eta<numMpi0etaBins; ++ bin_Mpi0eta){
		cout << "  Beginning Mpi0eta bin: " << bin_Mpi0eta << endl;
		double graph_integralBKGs[2][numRegions_ChiSq];
		double graph_integralSIGs[2][numRegions_ChiSq];
		double graph_puritys[2][numRegions_ChiSq];
		double graph_chiSqs[2][numRegions_ChiSq];
        	for (int bin=0; bin<numRegions_ChiSq; ++bin){
			cout << "    Begin work in ChiSq Bin: " << bin << endl;
			int particleCounter = 0;
			TCanvas *allCanvases_fits = new TCanvas("allCanvases_fits","",1440,900);
			allCanvases_fits->Divide(2,1);
			for ( auto particle : particles){
				cumulativeBinNum = bin_Mpi0eta*numRegions_ChiSq+bin;
				cout << "      Begin work with with: " << particle << " at cumulativeBinNum:" << cumulativeBinNum << endl;
				string name = particle+baseName+to_string(cumulativeBinNum);
				cout << "Getting object: " << name << endl;
				dataFile->GetObject(name.c_str(), anyHist);
				if (particleCounter==0){
					fit = new TF1("fit",fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFsig+numDOFbkg); 
					bkgFit = new TF1("bkgFit",background,fitRangeEta[0],fitRangeEta[1],numDOFbkg); 
					sigFit = new TF1("sigFit",signal,fitRangeEta[0],fitRangeEta[1],numDOFsig); 
				}
				else{
					fit = new TF1("fit",fitFunc,fitRangePi0[0],fitRangePi0[1],numDOFsig+numDOFbkg); 
					bkgFit = new TF1("bkgFit",background,fitRangePi0[0],fitRangePi0[1],numDOFbkg); 
					sigFit = new TF1("sigFit",signal,fitRangePi0[0],fitRangePi0[1],numDOFsig); 
				}

				double par0;
				double par1;
				double par2;

				int maxAttempts = 200;
				int nSig = 3;
				double weightedSigma;
				TLine *drawFitBounds = new TLine();
				for (int iAttempt=0; iAttempt < maxAttempts; ++iAttempt){
					//if (anyHist->GetEntries() < 40){ 
					//	cout << "------------- TOO LITTLE EVENTS SETTING ALL RESULTS TO ZERO" << endl;
					//	integralBKG=0;
					//	integralSIG=0;
					//	purity=0;
					//	chiSq=0;
					//	break; 
					//}
					int totalEntries = anyHist->GetEntries();
					int scaledEntries = totalEntries/anyHist->GetNbinsX();
					par0 = rand()%scaledEntries;
					par1 = rand()%scaledEntries;
					par2 = rand()%scaledEntries;
					
					if (particleCounter==0){
						fit->SetParameters(par0,par1,par2,0.55,0.01,2,2);
						fit->SetParLimits(3,0.525,0.575);
						fit->SetParLimits(4,0.005,0.02);
					}
					else{
						fit->SetParameters(par0,par1,par2,0.134,0.01,3,2);
						fit->SetParLimits(3,0.1275,0.1425);
						fit->SetParLimits(4,0.005,0.015);
					}
					fit->SetParLimits(0,0,totalEntries);
					fit->SetParLimits(2,0,totalEntries);
					fit->SetParLimits(5,0,3);
					fit->SetParLimits(6,0,3);

					Int_t fitStatus = anyHist->Fit(fit,"RLQ0");
					
					fit->GetParameters(par);
					bkgFit->SetParameters(par);
					sigFit->SetParameters(&par[numDOFbkg]);
					bkgFit->SetLineColor(kMagenta);
					sigFit->SetLineColor(kGreen+2);



					double w1 = sqrt(1/(1+par[3+numDOFbkg]));
					double w2 = sqrt(par[3+numDOFbkg]/(1+par[3+numDOFbkg]));
					weightedSigma = w1*par[2+numDOFbkg]+w2*par[4+numDOFbkg]*par[2+numDOFbkg];

					if(particleCounter==0){
						integralBKG = bkgFit->Integral(par[1+numDOFbkg]-nSig*weightedSigma,par[1+numDOFbkg]+nSig*weightedSigma);
						integralSIG = sigFit->Integral(par[1+numDOFbkg]-nSig*weightedSigma,par[1+numDOFbkg]+nSig*weightedSigma);
					}
					else {
						integralBKG = bkgFit->Integral(par[1+numDOFbkg]-nSig*weightedSigma,par[1+numDOFbkg]+nSig*weightedSigma);
						integralSIG = sigFit->Integral(par[1+numDOFbkg]-nSig*weightedSigma,par[1+numDOFbkg]+nSig*weightedSigma);
					}

					purity = integralSIG/(integralBKG+integralSIG);			

					if(fitStatus == 0 && fit->GetChisquare()/(fit->GetNDF()) < 30 && purity>=0 && purity<=1 ){
						break;
						cout << "------------- FOUND CONVERGENCE CRITERIA AFTER " << iAttempt << " ITERATIONS" << endl;
					}
					if (iAttempt==maxAttempts-1) {
						cout << "------------- DID NOT FIND CONVERGENCE CRITERIA AFTER " << maxAttempts << " ITERATIONS" << endl;
						integralBKG=0;
						integralSIG=0;
						purity=0;
						chiSq=0;
					}
				}



				cout << "weightedSigma, sigma1, sigma2: "  << weightedSigma << ", " << par[2+numDOFbkg] << ", " << par[2+numDOFbkg]*par[4+numDOFbkg] << endl;

				allCanvases_fits->cd(particleCounter+1);
				anyHist->Draw();
				cout << "Drew fit and histogram for particle " << particleCounter << " (" << particle << ")" << endl;

				chiSq=fit->GetChisquare()/(fit->GetNDF());

				// Saving all the integrals and chiSqs and purities
				integralBKGs[particleCounter][cumulativeBinNum] = integralBKG;
				integralSIGs[particleCounter][cumulativeBinNum] = integralSIG;
				puritys[particleCounter][cumulativeBinNum] = purity;
				chiSqs[particleCounter][cumulativeBinNum] = chiSq;

				// Saving the results into arrays so I can plot them 
				graph_integralBKGs[particleCounter][bin] = integralBKG;
				graph_integralSIGs[particleCounter][bin] = integralSIG;
				graph_puritys[particleCounter][bin] = purity;
				graph_chiSqs[particleCounter][bin] = chiSq;

				fit->Draw("SAME");
				bkgFit->Draw("SAME");
				sigFit->Draw("SAME");
				
				drawFitBounds->SetLineColor(kRed);
				cout << "Fit bounds are between ("<<par[1+numDOFbkg]-nSig*weightedSigma<<","<<par[1+numDOFbkg]+nSig*weightedSigma<<")"<<endl;
				cout << "  Drawn with an amplitude of: " << par[0+numDOFbkg]/par[2+numDOFbkg]/sqrt(2*TMath::Pi()) << endl;
				drawFitBounds->DrawLine(par[1+numDOFbkg]-nSig*weightedSigma,0,par[1+numDOFbkg]-nSig*weightedSigma,par[0+numDOFbkg]/par[2+numDOFbkg]/sqrt(2*TMath::Pi()));
				drawFitBounds->DrawLine(par[1+numDOFbkg]+nSig*weightedSigma,0,par[1+numDOFbkg]+nSig*weightedSigma,par[0+numDOFbkg]/par[2+numDOFbkg]/sqrt(2*TMath::Pi()));

				// Adding some text to output the results
				TPaveText *histID_pt = new TPaveText(0.7,0.7,0.9,0.9,"trNDC");
				histID_pt->AddText(("IntegralBKG="+std::to_string(integralBKG)).c_str());
				histID_pt->AddText(("IntegralSIG="+std::to_string(integralSIG)).c_str());
				histID_pt->AddText(("purity="+std::to_string(purity)).c_str());
				histID_pt->AddText(("chiSq="+std::to_string(chiSq)).c_str());
				histID_pt->Draw();
				++particleCounter;
			}
			allCanvases_fits->SaveAs(("newGraphs/sigEffs/fit_Mass_Mpi0etaChiBin"+to_string(cumulativeBinNum)+".png").c_str());
		}
		// We will fill up all the results for all the chiSq bins for a specific Mpi0eta bin
		gr_eta.push_back(new TGraph(numRegions_ChiSq,chiSqThresholds,graph_puritys[0]));
		gr_eta[bin_Mpi0eta]->SetMarkerStyle(20+bin_Mpi0eta);
		gr_eta[bin_Mpi0eta]->SetTitle(("Mpi0eta Bin: "+to_string(bin_Mpi0eta)).c_str());
		mg_eta->Add(gr_eta[bin_Mpi0eta],"PL");

		// include a legend on just on of the pads and use that as a reference
		std::stringstream streamEta;
		std::stringstream streamPi0;
		double currentMin=minMpi0eta+bin_Mpi0eta*binStepMpi0eta;
		double currentMax=minMpi0eta+(bin_Mpi0eta+1)*binStepMpi0eta;
		streamEta << std::fixed << std::setprecision(2) << currentMin;
		std::string sCurrentMin = streamEta.str();
		streamPi0 << std::fixed << std::setprecision(2) << currentMax;
		std::string sCurrentMax = streamPi0.str();
		legend->AddEntry(gr_eta[bin_Mpi0eta],(sCurrentMin+"<Mpi0eta<"+sCurrentMax).c_str(),"lp");

		gr_pi0.push_back(new TGraph(numRegions_ChiSq,chiSqThresholds,graph_puritys[1]));
		gr_pi0[bin_Mpi0eta]->SetMarkerStyle(20+bin_Mpi0eta);
		gr_pi0[bin_Mpi0eta]->SetTitle(("Mpi0eta Bin: "+to_string(bin_Mpi0eta)).c_str());
		mg_pi0->Add(gr_pi0[bin_Mpi0eta],"PL");
	}

	allCanvases_puritys->cd(1);
	mg_eta->Draw("A pmc plc");
	mg_eta->GetYaxis()->SetRangeUser(0,1.3);
	mg_eta->GetXaxis()->SetTitle("#Chi^{2} Threshold");
	legend->Draw();
	allCanvases_puritys->cd(2);
	mg_pi0->Draw("A pmc plc");
	mg_pi0->GetYaxis()->SetRangeUser(0,1.3);
	mg_pi0->GetXaxis()->SetTitle("#Chi^{2} Threshold");
	allCanvases_puritys->SaveAs("newGraphs/sigEffs/purities_Mpi0etaBin.png");

	//for ( auto particle : particles){
	//	logFile << "================== " << particle << " ==================" << endl;
	//	for (int bin_Mpi0eta=0; bin_Mpi0eta<numRegions_Mpi0eta; ++ bin_Mpi0eta){
	//		logFile << "------------------ Mpi0eta Bin: " << bin_Mpi0eta << " ------------------" << endl;
        //		for (int bin=0; bin<numRegions_ChiSq; ++bin){
	//			cumulativeBinNum = bin_Mpi0eta*numRegions_ChiSq+bin;
	//			logFile << "( ChiBin:" << bin << " ) IntegralBKG: " << integralBKGs[cumulativeBinNum] << "   IntegralSIG: " << integralSIGs[cumulativeBinNum] << "   Purity: " << puritys[cumulativeBinNum] << "   ChiSqPerDOF: " << chiSqs[cumulativeBinNum] << endl;
	//		}
	//		logFile << "--------------------------------------------------------------" << endl;
	//		++particleCounter;
	//	}
	//}
}



































