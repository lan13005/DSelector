

void plotDiagnostics(){
	double minMass=0.7;
	double maxMass=2.0;
	int fullNumBins=65;
	int numBins=65;
	int bins=65;
	double binWidth=(maxMass-minMass)/numBins;

	string fitDir = "EtaPi_fit";
	TFile* diagnosticFile;
	TCanvas* c1 = new TCanvas("","",1440,900);
	TH1F* dataHist;
	TH1F* genHist_neg;
	TH1F* genHist_pos;
	TH1F* accHist_neg;
	TH1F* accHist_pos;
	string diagnosticFileLoc;
	std::vector<string> histFileName={"cosTheta","phi","Phi","psi","t"};
	for (int iHist=0; iHist < histFileName.size(); ++iHist){
		for (int iBin=0; iBin<bins; ++iBin){
			THStack* stackedHist = new THStack("dataVsgen","");
			diagnosticFileLoc=fitDir+"/bin_"+to_string(iBin)+"/diagnostic.root";
			cout << "Loading: " << diagnosticFileLoc << endl;
		 	diagnosticFile = TFile::Open(diagnosticFileLoc.c_str());
			diagnosticFile->GetObject((histFileName[iHist]+"dat_Negative").c_str(),dataHist);
			diagnosticFile->GetObject((histFileName[iHist]+"gen_Negative").c_str(),genHist_neg);
			diagnosticFile->GetObject((histFileName[iHist]+"acc_Negative").c_str(),accHist_neg);
			diagnosticFile->GetObject((histFileName[iHist]+"gen_Positive").c_str(),genHist_pos);
			diagnosticFile->GetObject((histFileName[iHist]+"acc_Positive").c_str(),accHist_pos);
			dataHist->SetLineColor(kRed+1);
			stackedHist->Add(dataHist);
			genHist_neg->Add(genHist_pos);
			accHist_neg->Add(accHist_pos);
			genHist_neg->SetLineColorAlpha(kYellow+1,0.5);
			genHist_neg->SetFillColorAlpha(kYellow+1,0.5);
			accHist_neg->SetLineColorAlpha(kBlue+1,0.3);
			accHist_neg->SetFillColorAlpha(kBlue+1,0.3);
			stackedHist->Add(genHist_neg,"HIST");
			stackedHist->Add(accHist_neg,"HIST");
			stackedHist->SetTitle(("Mass = "+to_string(minMass+(0.5+iBin)*binWidth)+" GeV").c_str());
			stackedHist->Draw("nostack");
			//stackedHist->GetXaxis()->SetTitle(dataHist->GetXaxis()->GetTitle());
			//stackedHist->GetYaxis()->SetTitle(dataHist->GetYaxis()->GetTitle());
			c1->Update();
			c1->SaveAs(("diagnosticPlots/"+histFileName[iHist]+"_"+to_string(iBin)+".png").c_str());
			c1->Clear();
		}
	}	
}



