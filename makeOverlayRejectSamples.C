void makeOverlayRejectSamples(){
	TH1F* genHist;
	TH1F* accHist;
	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	allCanvases->Divide(2,1);
	gStyle->SetPalette(kOcean);

	string orientations[1] = {"000"};//,"045","090","135","AMO"};
	string dataTypes[2] = {"gen","acc"};
	string titles[2] = {"(gen) Thrown","(acc) Thrown"};

	int dTypeCounter=1;
	for ( auto dataType : dataTypes ){
		for ( auto orientation : orientations){
			auto hs = new THStack("hs", titles[dTypeCounter-1].c_str());
			allCanvases->cd(dTypeCounter);
			TFile* dataFile = TFile::Open(("deg000_"+dataType+"_hists_DSelector.root").c_str());
			dataFile->GetObject(("prodPlanePS_"+orientation).c_str(),genHist);	
			dataFile->GetObject(("prodPlanePS_"+orientation+"_rejSamp").c_str(),accHist);	
		
			genHist->SetFillColorAlpha(kRed+2,0.8);
			accHist->SetFillColorAlpha(kYellow,0.8);
			hs->Add(genHist);
			hs->Add(accHist);
			hs->Draw("NOSTACK");
			hs->GetXaxis()->SetTitle("#phi production plane");
			++dTypeCounter;
		}
	}
	allCanvases->SaveAs("rejectionSampling.pdf");
}
