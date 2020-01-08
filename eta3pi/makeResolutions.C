void makeResolutions(){
	TCanvas* canvas=new TCanvas("","",1440,900);
	canvas->Divide(2,1);
	gStyle->SetOptStat(0);
	TFile* geneFile=new TFile("eta3pi_a2_gen_hists_DSelector_pi0eta.root");
	TFile* recoFile=new TFile("pi0eta_eta3pi0_hists_DSelector.root");
	TH1F *recoHist;
	TH1F *geneHist;
	recoFile->GetObject("pi0eta1D_mMandelstamT_mBeamE8GeVPlus", recoHist);
	geneFile->GetObject("pi0eta1D", geneHist);

	Double_t recoNorm = recoHist->GetEntries();
	cout << "reco has " << recoNorm << " entries" << endl;
	recoNorm = 1;
	recoHist->Scale(recoNorm/recoHist->Integral(),"width");
	recoHist->Scale(1/recoNorm);

	Double_t geneNorm = geneHist->GetEntries();
	cout << "gen has " << geneNorm << " entries" << endl;
	geneNorm = 1;
	geneHist->Scale(geneNorm/geneHist->Integral(),"width");
	
	canvas->cd(1);
	geneHist->SetLineColor(kRed);
	geneHist->Draw("HIST");
	geneHist->SetTitle("Reconstructed Vs Thrown Normalized");
	recoHist->Draw("HIST SAME");

	auto legend = new TLegend(0.6,0.7,0.9,0.9);
	legend->AddEntry(recoHist,"reconstructed","l");
	legend->AddEntry(geneHist,"thrown","l");
	legend->Draw();

	canvas->Update();

	canvas->cd(2);
	TH1F* resolution = (TH1F*)recoHist->Clone();
	resolution->Add(geneHist,-1);
	resolution->Draw("HIST");
	resolution->SetTitle("Resolution (Reconstructed-Thrown)");


	canvas->SaveAs("resolutionComparison.png");
}
