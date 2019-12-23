void makeResolutions(){
	TCanvas* canvas=new TCanvas("","",1440,900);
	canvas->Divide(2,1);
	TFile* geneFile=new TFile("eta3pi_a2_gen_hists_DSelector_pi0eta.root");
	TFile* recoFile=new TFile("pi0eta_eta3pi0_hists_DSelector.root");
	TH1F *recoHist;
	TH1F *geneHist;
	recoFile->GetObject("pi0eta1D_mMandelstamT_mBeamE8GeVPlus", recoHist);
	geneFile->GetObject("pi0eta1D", geneHist);

	Double_t recoNorm =recoHist->GetEntries();
	recoHist->Scale(1/recoNorm);
	Double_t geneNorm =geneHist->GetEntries();
	geneHist->Scale(1/geneNorm);
	
	canvas->cd(1);
	recoHist->Draw();
	recoHist->SetTitle("Reconstructed Vs Thrown");
	geneHist->SetLineColor(kRed);
	geneHist->Draw("SAME");

	auto legend = new TLegend(0.1,0.7,0.48,0.9);
	legend->AddEntry(recoHist,"reconstructed","l");
	legend->AddEntry(geneHist,"thrown","l");
	legend->Draw();

	canvas->cd(2);
	recoHist->Add(geneHist,-1);
	recoHist->Draw();
	recoHist->SetTitle("Reconstructed-Thrown Normalized");


	canvas->SaveAs("resolutionComparison.png");
}
