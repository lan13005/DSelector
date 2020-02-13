void makeCompare2017_2018(){
	gStyle->SetOptStat(0);
	TFile* data_2018_8 = TFile::Open("pi0eta_all_tLT1_2018_8_hists_DSelector.root");
	TFile* data_2018_1 = TFile::Open("pi0eta_all_tLT1_2018_1_hists_DSelector.root");
	TFile* data_2017 = TFile::Open("pi0eta_all_tLT1_hists_DSelector.root");

	TH1F* hist_2018_8;
	TH1F* hist_2018_1;
	TH1F* hist_2017;
	TH1F* hist_2017_tCut;

	TLegend* leg1 = new TLegend(0.7,0.6,0.9,0.9);

	data_2018_8->GetObject("pi0eta1D_mMandelstamT",hist_2018_8);
	data_2018_1->GetObject("pi0eta1D_mMandelstamT",hist_2018_1);
	data_2017->GetObject("pi0eta1D_mMandelstamT",hist_2017);
	data_2017->GetObject("pi0eta1D_Cut",hist_2017_tCut);

	hist_2018_8->SetLineColor(kRed);
	hist_2018_1->SetLineColor(kGreen+1);
	hist_2017_tCut->SetLineColor(kMagenta);
	TCanvas* allCanvases = new TCanvas("","",1440,900);

	leg1->AddEntry(hist_2018_1,"2018_1 Data");
	leg1->AddEntry(hist_2018_8,"2018_8 Data");
	leg1->AddEntry(hist_2017,"2017 Data");
	leg1->AddEntry(hist_2017_tCut, "2017 Data - t'< 1 GeV^2");


	hist_2018_1->Draw();
	hist_2018_8->Draw("SAME");
	hist_2017->Draw("SAME");
	hist_2017_tCut->Draw("SAME");

	leg1->Draw();

	//TText *t = new TText();
	//t->SetTextFont(43);
	//t->SetTextSize(40);
	//t->SetTextColor(kGreen);
	//t->DrawTextNDC(0.6,0.5,("# Events: "+to_string(hist_2018_1->GetEntries())).c_str());
	//t->SetTextColor(kRed);
	//t->DrawTextNDC(0.6,0.45,("# Events: "+to_string(hist_2018_8->GetEntries())).c_str());
	//t->SetTextColor(kBlue);
	//t->DrawTextNDC(0.6,0.40,("# Events: "+to_string(hist_2017->GetEntries())).c_str());

	allCanvases->SaveAs("newGraphs/Mpi0eta_2017_2018.png");
}
