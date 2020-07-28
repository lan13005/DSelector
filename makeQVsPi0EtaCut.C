void makeQVsPi0EtaCut(){
	gStyle->SetOptStat(0);
	TFile* dataFile = TFile::Open("degALL_data_2017_newProspectusGraphs/degALL_data_2017_hists_DSelector.root");
	TFile* qValueFile = TFile::Open("/d/grid15/ln16/pi0eta/q-values/diagnosticPlots/postQVal.root");
	TCanvas *allCanvases = new TCanvas("","",1440,900);
	TH1F *signalRegion_mEllipse;
	TH1F *signalRegion;
	TH1F *pi0eta_Qsubbed;

	dataFile->GetObject("pi0gamma_Cut",signalRegion);
	dataFile->GetObject("pi0gamma_mEllipse_pre",signalRegion_mEllipse);
	qValueFile->GetObject("Mpi0g_sig",pi0eta_Qsubbed);

	signalRegion_mEllipse->SetTitle("");
	
	TLegend* leg1 = new TLegend(0.65,0.7,0.9,0.9);
	leg1->AddEntry(pi0eta_Qsubbed,"Q_{AS}","l");
	leg1->AddEntry(signalRegion,"Signal ( Cut M(#pi^{0}) vs M(#eta) )","l");
	leg1->AddEntry(signalRegion_mEllipse,"No Cut on M(#pi^{0}) vs M(#eta)","l");

	pi0eta_Qsubbed->SetLineColor(kRed);
	signalRegion_mEllipse->SetLineColor(kBlue);
	signalRegion->SetLineColor(kGreen);
	signalRegion_mEllipse->Draw();
	signalRegion->Draw("SAME");	
	pi0eta_Qsubbed->SetStats(0); // for some reason the stats wont get removed with gStyle
	pi0eta_Qsubbed->Draw("SAME");

	TH1F *rectSB_4;
	TH1F *rectSB_0268;
	TH1F *rectSB_17;
	TH1F *rectSB_35;
	dataFile->GetObject("pi0gamma_RectSBSubRegion17",rectSB_17);
	dataFile->GetObject("pi0gamma_RectSBSubRegion35",rectSB_35);
	dataFile->GetObject("pi0gamma_RectSBSubRegion0268",rectSB_0268);
	dataFile->GetObject("pi0gamma_RectSBSubRegion4",rectSB_4);
	TH1F *rectSB_total = (TH1F*)rectSB_4->Clone();
	rectSB_total->SetName("rectSB_total");
	TH1F *rectSB_background = (TH1F*)rectSB_0268->Clone();
	rectSB_background->SetName("rectSB_background");
	rectSB_background->Scale(0.5);
	rectSB_background->Add(rectSB_17,-0.5);
	rectSB_background->Add(rectSB_35,-1);
	TH1F *rectSB_signal = (TH1F*)rectSB_total->Clone();
	rectSB_signal->SetName("rectSB_signal");
	rectSB_signal->Add(rectSB_background,1);

	rectSB_signal->SetLineColor(kMagenta);
	rectSB_signal->Draw("SAME");
	leg1->AddEntry(rectSB_signal,"Sideband Subtracted","l");

	leg1->Draw("SAME");
	allCanvases->SaveAs("newGraphs/pi0gamma_QVsSignal.png");
}
