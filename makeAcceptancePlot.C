void makeAcceptancePlot(){
	TFile* infile_acc = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_flat_2017_mEllipsePre_hists_DSelector.root");
	TFile* infile_gen = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_gen_2017_hists_DSelector.root");
	TFile* infile_dat = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_data_2017_mEllipsePre_hists_DSelector.root");
	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);

	std::vector<string> names2D_gen = {"cosThetaVsMass_tpAll"};
	std::vector<string> names2D_acc = {"eta_cosTheta_GJvsM_baseAsymCut"};

	TH2F *any2DHist_acc;
	TH2F *any2DHist_gen;
	TH2F *any2DHist_dat;

	TH1F *any1DHist_acc;
	TH1F *any1DHist_gen;
	TH1F *any1DHist_dat;

	for (int i=0; i<names2D_gen.size(); ++i){
		cout << endl;
		infile_acc->GetObject(names2D_acc[i].c_str(),any2DHist_acc);
		infile_dat->GetObject(names2D_acc[i].c_str(),any2DHist_dat);
		infile_gen->GetObject(names2D_gen[i].c_str(),any2DHist_gen);

		cout << "MAKING " << names2D_gen[i] << " WITH " << names2D_acc[i] <<  " RIGHT NOW!" << endl;
		TH2F *acceptance = new TH2F(*any2DHist_acc);
		acceptance->SetStats(0);
		acceptance->Divide(any2DHist_gen);
		acceptance->Draw("COLZ");
		cout << "Saving acceptance plot to: " << "acceptance/acceptance-"+names2D_acc[i]+".png" << endl;
		allCanvases->SaveAs(("acceptance/acceptance-"+names2D_acc[i]+".png").c_str());

		TH2F *any2DHist_dat_corrected = new TH2F(*any2DHist_dat);
		any2DHist_dat_corrected->Divide(acceptance);
		allCanvases->Clear();
		any2DHist_dat_corrected->Draw("COLZ");
		cout << "Saving corrected plot to: " << "acceptance/corrected-"+names2D_acc[i]+".png" << endl;
		allCanvases->SaveAs(("acceptance/corrected-"+names2D_acc[i]+".png").c_str());

		allCanvases->Clear();
		any2DHist_dat->Draw("COLZ");
		cout << "Saving corrected plot to: " << "acceptance/uncorrected-"+names2D_acc[i]+".png" << endl;
		allCanvases->SaveAs(("acceptance/uncorrected-"+names2D_acc[i]+".png").c_str());

		allCanvases->Clear();
		infile_acc->GetObject("mandelstam_t",any1DHist_acc);	
		infile_dat->GetObject("mandelstam_t",any1DHist_dat);	
		any1DHist_acc->SetLineColor(kRed);
		any1DHist_dat->Draw("HIST");
		any1DHist_acc->Draw("HIST SAME");	
		allCanvases->SaveAs("acceptance/t-slope_datVsMC.png");

		allCanvases->Clear();
		infile_acc->GetObject("mandelstam_tp",any1DHist_acc);	
		infile_dat->GetObject("mandelstam_tp",any1DHist_dat);	
		any1DHist_acc->SetLineColor(kRed);
		any1DHist_dat->Draw("HIST");
		any1DHist_acc->Draw("HIST SAME");	
		allCanvases->SaveAs("acceptance/tp-slope_datVsMC.png");
	}
}






















