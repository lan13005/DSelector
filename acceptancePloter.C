void getAcceptance(){
	TFile* infile_acc = TFile::Open("pi0eta_flat_hists_DSelector.root");
	TFile* infile_gen = TFile::Open("v20_flat_gen_hists_DSelector_pi0eta.root");
	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);

	std::vector<string> names2D_gen = {"cosThetaVsMass_tpLT1"};
	std::vector<string> names2D_acc = {"eta_cosTheta_GJvsM_Cut"};

	TH2F *any2DHist_acc;
	TH2F *any2DHist_gen;
	infile_acc->GetObject(names2D_acc[0].c_str(),any2DHist_acc);
	infile_gen->GetObject(names2D_gen[0].c_str(),any2DHist_gen);

	cout << "MAKING " << names2D_gen[0] << " WITH " << names2D_acc[0] <<  " RIGHT NOW!" << endl;
	any2DHist_acc->SetStats(0);
	any2DHist_acc->Divide(any2DHist_gen);
	any2DHist_acc->Draw("COLZ");
	cout << "Saving to: " << "acceptance/"+names2D_gen[0]+".png" << endl;
	allCanvases->SaveAs(("acceptance/"+names2D_gen[0]+".png").c_str());

}
