int numDOFsig=8;
Double_t g2(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	Double_t r2 = Double_t((x[1]-par[3])/par[4]);
	return par[5]+par[6]*(x[0]-par[1])+par[7]*(x[1]-par[3])+par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}

int numDOFsig1D=5;
Double_t gaus(Double_t *x, Double_t *par){
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	return par[3]+par[4]*(x[0]-par[1])+par[0]/sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-0.5*r1*r1 );
}

void makeDeckPlots(){
	TFile* dataFile = TFile::Open("pi0eta_datatreeFlat_DSelector.root");
	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	TTree *dataTree;
	dataFile->GetObject("pi0eta_datatree_flat",dataTree);

	double new_t;
	double new_m;
	double new_accWeight;
	double new_meta;
	double new_mpi0;
        dataTree->SetBranchAddress("AccWeight",&new_accWeight);
	dataTree->SetBranchAddress("Meta_meas",&new_meta);
	dataTree->SetBranchAddress("Mpi0_meas",&new_mpi0);
	dataTree->SetBranchAddress("Mpi0eta_meas",&new_m);
	dataTree->SetBranchAddress("mandelstam_teta_meas",&new_t);

	int num_tBins=14;
	double tMin=0;
	double tMax=2.8;
	int num_massBins=12;
	double mMin=1.7;
	double mMax=2.9;
	double tStep=(tMax-tMin)/num_tBins;
	double mStep=(mMax-mMin)/num_massBins;

	// So the following indicies would describe the bin for a specific variable. We have to multiplex them into a single array. 
	// The obvious way to do that is to use the array index = 3*m+t
	int idx_t;
	int idx_m;
	const int numHists = (const int)num_tBins*num_massBins;
	cout << "numHists: " << numHists << endl;
	TH2F* hists_meta_mpi0[numHists];
	TH2F* hists_mpi0eta_t[numHists];
	TH2F* full_meta_mpi0 = new TH2F("full_meta_mpi0", "Cuts=GeneralCuts;M(#pi_{0}) GeV;M(#eta)(GeV)", 100,0.05,0.25,100,0.25,0.85 );

	double xfitMin=0.11;
	double xfitMax=0.165;
	double yfitMin=0.46;
	double yfitMax=0.62;
	TH1F* full_meta = new TH1F("full_meta", "Cuts=GeneralCuts;M(#eta)(GeV)",100,0.25,0.85 );
	TH1F* full_mpi0 = new TH1F("full_mpi0", "Cuts=GeneralCuts;M(#pi_0)(GeV)",100,0.05,0.25 );
	for (int i=0; i< numHists; ++i){
		hists_meta_mpi0[i] = new TH2F(("meta_mpi0"+to_string(i)).c_str(), "Cuts=GeneralCuts;M(#pi_{0}) GeV;M(#eta)(GeV)", 100,0.05,0.25,100,0.25,0.85 ); 
		hists_mpi0eta_t[i] = new TH2F(("mpi0eta_t"+to_string(i)).c_str(), "Cuts=mMandelstamT_eta;M(#pi_{0}#eta) (GeV);t_{#eta} (GeV^2)", 260, 0.6, 3.2, 80,0,8);
	}
	int histIdx;

	Long64_t nentries = dataTree->GetEntries();

	for (int ientry=0; ientry<nentries; ++ientry){
		dataTree->GetEntry(ientry);
		if ( new_t < tMin || new_m < mMin || new_t>tMax || new_m>mMax) {continue; }
		idx_t = (int)( (new_t-tMin)/tStep ); 
		idx_m = (int)( (new_m-mMin)/mStep );
		//cout << idx_t << ", " << idx_m << endl;
		histIdx=num_tBins*idx_m+idx_t;
		hists_meta_mpi0[histIdx]->Fill(new_mpi0,new_meta,new_accWeight);
		hists_mpi0eta_t[histIdx]->Fill(new_m,new_t,new_accWeight);
		full_meta_mpi0->Fill(new_mpi0,new_meta,new_accWeight);
		full_meta->Fill(new_meta,new_accWeight);
		full_mpi0->Fill(new_mpi0,new_accWeight);
	}

	int fit_nentries=full_meta_mpi0->GetEntries();
	cout << "TOTAL ENTRIES: " << fit_nentries << endl;
	gStyle->SetOptFit();
	//for(int i=0; i<numHists; ++i){
	TF2 * f2 = new TF2("f2",g2,xfitMin,xfitMax,yfitMin,yfitMax,numDOFsig); 
	double initParams[8]={18,0.1348, 0.0952, 0.5393, 0.0215, 9, 100, 23}
	f2->SetParameters(initParams[0],initParams[1],initParams[2],initParams[3],initParams[4],initParams[5],initParams[6],initParams[7]); 
	f2->SetParLimits(0,0,500);
	f2->SetParLimits(1,0.13,0.14);
	f2->SetParLimits(3,0.52,0.585);
	f2->SetParLimits(2,0.01,0.03);
	f2->SetParLimits(4,0,0.1);
	//f2->SetParLimits(5,-full_meta_mpi0->GetEntries()/10000,full_meta_mpi0->GetEntries()/10000);
	//f2->SetParLimits(6,-full_meta_mpi0->GetEntries()/10000/0.02/100, full_meta_mpi0->GetEntries()/10000/0.02/100);
	//f2->SetParLimits(7,-full_meta_mpi0->GetEntries()/10000/0.09/100, full_meta_mpi0->GetEntries()/10000/0.09/100);
	allCanvases->Clear();
	full_meta_mpi0->Fit(f2,"RL");
	full_meta_mpi0->Draw("COLZ");
	//full_meta_mpi0->Draw("SURF2");
	//f2->Draw("SURF SAME");
	allCanvases->SaveAs("deckPlots/full_meta_mpi0.png");	
	allCanvases->Clear();
	TF1 * feta = new TF1("meta_fit",gaus,yfitMin,yfitMax,numDOFsig1D); 
	full_meta->Draw();
	feta->SetParameters(2000,0.55,0.03,0, 0);
	feta->SetParLimits(0,0,5000);
	feta->SetParLimits(1,0.52,0.585);
	feta->SetParLimits(2,0,0.1);
	full_meta->Fit(feta,"RL");
	allCanvases->SaveAs("deckPlots/project_full_meta.png");
	allCanvases->Clear();
	TF1 * fpi0 = new TF1("mpi0_fit",gaus,xfitMin,xfitMax,numDOFsig1D); 
	fpi0->SetParameters(500,0.135,0.015,0, 0);
	fpi0->SetParLimits(0,0,1000);
	fpi0->SetParLimits(1,0.125,0.145);
	fpi0->SetParLimits(2,0.005,0.02);
	f2->SetParLimits(0,0,500);
	f2->SetParLimits(1,0.13,0.14);
	f2->SetParLimits(3,0.52,0.585);
	f2->SetParLimits(2,0.01,0.03);
	f2->SetParLimits(4,0,0.1);
	full_mpi0->Fit(fpi0,"RL");
	allCanvases->SaveAs("deckPlots/project_full_mpi0.png");


	for(int i=0; i<10; ++i){
		f2->SetParameters(2000 , 0.1355, 0.02, 0.55, 0.08, 0, 0, 0); 
		f2->SetParLimits(0,0,hists_meta_mpi0[i]->GetEntries());
		f2->SetParLimits(7,0,hists_meta_mpi0[i]->GetEntries());
		f2->SetParLimits(1,0.13,0.14);
		f2->SetParLimits(3,0.54,0.57);
		hists_meta_mpi0[i]->Fit(f2,"RL"); 

		allCanvases->Clear();
		hists_meta_mpi0[i]->Draw("COLZ");
		allCanvases->SaveAs(("deckPlots/meta_mpi0-"+to_string(i)+".png").c_str());	
		allCanvases->Clear();
		hists_mpi0eta_t[i]->Draw("COLZ");
		allCanvases->SaveAs(("deckPlots/mpi0eta_t-"+to_string(i)+".png").c_str());	
	}


	//double currLow_m;
	//double currUp_m;
	//double currLow_t;
	//double currUp_t;
	//for (int iMass=0; iMass<num_massBins; iMass++){
	//	currLow_m = mMin+iMass*mStep;
	//	currUp_m = mMin+(iMass+1)*mStep;
	//	//cout << " #### " << currLow_m << " < M(pi0eta) < " << currUp_m << " #### " << endl;
	//	for (int iT=0; iT<num_tBins; iT++){
	//		currLow_t = tMin+iT*tStep;
	//		currUp_t = tMin+(iT+1)*tStep;
	//		//cout << currLow_t << " < t < " << currUp_t << endl;
	//	}
	//}


	//allCanvases->SaveAs(("acceptance/"+names2D_gen[0]+".png").c_str());

}
