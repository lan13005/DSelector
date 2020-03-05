int numDOFsig = 3;
Double_t shiftedCos(Double_t *x, Double_t *par){
	return par[0]*(1.0 + par[1]*TMath::Cos(2*(x[0]-par[2])));
}


void fitAsymmetryPlots(){
	gStyle->SetOptFit();
	gStyle->SetStatY(0.95);
	gStyle->SetStatX(0.95);
	gStyle->SetStatW(0.1);
	gStyle->SetStatH(0.1);
	TFile *dataFile = new TFile("deg000_data_hists_DSelector.root");

	TCanvas *allCanvases = new TCanvas("","",1440,900);
	allCanvases->Divide(3,2);

	TH1F *phi000;
	TH1F *phi045;
	TH1F *phi090;
	TH1F *phi135;
	TH1F *phiAMO;
	string names[5] = {"phi000","phi045","phi090","phi135", "phiAMO"};
	
	dataFile->GetObject("prodPlanePSphi_000_cut", phi000);
	dataFile->GetObject("prodPlanePSphi_045_cut", phi045);
	dataFile->GetObject("prodPlanePSphi_090_cut", phi090);
	dataFile->GetObject("prodPlanePSphi_135_cut", phi135);
	dataFile->GetObject("prodPlanePSphi_AMO_cut", phiAMO);

	TH1F *phis[5] = {phi000, phi045, phi090, phi135, phiAMO};

	TF1 * fit = new TF1("fit",shiftedCos,-180,180,numDOFsig); 
	double p0;
	int counter=0;
	for (auto phi: phis){
		allCanvases->cd(counter+1);
		phi->SetTitle((names[counter]).c_str());
		p0 = phi->GetEntries()/phi->GetNbinsX();
		fit->SetParameters(p0,p0/2,0);
		Int_t fitStatus = phi->Fit(fit,"RLQE");
		phi->Draw("SAME");
		++counter;
	}
	allCanvases->SaveAs("asymmetryPlots/asymmetry.png");

}
