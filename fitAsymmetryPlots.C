double degToRad=TMath::Pi()/180;
int numDOFsig_sc = 3;
Double_t shiftedCos(Double_t *x, Double_t *par){
	return par[0]*(1.0 + par[1]*TMath::Cos(2*degToRad*(x[0]-par[2])));
}

int numDOFsig_flat = 1;
Double_t flat(Double_t *x, Double_t *par){
	return par[0];
}

int numDOFsig_asym=4;
Double_t asymmetry(Double_t *x, Double_t *par){
	return ((par[0]+par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3])))/(2+(par[0]-par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3])));
}

void fitAsymmetryPlots(){
	static const int num_tBins=5;

	gStyle->SetOptFit();
	gStyle->SetStatY(0.95);
	gStyle->SetStatX(0.95);
	gStyle->SetStatW(0.1);
	gStyle->SetStatH(0.1);
	TFile *dataFile = new TFile("deg000_data_hists_DSelector.root");

	TCanvas *allCanvases = new TCanvas("","",1440,900);

	TH1F *phi000;
	TH1F *phi045;
	TH1F *phi090;
	TH1F *phi135;
	TH1F *phiAMO;
	string names[4] = {"phi000","phi045","phi090","phi135"};
	double asymmetries_000[num_tBins];
	double asymmetries_045[num_tBins];
	double asymmetries_000_err[num_tBins];
	double asymmetries_045_err[num_tBins];
	double tBins[num_tBins];
	double tBins_err[num_tBins];
	double tBinSize=0.2;
	
	for (int iteta=0; iteta<num_tBins; ++iteta){
		dataFile->GetObject(("prodPlanePSphi_000_tetaBin"+to_string(iteta)).c_str(), phi000);
		dataFile->GetObject(("prodPlanePSphi_045_tetaBin"+to_string(iteta)).c_str(), phi045);
		dataFile->GetObject(("prodPlanePSphi_090_tetaBin"+to_string(iteta)).c_str(), phi090);
		dataFile->GetObject(("prodPlanePSphi_135_tetaBin"+to_string(iteta)).c_str(), phi135);
		dataFile->GetObject(("prodPlanePSphi_AMO_tetaBin"+to_string(iteta)).c_str(), phiAMO);

		TH1F *phis[4] =  {phi000, phi045, phi090, phi135};
		double orientation[4] = {0,45,90,135};

		TH1* asymmetry000_090 = phi000->GetAsymmetry(phi090);
		TH1* asymmetry045_135 = phi045->GetAsymmetry(phi135);
		asymmetry000_090->SetTitle("0/90 Asymmetry");
		asymmetry045_135->SetTitle("45/135 Asymmetry");

		// We set P_perp = P_para = 0.4 which is close to the expected. We initialize asymmetry to be 0 since it can vary from [-1,1]. 
		TF1 * fit_asym = new TF1("fit_asym",asymmetry,-180,180,numDOFsig_asym); 
		fit_asym->SetParameters(0.4,0.4,0.5,0);
		fit_asym->FixParameter(0,0.4);
		fit_asym->FixParameter(1,0.4);
		fit_asym->FixParameter(3,0);
		//fit_asym->SetParLimits(3,0,90); // want to constrain phase to be between -pi/2 and pi/2 so that with the factor of 2 it will not allow for arbitrary (-) signs
		allCanvases->Clear();
		allCanvases->Divide(2,1);
		allCanvases->cd(1);
		Int_t fitStatus = asymmetry000_090->Fit(fit_asym,"RLQE");
		asymmetries_000[iteta] = fit_asym->GetParameter(2);
		asymmetries_000_err[iteta] = fit_asym->GetParError(2);
		asymmetry000_090->Draw("SAME");
		asymmetry000_090->SetAxisRange(-0.5,0.6,"Y");

		allCanvases->cd(2);
		fit_asym->SetParameters(0.4,0.4,0.5,45);
		fit_asym->FixParameter(0,0.4);
		fit_asym->FixParameter(1,0.4);
		fit_asym->FixParameter(3,45);
		fitStatus = asymmetry045_135->Fit(fit_asym,"RLQE");
		asymmetry045_135->Draw("SAME");
		asymmetry045_135->SetAxisRange(-0.5,0.6,"Y");
		asymmetries_045[iteta] = fit_asym->GetParameter(2);
		asymmetries_045_err[iteta] = fit_asym->GetParError(2);
		tBins[iteta] = iteta*tBinSize;
		tBins_err[iteta] = tBinSize/2;
		allCanvases->SaveAs(("asymmetryPlots/asymmetry_tetaBin"+to_string(iteta)+".png").c_str());


		allCanvases->Clear();
		allCanvases->Divide(3,2);
		TF1 * fit_sc = new TF1("fit_sc",shiftedCos,-180,180,numDOFsig_sc); 
		double p0;
		int counter=0;
		for (auto phi: phis){
			allCanvases->cd(counter+1);
			phi->SetTitle((names[counter]).c_str());
			p0 = phi->GetEntries()/phi->GetNbinsX();
			fit_sc->SetParameters(p0,p0/2,orientation[counter]);
			fitStatus = phi->Fit(fit_sc,"RLQE");
			phi->Draw("SAME");
			++counter;
		}
		TF1 * fit_flat = new TF1("fit_flat",flat,-180,180,numDOFsig_flat); 
		allCanvases->cd(5);
		phiAMO->SetTitle("phiAMO");
		p0 = phiAMO->GetEntries()/phiAMO->GetNbinsX();
		fit_sc->SetParameters(p0,p0/2,orientation[counter]);
		fitStatus = phiAMO->Fit(fit_flat,"RLQE");
		phiAMO->Draw("SAME");

		allCanvases->SaveAs(("asymmetryPlots/phiYieldFits_tetaBin"+to_string(iteta)+".png").c_str());
	}

	allCanvases->Clear();
	allCanvases->Divide(2,1);
	allCanvases->cd(1);
	auto gr_000 = new TGraphErrors(num_tBins,tBins,asymmetries_000,tBins_err,asymmetries_000_err);
	gr_000->SetTitle("Beam Asymmetry 0/90 Orientation");
	gr_000->SetMarkerColor(4);
	gr_000->SetMarkerStyle(21);
	gr_000->Draw("AP");
	gr_000->GetHistogram()->SetMaximum(1.2);
	gr_000->GetHistogram()->SetMinimum(-1);

	allCanvases->cd(2);
	auto gr_045 = new TGraphErrors(num_tBins,tBins,asymmetries_045,tBins_err,asymmetries_045_err);
	gr_045->SetTitle("Beam Asymmetry 45/135 Orientation");
	gr_045->SetMarkerColor(4);
	gr_045->SetMarkerStyle(21);
	gr_045->Draw("AP");
	gr_045->GetHistogram()->SetMaximum(1.2);
	gr_045->GetHistogram()->SetMinimum(-1);
	allCanvases->SaveAs("asymmetryPlots/asymVsteta.png");
}
















