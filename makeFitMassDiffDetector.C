int numDOFbkg = 2;
Double_t background(Double_t *x, Double_t *par){
	return par[0]+par[1]*x[0];
}

int numDOFsig = 3;
Double_t signal(Double_t *x, Double_t *par){
	//return par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));
	return par[0]/par[2]/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]));// + par[3]*exp(-0.5*((x[0]-par[1])/par[4])*((x[0]-par[1])/par[4]));
}

Double_t fitFunc(Double_t *x, Double_t *par){
	return background(x,par)+signal(x,&par[numDOFbkg]);
}

void makeFitMassDiffDetector(){
	TFile* dataFile = TFile::Open("pi0eta_data_hists_DSelector.root");
	TCanvas *allCanvases = new TCanvas("","",1440,900);
	auto legend = new TLegend(0.1,0.7,0.48,0.9);
	TH1F *bcal_hist;
	TH1F *fcal_hist;
	TH1F *split_hist;
	TH1F *total_hist;

	dataFile->GetObject("etaMass_Kin_mEllipsePre",total_hist);
	dataFile->GetObject("etaMassFCAL_Kin_mEllipsePre",fcal_hist);
	dataFile->GetObject("etaMassBCAL_Kin_mEllipsePre",bcal_hist);
	dataFile->GetObject("etaMassSPLIT_Kin_mEllipsePre",split_hist);

	total_hist->Draw();

	TF1 *etaSPLITFit = new TF1("etaSPLITFit",fitFunc,0.4,0.7,numDOFsig+numDOFbkg);
	etaSPLITFit->SetParameters(-700,2800,225,0.540393, 0.018);
	//etaSPLITFit->FixParameter(0,-700);
	//etaSPLITFit->FixParameter(1,3100);
	//etaSPLITFit->FixParameter(2,225);
	//etaSPLITFit->FixParameter(3,0.540393);
	//etaSPLITFit->FixParameter(4,0.018);
	//etaSPLITFit->SetParLimits(3,0.53,0.55);
	//etaSPLITFit->SetParLimits(4,0.02,0.03);
	split_hist->Fit("etaSPLITFit");
	split_hist->Draw("SAME");
	//etaSPLITFit->Draw("SAME");
	//etaSPLITFit->SetLineColor(kMagenta+1);
	split_hist->SetLineColor(kMagenta+1);

	TF1 *etaFCALFit = new TF1("etaFCALFit",fitFunc,0.4,0.7,numDOFsig+numDOFbkg);
	etaFCALFit->SetParameters(700,0,5000,0.540393, 0.0233083);
	fcal_hist->Fit("etaFCALFit","N");
	fcal_hist->Draw("SAME");
	etaFCALFit->SetLineColor(kRed+1);
	etaFCALFit->Draw("SAME");
	fcal_hist->SetLineColor(kRed+1);

	TF1 *etaBCALFit = new TF1("etaBCALFit",fitFunc,0.4,0.7,numDOFsig+numDOFbkg);
	etaBCALFit->SetParameters(200,0,1000,0.540393, 0.0233083);
	bcal_hist->Fit("etaBCALFit","N");
	bcal_hist->Draw("SAME");
	etaBCALFit->SetLineColor(kGreen+1);
	etaBCALFit->Draw("SAME");
	bcal_hist->SetLineColor(kGreen+1);

	legend->AddEntry(total_hist,"Total");
	legend->AddEntry(fcal_hist,"FCAL");
	legend->AddEntry(bcal_hist,"BCAL");
	legend->AddEntry(split_hist,"SPLIT");

	allCanvases->SaveAs("gausFitEtaMassDiffDetector.png");

	
		
	
}
