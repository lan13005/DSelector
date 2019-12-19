#include "/d/grid15/ln16/pi0eta/q-values/main.h"

double calculateEfficiency( double var1s[], Long64_t nentries, double chiSqs[], double unusedEnergies[], int iter ) {
	TH1F *hist_var;
	TCanvas *canvas = new TCanvas("","",1440,900);
        double binRange[3] = {100,0.35,0.8};
        double fitRange[2] = {0.4,0.65};
	double par[numDOFbkg+numDOFsig];
	double parBkg[numDOFbkg];

	hist_var = new TH1F("","",binRange[0],binRange[1],binRange[2]);
	TF1* fit = new TF1("fit",fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
	TF1* bkgFit = new TF1("bkgFit",fitFunc,fitRange[0],fitRange[1],numDOFbkg);
	TF1* sigFit = new TF1("sigFit",fitFunc,fitRange[0],fitRange[1],numDOFsig);
	fit->SetParameters(4600,8600,20000,0.5406,0.0237);

	for ( auto ivar=0; ivar<nentries; ++ivar ){
		hist_var->Fill(var1s[ivar]);
	}

	hist_var->Fit("fit","RQB");
	fit->GetParameters(par);
	bkgFit->SetParameters(par);
	sigFit->SetParameters(&par[numDOFbkg]);

	double nSig = sigFit->Integral(binRange[0],binRange[1]);
	double nBkg = bkgFit->Integral(binRange[0],binRange[1]);

	hist_var->Draw();
	fit->SetLineColor(kRed);
	bkgFit->SetLineColor(kMagenta);
	bkgFit->SetFillColor(kMagenta);
	bkgFit->SetFillStyle(3344);
	sigFit->SetLineColor(kBlue);
	sigFit->SetFillColor(kBlue);
	sigFit->SetFillStyle(3344);
	fit->Draw("same");
	bkgFit->Draw("same");
	sigFit->Draw("same");

	canvas->SaveAs(("/d/grid15/ln16/pi0eta/chiSqUEplots/fit"+to_string(iter)+".png").c_str());

	return nSig/nBkg; 
}


void calcSigEfficiency(){
	cout << "Initializing" << endl;
	TFile* dataFile=new TFile("/d/grid15/ln16/pi0eta/092419/zSelectedmEllipseUEChiSq/pi0eta_mEllipseUEChiSq_pretreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_mEllipseUEChiSq_pretree_flat",dataTree);

	Long64_t nentries=dataTree->GetEntries();
	cout << "There are " << nentries << " nentries" << endl;
	double Metas[nentries];
	double chiSqs[nentries];
	double unusedEnergies[nentries];
	double Meta;
	double chiSq;
	double unusedEnergy;
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("chiSq",&chiSq);
	dataTree->SetBranchAddress("unusedEnergy",&unusedEnergy);

	cout << "Loading the data" << endl;
	for (Long64_t ientry=0; ientry<nentries; ientry++)
	{
        	dataTree->GetEntry(ientry);
		Metas[ientry] = Meta;
		unusedEnergies[ientry] = unusedEnergy;
		chiSqs[ientry] = chiSq;
	}

	cout << "Doing the fit" << endl;
	int iter=0;
	calculateEfficiency( Metas, nentries, chiSqs, unusedEnergies, iter );



}

