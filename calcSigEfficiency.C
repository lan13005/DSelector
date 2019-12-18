#include "/d/grid15/ln16/pi0eta/q-values/main.h"

double calculateEfficiency( double var1s[], double chiSqs[], double unusedEnergies[] ) {
	TH1F *hist_var;
        double binRange[3] = {100,0.35,0.8};
        double fitRange[2] = {0.4,0.65};
	double par[numDOFbkg+numDOFsig];
	double parBkg[numDOFbkg];
	Long64_t nentries = sieof(var1)/sizeof(var1[0]);

	hist_var = new TH1F("","",binRange[0],binRange[1],binRange[2]);
	fit_eta = new TF1("fit_eta",fitFunc,fitRange[0],fitRange[1],numDOFbkg+numDOFsig);
	fit_eta->SetParameters(4600,8600,20000,0.5406,0.0237);

	for ( auto var1 : var1s ){
		hist_var->Fill(var1s);
	}

	fit_eta->Fit("fit_eta","RQB");
	fit->GetParameters(par);
	bkgFit->SetParameters(par);
	sigFit->SetParameters(&par[numDOFbkg]);

	double nSig = sigFit->Integral(binRange[0],binRange[1]);
	double nBkg = bjgFit->Integral(binRange[0],binRange[1]);

	return nSig/nBkg; 
}


void calcSigEfficiency(){
	TFile* dataFile=new TFile("pi0eta_datatreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_datatree_flat",dataTree);

	Long64_t nentries=dataTree->GetEntries();
	double Metas[nentries];
	double chiSqs[nentries];
	double unusedEnergies[nentries];
	double Meta;
	double chiSq;
	double unusedEnergy;
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("chiSq",&chiSq);
	dataTree->SetBranchAddress("unusedEnergy",&unusedEnergy);

	for (Long64_t ientry=0; ientry<nentries; ientry++)
	{
        	dataTree->GetEntry(ientry);
		Metas[ientry] = Meta;
		unusedEnergies[ientry] = unusedEnergy;
		chiSqs[ientry] = chiSq;
	}
}

