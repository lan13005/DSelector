// I can use this to do t-slope matching 
#include "/d/grid13/gluex/gluex_top/gluex_style.C"

void makeTSlope(){
        gluex_style();
	TFile* dataFile = TFile::Open("pi0eta_data_hists_DSelector.root");
	TFile* recoFile = TFile::Open("pi0eta_flat8GeVPlus_hists_DSelector.root");

	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();
	allCanvases->Divide(2,1);

	TTree *dataTree;
	dataFile->GetObject("pi0eta_datatree_flat",dataTree);
	TTree *recoTree;
	recoFile->GetObject("pi0eta_8GeVPlustree_flat",recoTree);

	TH1F *data_tHist;
	TH1F *reco_tHist;

	dataFile->GetObject( "mandelstam_tp", data_tHist);
	recoFile->GetObject( "mandelstam_tp", reco_tHist);
	data_tHist->SetTitle("Data");
	reco_tHist->SetTitle("MC");

	allCanvases->cd(1);
	gPad->SetLogy(1);
	TF1* expFit_t;
	expFit_t = new TF1("expFit_t","expo",0.3,1);//+pol0(2)) ,tMin+tStep,tMax);
	expFit_t->SetLineStyle(2);
	data_tHist->Fit("expFit_t","RB");
	data_tHist->Draw();
	expFit_t->SetLineColor(kRed);
	expFit_t->Draw("SAME");
	
	allCanvases->cd(2);
	gPad->SetLogy(1);
	expFit_t = new TF1("expFit_t","expo",0.3,1);//+pol0(2)) ,tMin+tStep,tMax);
	expFit_t->SetLineStyle(2);
	reco_tHist->Fit("expFit_t","RB");
	reco_tHist->Draw();
	expFit_t->SetLineColor(kRed);
	expFit_t->Draw("SAME");

	allCanvases->SaveAs("tSlopePlot/tSlope.png");

}
