
#include "/d/grid15/ln16/pi0eta/q-values/makeDiagnosticHists.h"

void makeRectSBGraphs(){
	TFile* dataFile = TFile::Open("pi0eta_data_hists_DSelector.root");
	TCanvas *allCanvases = new TCanvas("","",1440,900);
	TH1F *rectSB_4;
	TH1F *rectSB_0268;
	TH1F *rectSB_17;
	TH1F *rectSB_35;

	dataFile->GetObject("pi0eta1D_RectSBSubRegion17",rectSB_17);
	dataFile->GetObject("pi0eta1D_RectSBSubRegion35",rectSB_35);
	dataFile->GetObject("pi0eta1D_RectSBSubRegion0268",rectSB_0268);
	dataFile->GetObject("pi0eta1D_RectSBSubRegion4",rectSB_4);

	// Save the total histogram: section 4
	TH1F *rectSB_total = (TH1F*)rectSB_4->Clone();
	rectSB_total->SetName("rectSB_total");
	
	// Scale, and subtract the sidebands all together into 1 histogram. This is what we will add (+1 weight) to the total (region 4) to get the signal.
	TH1F *rectSB_background = (TH1F*)rectSB_0268->Clone();
	rectSB_background->SetName("rectSB_background");

	rectSB_background->Scale(0.5);
	rectSB_background->Add(rectSB_17,-0.5);
	rectSB_background->Add(rectSB_35,-1);

	// Simply add the background to the total since the background is properly weighted
	TH1F *rectSB_signal = (TH1F*)rectSB_total->Clone();
	rectSB_signal->SetName("rectSB_signal");

	rectSB_signal->Add(rectSB_background,1);

	// For illustration we will just invert the background to see what we a removing from the total
	rectSB_background->Scale(-1);

	makeStackedHist(rectSB_total, rectSB_signal, rectSB_background, "pi0eta1D_total","newGraphs");

}
