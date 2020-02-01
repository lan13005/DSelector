
#include "/d/grid15/ln16/pi0eta/q-values/makeDiagnosticHists.h"

void makeRectSBGraphs(){
	gStyle->SetOptStat(0);
	TFile* dataFile = TFile::Open("pi0eta_all_tLT1_hists_DSelector.root");
	TFile* qValueFile = TFile::Open("/d/grid15/ln16/pi0eta/q-values/postQVal.root");
	TCanvas *allCanvases = new TCanvas("","",1440,900);
	TH1F *rectSB_4;
	TH1F *rectSB_0268;
	TH1F *rectSB_17;
	TH1F *rectSB_35;
	TH1F *pi0eta_Qsubbed;

	dataFile->GetObject("pi0eta1D_RectSBSubRegion17",rectSB_17);
	dataFile->GetObject("pi0eta1D_RectSBSubRegion35",rectSB_35);
	dataFile->GetObject("pi0eta1D_RectSBSubRegion0268",rectSB_0268);
	dataFile->GetObject("pi0eta1D_RectSBSubRegion4",rectSB_4);

	qValueFile->GetObject("Mpi0eta_sig_kin",pi0eta_Qsubbed);

	// Save the total histogram: section 4
	TH1F *rectSB_total = (TH1F*)rectSB_4->Clone();
	rectSB_total->SetTitle("Rectangular SBSub of M(#pi^{0}#eta)- ");
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

	allCanvases->Clear();	
	TLegend* leg1 = new TLegend(0.65,0.7,0.9,0.9);
	leg1->AddEntry(pi0eta_Qsubbed,"Q(as,sb)+SB(pi0)","l");
	leg1->AddEntry(rectSB_signal,"SB(pi0,eta)","l");
	rectSB_signal->SetTitle("Q(as,sb)+SB(pi0) vs SB(pi0,eta)");
	pi0eta_Qsubbed->SetLineColor(kBlue);
	rectSB_signal->SetLineColor(kRed);
	rectSB_signal->Draw("HIST");
	pi0eta_Qsubbed->Draw("HIST SAME");	
	pi0eta_Qsubbed->SetAxisRange(rectSB_signal->GetMinimum(),rectSB_signal->GetMaximum(),"Y");
	leg1->Draw();
	allCanvases->SaveAs("newGraphs/pi0eta_sbSubQValue.png");

	allCanvases->Clear();
	pi0eta_Qsubbed->SetLineColor(kBlue);
	pi0eta_Qsubbed->SetTitle("Q(as,sb)+SB(pi0)");
	pi0eta_Qsubbed->Draw("HIST");
	allCanvases->SaveAs("newGraphs/pi0eta_qValue.png");

}
