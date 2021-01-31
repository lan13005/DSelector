#include "/d/grid15/ln16/pi0eta/092419/makeStackedHists.h"
//#include "/d/grid13/gluex/gluex_top/gluex_style.C"

TH1F combineSBRegions(TFile* dataFile, string reg1Name, string reg2Name, string reg3Name, string reg4Name, string histName, string fileSaveName,
        string totLabel, string sigLabel, string bkgLabel){
	TH1F* rectSB_17;
	TH1F* rectSB_35;
	TH1F* rectSB_0268;
	TH1F* rectSB_4;
	dataFile->GetObject(reg1Name.c_str(),rectSB_17);
	dataFile->GetObject(reg2Name.c_str(),rectSB_35);
	dataFile->GetObject(reg3Name.c_str(),rectSB_0268);
	dataFile->GetObject(reg4Name.c_str(),rectSB_4);

	// Save the total histogram: section 4
	TH1F *rectSB_total = (TH1F*)rectSB_4->Clone();
	rectSB_total->SetTitle("Rectangular SBSub of M(#pi^{0}#eta) ");
	rectSB_total->SetName("rectSB_total");
	
	// Scale, and subtract the sidebands all together into 1 histogram. This is what we will add (+1 weight) to the total (region 4) to get the signal.
	TH1F *rectSB_background = (TH1F*)rectSB_0268->Clone();
	rectSB_background->SetName("rectSB_background");

	rectSB_background->Scale(0.5);
	rectSB_background->Add(rectSB_17,-0.5);
	rectSB_background->Add(rectSB_35,-1);

	// Simply add the background to the total since the background is properly weighted
	TH1F *rectSB_signal = (TH1F*)rectSB_total->Clone();
	rectSB_signal->SetName(histName.c_str());

	rectSB_signal->Add(rectSB_background,1);

	// For illustration we will just invert the background to see what we a removing from the total
	rectSB_background->Scale(-1);

	makeStackedHist(rectSB_total, rectSB_signal, rectSB_background, fileSaveName.c_str(),"newGraphs",totLabel,sigLabel,bkgLabel);
	return *rectSB_signal;
}

void makeRectSBGraphs(string fileLoc){
        //gluex_style();
	gStyle->SetOptStat(0);
	TFile* dataFile = TFile::Open(fileLoc.c_str());
	TFile* qValueFile = TFile::Open("/d/grid15/ln16/pi0eta/q-values/diagnosticPlots/postQVal.root");
	TCanvas *allCanvases = new TCanvas("","",1440,900);
	TH1F *pi0eta_Qsubbed;
	TH1F objRectSBSignal;
	TH1F* rectSB_signal;
	TH1F objRectSBSignal_fixed;
	TH1F* rectSB_signal_fixed;
	//TH1F *rectSB_4;
	//TH1F *rectSB_0268;
	//TH1F *rectSB_17;
	//TH1F *rectSB_35;

	//dataFile->GetObject("pi0eta1D_RectSBSubRegion17",rectSB_17);
	//dataFile->GetObject("pi0eta1D_RectSBSubRegion35",rectSB_35);
	//dataFile->GetObject("pi0eta1D_RectSBSubRegion0268",rectSB_0268);
	//dataFile->GetObject("pi0eta1D_RectSBSubRegion4",rectSB_4);
	objRectSBSignal = combineSBRegions( dataFile, "pi0eta1D_RectSBSubRegion17", "pi0eta1D_RectSBSubRegion35", "pi0eta1D_RectSBSubRegion0268", "pi0eta1D_RectSBSubRegion4", "rectSB_signal", "pi0eta1D_total", "Total","Signal","SB");
	rectSB_signal = &objRectSBSignal;
	objRectSBSignal_fixed = combineSBRegions( dataFile, "pi0eta1D_RectSBSubRegion17_fixed", "pi0eta1D_RectSBSubRegion35_fixed", "pi0eta1D_RectSBSubRegion0268_fixed", "pi0eta1D_RectSBSubRegion4_fixed", "rectSB_signal_fixed", "pi0eta1D_total_fixed","Total","Signal","SB");
	rectSB_signal_fixed = &objRectSBSignal_fixed;
	cout << "Made stacked histograms" << endl;

	qValueFile->GetObject("Mpi0eta_sig_kin",pi0eta_Qsubbed);

	//// Save the total histogram: section 4
	//TH1F *rectSB_total = (TH1F*)rectSB_4->Clone();
	//rectSB_total->SetTitle("Rectangular SBSub of M(#pi^{0}#eta)- ");
	//rectSB_total->SetName("rectSB_total");
	//
	//// Scale, and subtract the sidebands all together into 1 histogram. This is what we will add (+1 weight) to the total (region 4) to get the signal.
	//TH1F *rectSB_background = (TH1F*)rectSB_0268->Clone();
	//rectSB_background->SetName("rectSB_background");

	//rectSB_background->Scale(0.5);
	//rectSB_background->Add(rectSB_17,-0.5);
	//rectSB_background->Add(rectSB_35,-1);

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////// DO IT AGAIN WITH THE FIXED MASS ////////////////////////////////////


	allCanvases->Clear();	
	TLegend* leg1 = new TLegend(0.75,0.7,0.9,0.9);
	leg1->AddEntry(pi0eta_Qsubbed,"Q_{AS}","l");
	leg1->AddEntry(rectSB_signal,"SB_{#pi^{0},#eta}","l");
	//leg1->AddEntry(rectSB_signal_fixed,"SB_{#pi^{0},#eta} fixed","l");
	rectSB_signal->SetTitle("Q_{AS} vs SB_{#pi^{0},#eta}");
	pi0eta_Qsubbed->SetLineColor(kBlue);
	rectSB_signal->SetLineColor(kRed);
	rectSB_signal_fixed->SetLineColor(kGreen+2);
	cout << "Initialized legend and reset canvas" << endl;
	rectSB_signal->Draw("HIST");
	//rectSB_signal_fixed->Draw("HIST SAME");
	cout << "Drawing rectSB_signal" << endl;

	pi0eta_Qsubbed->Draw("HIST SAME");	
	cout << "Drawing pi0eta_Qsubbed" << endl;
	pi0eta_Qsubbed->SetAxisRange(rectSB_signal->GetMinimum(),rectSB_signal->GetMaximum(),"Y");
	leg1->Draw();
	allCanvases->SaveAs("newGraphs/pi0eta_sbSubQValue.png");

	allCanvases->Clear();
	pi0eta_Qsubbed->SetLineColor(kBlue);
	pi0eta_Qsubbed->SetTitle("Q_{AS}");
	pi0eta_Qsubbed->Draw("HIST");
	cout << "Drawing pi0eta_Qsubbed" << endl;
	allCanvases->SaveAs("newGraphs/pi0eta_qValue.png");
}

