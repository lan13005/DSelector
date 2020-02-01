void makeBaryonPlots(){
        double etaProtonBaryonCut = 1.65;
        double pi0ProtonBaryonCut = 2;
	double tCut = 1;
	TLine *lineCut;

	TFile* dataFile = TFile::Open("pi0eta_all_tLT1_hists_DSelector.root");
	TCanvas *allCanvases = new TCanvas("","",1440,900);

	string baseNames[5] = {"vanHove","mandelstam_tp","pi0proton1D", "etaproton1D", "pi0eta1D"};
	string cutString="mMandelstamT_mdelta";
	
	TH2F *vanHove;
	TH1F *mandelstam_t;
	TH1F *pi0proton;
	TH1F *etaproton;
	TH1F *anyHist;
	TH2F *any2DHist;

	// These will be the base histograms with no extra cuts applied only the base
	for (int iHist=0; iHist< 5; ++iHist){
		cout << "Making base hist: " << baseNames[iHist] << "_" << cutString << endl;
		if ( iHist == 0 ) { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), any2DHist);
			any2DHist->Draw("COLZ");
		}
		else {
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), anyHist);
			anyHist->Draw("HIST");
			if ( iHist == 3 ) {
				lineCut = new TLine(etaProtonBaryonCut,0,etaProtonBaryonCut,4000);
				lineCut->SetLineColor(kRed);
				lineCut->Draw();
			}
			if ( iHist == 1 ) {
				lineCut = new TLine(1,0,1,40000);
				lineCut->SetLineColor(kRed);
				lineCut->Draw();
			}
			if ( iHist == 2 ) {
				lineCut = new TLine(pi0ProtonBaryonCut,0,pi0ProtonBaryonCut,4000);
				lineCut->SetLineColor(kRed);
				lineCut->Draw();
			}
		}
		allCanvases->SaveAs(("baryonPlots/"+baseNames[iHist]+"_"+cutString+".pdf").c_str());
	}

	// mMandelstamT_mdelta = pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        // mMandelstamT = !pMPi0P14*mMandelstamT_mdelta; // with delta cut applied
	// mMandelstamT_mdelta_petaProton = pEtaProtonBaryonCut*mMandelstamT_mdelta; 
	// mMandelstamT_mdelta_pvanHove = pvanHove*mMandelstamT_mdelta; 
	// mDelta = ptpLT1*mMandelstamT_delta; // t Cut applied
	//
	
	allCanvases->Clear();
	TPaveText *pt = new TPaveText(.44,.46,.56,.52,"brNDC");
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	TText *t1 = pt->AddText("VANHOVE");
	t1->SetTextColor(kRed);
	allCanvases->Divide(2,2);

	int ignoreIdx;
	int padIdx;
	
	// First up would be the vanHove as the primary
	ignoreIdx=0;
	padIdx=1; 
	cutString = "mMandelstamT_mdelta_pVanHove"; 
	for (int iHist=0; iHist<5; ++iHist){
		allCanvases->cd(padIdx);
		gStyle->SetOptStat(0);
		if (iHist == ignoreIdx) { continue; } 	
		cout << "Making hists with vanHove primary: " << baseNames[iHist] << "_" << cutString << endl;
		if (iHist !=0) { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), anyHist);
			anyHist->GetXaxis()->SetLabelSize(0.05);
			anyHist->GetYaxis()->SetLabelSize(0.05);
			anyHist->GetXaxis()->SetTitleSize(0.05);
			anyHist->GetYaxis()->SetTitleSize(0.05);
			anyHist->SetBit(TH1::kNoTitle);
			anyHist->Draw("HIST");
		}
		else { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), any2DHist);
			any2DHist->GetXaxis()->SetLabelSize(0.05);
			any2DHist->GetYaxis()->SetLabelSize(0.05);
			any2DHist->GetXaxis()->SetTitleSize(0.05);
			any2DHist->GetYaxis()->SetTitleSize(0.05);
			any2DHist->SetBit(TH1::kNoTitle);
			any2DHist->Draw("COLZ");
		}
		++padIdx;
	}
	allCanvases->cd();
	pt->Draw();
	allCanvases->SaveAs(("baryonPlots/baryonCheck-"+baseNames[ignoreIdx]+".pdf").c_str());
	
	// mandelstam_tp would be the new primary
	allCanvases->Clear();
	pt = new TPaveText(.44,.46,.56,.52,"brNDC");
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	t1 = pt->AddText("t");
	t1->SetTextColor(kRed);
	allCanvases->Divide(2,2);
	ignoreIdx=1;
	padIdx=1; 
	cutString = "mDelta"; 
	for (int iHist=0; iHist<5; ++iHist){
		allCanvases->cd(padIdx);
		gStyle->SetOptStat(0);
		if (iHist == ignoreIdx) { continue; } 	
		if (iHist !=0) { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), anyHist);
			anyHist->GetXaxis()->SetLabelSize(0.05);
			anyHist->GetYaxis()->SetLabelSize(0.05);
			anyHist->GetXaxis()->SetTitleSize(0.05);
			anyHist->GetYaxis()->SetTitleSize(0.05);
			anyHist->SetBit(TH1::kNoTitle);
			anyHist->Draw("COLZ");
		}
		else {
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), any2DHist);
			any2DHist->GetXaxis()->SetLabelSize(0.05);
			any2DHist->GetYaxis()->SetLabelSize(0.05);
			any2DHist->GetXaxis()->SetTitleSize(0.05);
			any2DHist->GetYaxis()->SetTitleSize(0.05);
			any2DHist->SetBit(TH1::kNoTitle);
			any2DHist->Draw("COLZ");
		}
		++padIdx;
	}
	allCanvases->cd();
	pt->Draw();
	allCanvases->SaveAs(("baryonPlots/baryonCheck-"+baseNames[ignoreIdx]+".pdf").c_str());

	// pi0proton would be the new primary
	allCanvases->Clear();
	pt = new TPaveText(.44,.46,.56,.52,"brNDC");
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	t1 = pt->AddText("M(#pi_{0} p)");
	t1->SetTextColor(kRed);
	allCanvases->Divide(2,2);
	ignoreIdx=2;
	padIdx=1; 
	cutString = "mMandelstamT"; 
	for (int iHist=0; iHist<5; ++iHist){
		allCanvases->cd(padIdx);
		gStyle->SetOptStat(0);
		if (iHist == ignoreIdx) { continue; } 	
		if (iHist !=0) { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), anyHist);
			anyHist->GetXaxis()->SetLabelSize(0.05);
			anyHist->GetYaxis()->SetLabelSize(0.05);
			anyHist->GetXaxis()->SetTitleSize(0.05);
			anyHist->GetYaxis()->SetTitleSize(0.05);
			anyHist->SetBit(TH1::kNoTitle);
			anyHist->Draw("COLZ");
		}
		else {
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), any2DHist);
			any2DHist->GetXaxis()->SetLabelSize(0.05);
			any2DHist->GetYaxis()->SetLabelSize(0.05);
			any2DHist->GetXaxis()->SetTitleSize(0.05);
			any2DHist->GetYaxis()->SetTitleSize(0.05);
			any2DHist->SetBit(TH1::kNoTitle);
			any2DHist->Draw("COLZ");
		}
		++padIdx;
	}
	allCanvases->cd();
	pt->Draw();
	allCanvases->SaveAs(("baryonPlots/baryonCheck-"+baseNames[ignoreIdx]+".pdf").c_str());

	// etaproton would be the new primary
	allCanvases->Clear();
	pt = new TPaveText(.44,.46,.56,.52,"brNDC");
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	t1 = pt->AddText("M(#eta p)");
	t1->SetTextColor(kRed);

	allCanvases->Divide(2,2);
	ignoreIdx=3;
	padIdx=1; 
	cutString = "mMandelstamT_mdelta_petaProton"; 
	for (int iHist=0; iHist<5; ++iHist){
		allCanvases->cd(padIdx);
		gStyle->SetOptStat(0);
		if (iHist == ignoreIdx) { continue; } 	
		if (iHist !=0) { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), anyHist);
			anyHist->GetXaxis()->SetLabelSize(0.05);
			anyHist->GetYaxis()->SetLabelSize(0.05);
			anyHist->GetXaxis()->SetTitleSize(0.05);
			anyHist->GetYaxis()->SetTitleSize(0.05);
			anyHist->SetBit(TH1::kNoTitle);
			anyHist->Draw("COLZ");
		}
		else {
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), any2DHist);
			any2DHist->GetXaxis()->SetLabelSize(0.05);
			any2DHist->GetYaxis()->SetLabelSize(0.05);
			any2DHist->GetXaxis()->SetTitleSize(0.05);
			any2DHist->GetYaxis()->SetTitleSize(0.05);
			any2DHist->SetBit(TH1::kNoTitle);
			any2DHist->Draw("COLZ");
		}
		++padIdx;
	}
	allCanvases->cd();
	pt->Draw();
	allCanvases->SaveAs(("baryonPlots/baryonCheck-"+baseNames[ignoreIdx]+".pdf").c_str());
}
