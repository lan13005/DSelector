void makeBaryonPlots(){
	TFile* dataFile = TFile::Open("pi0eta_data_hists_DSelector.root");
	TCanvas *allCanvases = new TCanvas("","",1440,900);

	string baseNames[5] = {"vanHove","mandelstam_tp","pi0proton1D", "etaproton1D", "pi0eta1D"};
	string cutString="mMandelstamT_delta";
	
	TH2F *vanHove;
	TH1F *mandelstam_t;
	TH1F *pi0proton;
	TH1F *etaproton;
	TH1F *anyHist;
	TH2F *any2DHist;

	// These will be the base histograms with no extra cuts applied only the base
	for (int iHist=0; iHist< 5; ++iHist){
		dataFile->GetObject((baseNames[i]+"_"+cutString).c_str(), anyHist);
		anyHist->Draw("COLZ");
		allCanvases->SaveAs(("baryonPlots/"+baseNames[i]+"_"+cutString+".pdf").c_str());
	}

	// mMandelstamT_mdelta = pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        // mMandelstamT = !pMPi0P14*mMandelstamT_mdelta; // with delta cut applied
	// mMandelstamT_mdelta_petaProton = pEtaProtonBaryonCut*mMandelstamT_mdelta; 
	// mMandelstamT_mdelta_pvanHove = pvanHove*mMandelstamT_mdelta; 
	// mDelta = ptpLT1*mMandelstamT_delta; // t Cut applied
	allCanvases->Divide(2,2);
	int ignoreIdx;
	int padIdx;
	
	// First up would be the vanHove as the primary
	ignoreIdx=0;
	padIdx=1; 
	allCanvases->cd(padIdx);
	cutString = "mMandelstamT_mdelta_pvanHove"; 
	for (int iHist=0; iHist<5; ++iHist){
		if (iHist == ignoreIdx) { continue; } 	
		if (iHist !=0) { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), anyHist);
			anyHist->Draw("COLZ");
		}
		else { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), any2DHist);
			any2DHist->Draw("COLZ");
		}
		++padIdx;
	}
	allCanvases->SaveAs(("baryonPlots/baryonCheck-"+baseNames[ignoreIdx]+".pdf").c_str());
	
	// mandelstam_tp would be the new primary
	ignoreIdx=1;
	padIdx=1; 
	allCanvases->cd(padIdx);
	cutString = "mDelta"; 
	for (int iHist=0; iHist<5; ++iHist){
		if (iHist == ignoreIdx) { continue; } 	
		if (iHist !=0) { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), anyHist);
			anyHist->Draw("COLZ");
		}
		else {
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), any2DHist);
			any2DHist->Draw("COLZ");
		}
		++padIdx;
	}
	allCanvases->SaveAs(("baryonPlots/baryonCheck-"+baseNames[ignoreIdx]+".pdf").c_str());

	// pi0proton would be the new primary
	ignoreIdx=2;
	padIdx=1; 
	allCanvases->cd(padIdx);
	cutString = "mMandelstamT"; 
	for (int iHist=0; iHist<5; ++iHist){
		if (iHist == ignoreIdx) { continue; } 	
		if (iHist !=0) { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), anyHist);
			anyHist->Draw("COLZ");
		}
		else {
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), any2DHist);
			any2DHist->Draw("COLZ");
		}
		++padIdx;
	}
	allCanvases->SaveAs(("baryonPlots/baryonCheck-"+baseNames[ignoreIdx]+".pdf").c_str());

	// etaproton would be the new primary
	ignoreIdx=3;
	padIdx=1; 
	allCanvases->cd(padIdx);
	cutString = "mMandelstamT_mdelta_petaProton"; 
	for (int iHist=0; iHist<5; ++iHist){
		if (iHist == ignoreIdx) { continue; } 	
		if (iHist !=0) { 
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), anyHist);
			anyHist->Draw("COLZ");
		}
		else {
			dataFile->GetObject((baseNames[iHist]+"_"+cutString).c_str(), any2DHist);
			any2DHist->Draw("COLZ");
		}
		++padIdx;
	}
	allCanvases->SaveAs(("baryonPlots/baryonCheck-"+baseNames[ignoreIdx]+".pdf").c_str());
}
