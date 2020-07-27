// 1d double gaussian
int numDOFsig=5;
Double_t signal(Double_t *x, Double_t *par){
	return par[0]/par[2]/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2]))
	     + par[3]*par[0]/(par[4]*par[2])/TMath::Sqrt( 2*TMath::Pi() )*exp(-0.5*((x[0]-par[1])/(par[4]*par[2]))*((x[0]-par[1])/(par[4]*par[2])));

}

int numDOFbkg = 2;
Double_t background(Double_t *x, Double_t *par){
	return par[0]+par[1]*x[0];
}

Double_t fitFunc(Double_t *x, Double_t *par){
	return background(x,par)+signal(x,&par[numDOFbkg]);
}


bool verbose=true;
int maxAttempts = 150;
int minEntries = 200;

void getResolutionCompareDataMC(){
	TFile* fileReco = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_compare_reco_2017_hists_DSelector.root");
	TFile* fileData = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_compare_data_2017_hists_DSelector.root");
        vector<TFile*> files={fileReco,fileData};
        vector<string> fileNames={"MC","Data"};
        Int_t color1[3]={kAzure-4,kOrange+6, kGreen+2};
        Int_t color2[3]={kAzure-4,kOrange+6, kGreen+2};
        Int_t color3[3]={kAzure-4,kOrange+6, kGreen+2};
        Int_t markerShapes[3]={20,21,22};

	gStyle->SetOptStat(0);

	TCanvas *allCanvases = new TCanvas("allCanvases","",1440,900);
	allCanvases->Divide(2,1);

	TH1F *anyHist;
	TH1F *pi0Hist;
	string particles[2] = {"eta","pi0"};
        string prepend = "Mass";
        string affix = "_Kin_mEllipsePre_thetaBin";
        static const int nDetectors=3;
        string detectors[nDetectors] = {"BCAL","FCAL","SPLIT"};
	gStyle->SetTitleSize(0.08,"t");

	std::vector<double> fitRangeEta;
	std::vector<double> fitRangePi0;
	fitRangeEta={0.44,0.64};
	fitRangePi0={0.10,0.17};

        double chiSq;
	TF1* fit;
	TF1* bkgFit;
	TF1* sigFit;
	Double_t par[numDOFbkg+numDOFsig];
    
	if(verbose)cout << "Starting to loop through all the histograms" << endl;
	TMultiGraph *mg_eta = new TMultiGraph("mg_eta","");
	TMultiGraph *mg_pi0 = new TMultiGraph("mg_pi0","");
	auto legend1 = new TLegend(0.05,0.1,0.9,0.3);
        TPaveText *pt = new TPaveText(0.05,0.1,.9,0.3,"NDC");
	//auto legend2 = new TLegend(0.425,0.75,0.575,0.9);

        static const int nThetaBins=5;
        double thetaBinWidth=40/nThetaBins;

        TRandom* jitter = new TRandom();
	std::vector<TGraphErrors*> gr_etaWidths[2];
	std::vector<TGraphErrors*> gr_pi0Widths[2];
	vector<double> widths[2][2][nDetectors]; // [file][particle][detector][nThetaBins]
        vector<double> widthErrs[2][2][nDetectors]; 
        vector<double> thetas[2][2][nDetectors];
        vector<double> thetaErrs[2][2][nDetectors];
        double theta;
        double thetaErr;
        
        int nFile=0;
        for (auto file: files){
            if(verbose)cout << " ---------------------------- " << endl;
            if(verbose)cout << " Starting file: " << fileNames[nFile] << endl;
            if(verbose)cout << " ---------------------------- " << endl;
	    for (int thetaBin=0; thetaBin<nThetaBins; ++ thetaBin){
	    	if(verbose)cout << "    Begin work in theta bin: " << thetaBin << endl;
	    	TCanvas *allCanvases_fits = new TCanvas("allCanvases_fits","",1440,900);
	    	allCanvases_fits->Divide(2,1);
                    int iDetector=0;
                    for ( auto detector : detectors) {
	    	        int particleCounter = 0;
	    	        if(verbose)cout << "      Begin work with detector: " << detector << endl;
	    	        for ( auto particle : particles){
                            allCanvases_fits->cd(particleCounter+1);
	    	            string name = particle+prepend+detector+affix+to_string(thetaBin);
	    		    if(verbose)cout << "      Getting object: " << name << endl;
	    		    file->GetObject(name.c_str(), anyHist);
	    		    if (particleCounter==0){
	    		    	fit = new TF1("fit",fitFunc,fitRangeEta[0],fitRangeEta[1],numDOFsig+numDOFbkg); 
	    		    	bkgFit = new TF1("bkgFit",background,fitRangeEta[0],fitRangeEta[1],numDOFbkg); 
	    		    	sigFit = new TF1("sigFit",signal,fitRangeEta[0],fitRangeEta[1],numDOFsig); 
	    		    }
	    		    else{
	    		    	fit = new TF1("fit",fitFunc,fitRangePi0[0],fitRangePi0[1],numDOFsig+numDOFbkg); 
	    		    	bkgFit = new TF1("bkgFit",background,fitRangePi0[0],fitRangePi0[1],numDOFbkg); 
	    		    	sigFit = new TF1("sigFit",signal,fitRangePi0[0],fitRangePi0[1],numDOFsig); 
	    		    }
                                gStyle->SetOptFit(1);
                                gStyle->SetStatY(.9);

	    		    double par0;
	    		    double par1;
	    		    double par2;
                                
	    		    double weightedSigma;
                            double weightedSigmaErr;
	    		    TLine *drawFitBounds = new TLine();
                            if(verbose)cout << "Beggining fit attempts" << endl;
	    		    for (int iAttempt=0; iAttempt < maxAttempts; ++iAttempt){
	    		    	int totalEntries = anyHist->GetEntries();
                                if(verbose)cout << " - entries: " << totalEntries << endl;
                                if ( totalEntries < minEntries ) {
                                    if(verbose)cout << " -- skipping since total entries < number of bins! Setting width=0, chiSq=0" << endl;
                                    weightedSigma=0;
                                    weightedSigmaErr=0;
                                    chiSq=0;
                                    break;
                                }
                                // in case there is very little number of events we can scale things up and back down to prevent integer rounding
	    		    	int scaledEntries = 1000.0*totalEntries/anyHist->GetNbinsX();
                                if(verbose)cout << " - scaledEntries: " << scaledEntries << endl;
	    		    	par0 = rand()%scaledEntries;
	    		    	par1 = rand()%scaledEntries;
	    		    	par2 = rand()%scaledEntries;
                                par0=par0/1000;
                                par1=par1/1000;
                                par2=par2/1000;
	    		    	
	    		    	if (particleCounter==0){
	    		    		fit->SetParameters(par0,par1,par2,0.55,0.01,2,2);
	    		    		fit->SetParLimits(3,0.525,0.575);
	    		    		fit->SetParLimits(4,0.01,0.05);
                                        fit->SetParLimits(6,0.25,4);
	    		    	}
	    		    	else{
	    		    		fit->SetParameters(par0,par1,par2,0.134,0.01,1,0.5);
	    		    		fit->SetParLimits(3,0.1275,0.1425);
	    		    		fit->SetParLimits(4,0.005,0.01);
                                        fit->SetParLimits(6,0.5,2);
	    		    	}
	    		    	fit->SetParLimits(0,0,par0*3);
	    		    	fit->SetParLimits(2,0,par2*3);
	    		    	fit->SetParLimits(5,0,3);

                                Int_t fitStatus;
	    		    	fitStatus = anyHist->Fit(fit,"RQE0");
	    		        chiSq=fit->GetChisquare()/(fit->GetNDF());

	    		    	fit->GetParameters(par);

                                double sig1=par[2+numDOFbkg];
                                double B=par[3+numDOFbkg];
                                double C=par[4+numDOFbkg];
                                double dsig1=fit->GetParError(2+numDOFbkg);
                                double dB=fit->GetParError(3+numDOFbkg);
                                double dC=fit->GetParError(4+numDOFbkg);
                                if(verbose)cout << "sig1/err: " << sig1 << "/" << dsig1 << endl;
                                if(verbose)cout << "B/err: " << B << "/" << dB << endl;
                                if(verbose)cout << "C/err: " << C << "/" << dC << endl;
                                weightedSigma=(1+B*C)/(1+B)*sig1;
                                double dwdB=weightedSigma*(C-1)/(1+B*C)/(1+B);
                                double dwdC=B*weightedSigma/(1+B*C);
                                double dwdsig1=weightedSigma/sig1;
                                weightedSigmaErr=TMath::Sqrt(dwdB*dwdB*dB*dB+dwdC*dwdC*dC*dC+dwdsig1*dwdsig1*dsig1*dsig1);

                                double startChi=10;
                                double chiStep=10;
                                int attemptWidth=20;
                                double chiTheshold=iAttempt/attemptWidth*chiStep+startChi;
                                if(verbose)cout << " - chiTheshold at " << chiTheshold << endl;

                                cout << "fitStatus " << fitStatus << " reducedChi " << fit->GetChisquare()/(fit->GetNDF()) << " weightedSigmaErr/weightedSigma " << weightedSigmaErr/weightedSigma << endl;
	    		    	if(fitStatus == 0 && fit->GetChisquare()/(fit->GetNDF()) < chiTheshold && weightedSigmaErr/weightedSigma < 0.4){
	    		    		cout << "------------- FOUND CONVERGENCE CRITERIA AFTER " << iAttempt << " ITERATIONS" << endl;
                                        theta=thetaBinWidth*(thetaBin+0.5)+jitter->Gaus(0,0.5);
                                        thetaErr=0;
	    		                widths[nFile][particleCounter][iDetector].push_back(weightedSigma);
	    		                widthErrs[nFile][particleCounter][iDetector].push_back(weightedSigmaErr);
                                        thetas[nFile][particleCounter][iDetector].push_back(theta);
                                        thetaErrs[nFile][particleCounter][iDetector].push_back(thetaErr);
	    		    		break;
	    		    	}
	    		    	if (iAttempt==maxAttempts-1) {
	    		    		cout << "------------- DID NOT FIND CONVERGENCE CRITERIA AFTER " << maxAttempts << " ITERATIONS" << endl;
                                        weightedSigma=0;
                                        weightedSigmaErr=0;
                                        chiSq=0;
                                        //exit(0);
	    		    	}
	    		    }
	    		    anyHist->Draw();

	    	            // include a legend on just on of the pads and use that as a reference
	    	            std::stringstream streamChi;
                            std::stringstream streamSig;
                            std::stringstream streamSigErr;
	    	            streamChi << std::fixed << std::setprecision(4) << chiSq;
	    	            std::string streamedChi = streamChi.str();

	    	            streamSig << std::fixed << std::setprecision(4) << weightedSigma;
	    	            std::string streamedSig = streamSig.str();

	    	            streamSigErr << std::fixed << std::setprecision(4) << weightedSigmaErr;
	    	            std::string streamedSigErr = streamSigErr.str();
                                
                            anyHist->SetTitle((particle+" "+detector+" chiSq/DOF: "+streamedChi+" weightedSigma/Err: "+streamedSig+"/"+streamedSigErr+
                                            "    "+to_string(thetaBin*thetaBinWidth)+" < #theta < "+to_string((thetaBin+1)*thetaBinWidth)).c_str());

	    		    fit->Draw("SAME");
                                
	    	            allCanvases_fits->SaveAs(("temp/"+fileNames[nFile]+"_"+detector+"_thetaBin"+to_string(thetaBin)+".png").c_str());
	    	            ++particleCounter;
                        }
                        ++iDetector;
	    	}
	    }
	    	
            if(verbose)cout << "Drawing the final resolution vs theta graphs" << endl;
            cout << "Num converged eta " << fileNames[nFile] << " " << detectors[0] << " " << thetas[nFile][0][0].size() << endl;
            cout << "Num converged eta " << fileNames[nFile] << " " << detectors[1] << " " << thetas[nFile][0][1].size() << endl;
            cout << "Num converged eta " << fileNames[nFile] << " " << detectors[2] << " " << thetas[nFile][0][2].size() << endl;
            cout << "Num converged pi0 " << fileNames[nFile] << " " << detectors[0] << " " << thetas[nFile][1][0].size() << endl;
            cout << "Num converged pi0 " << fileNames[nFile] << " " << detectors[1] << " " << thetas[nFile][1][1].size() << endl;
            cout << "Num converged pi0 " << fileNames[nFile] << " " << detectors[2] << " " << thetas[nFile][1][2].size() << endl;

	    // We will fill up all the results for all the chiSq bins for a specific Mpi0eta bin
	    gr_etaWidths[nFile].push_back(new TGraphErrors(thetas[nFile][0][0].size(),&(thetas[nFile][0][0][0]),&(widths[nFile][0][0][0]),&(thetaErrs[nFile][0][0][0]),&(widthErrs[nFile][0][0][0])));
	    gr_etaWidths[nFile].push_back(new TGraphErrors(thetas[nFile][0][1].size(),&(thetas[nFile][0][1][0]),&(widths[nFile][0][1][0]),&(thetaErrs[nFile][0][1][0]),&(widthErrs[nFile][0][1][0])));
	    gr_etaWidths[nFile].push_back(new TGraphErrors(thetas[nFile][0][2].size(),&(thetas[nFile][0][2][0]),&(widths[nFile][0][2][0]),&(thetaErrs[nFile][0][2][0]),&(widthErrs[nFile][0][2][0])));
	    gr_etaWidths[nFile][0]->SetMarkerStyle(markerShapes[0]+4*nFile);
            gr_etaWidths[nFile][0]->SetMarkerSize(1.5);
            gr_etaWidths[nFile][0]->SetMarkerColor(color1[0]);
	    gr_etaWidths[nFile][1]->SetMarkerStyle(markerShapes[1]+4*nFile);
            gr_etaWidths[nFile][1]->SetMarkerSize(1.5);
            gr_etaWidths[nFile][1]->SetMarkerColor(color2[1]);
	    gr_etaWidths[nFile][2]->SetMarkerStyle(markerShapes[2]+4*nFile);
            gr_etaWidths[nFile][2]->SetMarkerSize(1.5);
            gr_etaWidths[nFile][2]->SetMarkerColor(color3[2]);
	    mg_eta->Add(gr_etaWidths[nFile][0],"p");
	    mg_eta->Add(gr_etaWidths[nFile][1],"p");
	    mg_eta->Add(gr_etaWidths[nFile][2],"p");

	    gr_pi0Widths[nFile].push_back(new TGraphErrors(thetas[nFile][1][0].size(),&(thetas[nFile][1][0][0]),&(widths[nFile][1][0][0]),&(thetaErrs[nFile][1][0][0]),&(widthErrs[nFile][1][0][0])));
	    gr_pi0Widths[nFile].push_back(new TGraphErrors(thetas[nFile][1][1].size(),&(thetas[nFile][1][1][0]),&(widths[nFile][1][1][0]),&(thetaErrs[nFile][1][1][0]),&(widthErrs[nFile][1][1][0])));
	    gr_pi0Widths[nFile].push_back(new TGraphErrors(thetas[nFile][1][2].size(),&(thetas[nFile][1][2][0]),&(widths[nFile][1][2][0]),&(thetaErrs[nFile][1][2][0]),&(widthErrs[nFile][1][2][0])));
	    gr_pi0Widths[nFile][0]->SetMarkerStyle(markerShapes[0]+4*nFile);
            gr_pi0Widths[nFile][0]->SetMarkerColor(color1[0]);
            gr_pi0Widths[nFile][0]->SetMarkerSize(1.5);
	    gr_pi0Widths[nFile][1]->SetMarkerStyle(markerShapes[1]+4*nFile);
            gr_pi0Widths[nFile][1]->SetMarkerColor(color2[1]);
            gr_pi0Widths[nFile][1]->SetMarkerSize(1.5);
	    gr_pi0Widths[nFile][2]->SetMarkerStyle(markerShapes[2]+4*nFile);
            gr_pi0Widths[nFile][2]->SetMarkerColor(color3[2]);
            gr_pi0Widths[nFile][2]->SetMarkerSize(1.5);
	    mg_pi0->Add(gr_pi0Widths[nFile][0],"p");
	    mg_pi0->Add(gr_pi0Widths[nFile][1],"p");
	    mg_pi0->Add(gr_pi0Widths[nFile][2],"p");

            legend1->AddEntry(gr_etaWidths[nFile][0],(fileNames[nFile]+"_BCAL").c_str());
            legend1->AddEntry(gr_etaWidths[nFile][1],(fileNames[nFile]+"_FCAL").c_str());
            legend1->AddEntry(gr_etaWidths[nFile][2],(fileNames[nFile]+"_SPLIT").c_str());

            ++nFile; // increment the file index
        }
	allCanvases->cd(1);
        gPad->SetBottomMargin(0.4);
	mg_eta->Draw("A pmc pc");
	mg_eta->SetTitle("#eta mass resolution vs theta");
	mg_eta->GetXaxis()->SetTitle("#theta + Gaus(0,0.5) (degrees)");
        legend1->Draw("SAME");
	allCanvases->cd(2);
        gPad->SetBottomMargin(0.4);
        pt->AddText("Comparing MC/data resolutions extracted as a weighted sigma from a double gaussian fit");
        pt->AddText(("Skip resolution extraction if nentries<"+to_string(minEntries)+" or if fit failed "+to_string(maxAttempts)+" times").c_str());
        // SetAllWith - seems liek it uses the first argument to match to specific lines: http://151.100.123.6/studenti/cmp01/root-61804/tutorials/geom/csgdemo.C
        pt->SetAllWith("","size",0.02);
	mg_pi0->Draw("A pmc pc");
	mg_pi0->SetTitle("#pi mass resolution vs theta");
	mg_pi0->GetXaxis()->SetTitle("#theta + Gaus(0,0.5) (degrees)");
        //legend2->Draw("SAME");
        pt->Draw();
	allCanvases->SaveAs("temp/resolutionVsTheta.png");

}



