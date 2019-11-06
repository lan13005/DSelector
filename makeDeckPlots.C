int numDOFsig=8;
Double_t g2(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	Double_t r2 = Double_t((x[1]-par[3])/par[4]);
	return par[5]+par[6]*(x[0]-par[1])+par[7]*(x[1]-par[3])+par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}

int numDOFsigFlat=6;
Double_t gflat2(Double_t *x, Double_t *par) {
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	Double_t r2 = Double_t((x[1]-par[3])/par[4]);
	return par[5]+par[0]/2/TMath::Pi()/par[2]/par[4]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}

int numDOFsig1D=5;
Double_t gaus(Double_t *x, Double_t *par){
	Double_t r1 = Double_t((x[0]-par[1])/par[2]);
	return par[3]+par[4]*(x[0]-par[1])+par[0]/sqrt(2*TMath::Pi())/par[2]*TMath::Exp(-0.5*r1*r1 );
}

void makeDeckPlots(){
	gSystem->Exec("rm -rf deckPlots");
	gSystem->Exec("mkdir -p deckPlots/mandelstam_teta_meas");
	gSystem->Exec("mkdir deckPlots/mandelstam_tpi0_meas");
    	ofstream logFile;
    	logFile.open("deckPlots/failedFittingPlotIDs.txt");

	TFile* dataFile = TFile::Open("pi0eta_datatreeFlat_DSelector.root");
	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	TCanvas *allCanvases_yields = new TCanvas("anyHists_yields","",1440,900);
	allCanvases_yields->Divide(3,4,0,0);

	TTree *dataTree;
	dataFile->GetObject("pi0eta_datatree_flat",dataTree);

	double new_t;
	double new_m;
	double new_accWeight;
	Bool_t new_isUnique12B_1234B;
	Bool_t new_isUnique34B_1234B;
	Bool_t isUnique;
	double new_meta;
	double new_mpi0;
        dataTree->SetBranchAddress("AccWeight",&new_accWeight);
        dataTree->SetBranchAddress("isNotRepeated_pi0_pi0eta",&new_isUnique12B_1234B);
        dataTree->SetBranchAddress("isNotRepeated_eta_pi0eta",&new_isUnique34B_1234B);
	dataTree->SetBranchAddress("Meta_meas",&new_meta);
	dataTree->SetBranchAddress("Mpi0_meas",&new_mpi0);
	dataTree->SetBranchAddress("Mpi0eta_meas",&new_m);
	cout << "Loaded addresses" << endl;

	//string branchNames[2]={"mandelstam_tpi0_meas","mandelstam_teta_meas"};
	string branchNames[2]={"mandelstam_teta_meas","mandelstam_tpi0_meas"};
	int counter=-1;
	for ( string branchName: branchNames ) {
		logFile << branchName << endl;
		++counter;
		allCanvases->cd();
		dataTree->SetBranchAddress( branchName.c_str(),&new_t);

		int num_tBins=14;
		double tMin=0;
		double tMax=2.8;
		int num_massBins=12;
		double mMin=1.7;
		double mMax=2.9;
		double tStep=(tMax-tMin)/num_tBins;
		double mStep=(mMax-mMin)/num_massBins;

		// So the following indicies would describe the bin for a specific variable. We have to multiplex them into a single array. 
		// The obvious way to do that is to use the array index = num_tBins*m+t
		// This would order the array where the first [0,num_tBins] would belong in the smallest mass bin, [num_tBins,2*num_tBins] belongs to the second smallest.
		int idx_t;
		int idx_m;
		const int numHists = (const int)num_tBins*num_massBins;
		cout << "numHists: " << numHists << endl;
		TH2F* hists_meta_mpi0[numHists];
		TH2F* hists_mpi0eta_t[numHists];
		double yields[numHists];
		double yieldErrors[numHists];
		TH2F* full_meta_mpi0 = new TH2F(("full_meta_mpi0_"+std::to_string(counter)).c_str(), "Cuts=GeneralCuts;M(#pi_{0}) GeV;M(#eta)(GeV)", 100,0.05,0.25,100,0.25,0.85 );

		double xfitMin=0.085;
		double xfitMax=0.185;
		double yfitMin=0.4;
		double yfitMax=0.64;
		double chiSq;
		double dof;

		double yBins[2] = {0.25,0.85};
		double xBins[2] = {0.05,0.25};
		int nxBins=100;
		int nyBins=100;
		double xWidth = (xBins[1]-xBins[0])/nxBins;
		double yWidth = (yBins[1]-yBins[0])/nyBins;
		double differentialArea = 1/(xWidth*yWidth);
		cout << "Differential Area:"<< differentialArea << endl;
		TH1F* chiSqPerDOF = new TH1F(("chiSqPerDOF_"+std::to_string(counter)).c_str(), "#chi^{2}/n_{dof}; #chi^{2}/n_{dof}; Events / 0.1", 30,0,3 );
		TH1F* full_meta = new TH1F(("full_meta_"+std::to_string(counter)).c_str(), ";M(#eta)(GeV)", nxBins,yBins[1],yBins[0] );
		TH1F* full_mpi0 = new TH1F(("full_mpi0_"+std::to_string(counter)).c_str(), ";M(#pi_0)(GeV)",nyBins,xBins[1],xBins[0] );
		
		for (int i=0; i< numHists; ++i){
			hists_meta_mpi0[i] = new TH2F(("meta_mpi0"+to_string(i)+"_"+std::to_string(counter)).c_str(), "Cuts=GeneralCuts;M(#pi_{0}) GeV;M(#eta)(GeV)", 100,0.05,0.25,100,0.25,0.85 ); 
			hists_mpi0eta_t[i] = new TH2F(("mpi0eta_t"+to_string(i)+"_"+std::to_string(counter)).c_str(), "Cuts=mMandelstamT_eta;M(#pi_{0}#eta) (GeV);t_{#eta} (GeV^2)", 260, 0.6, 3.2, 80,0,8);
		}
		int histIdx;
		cout << "Initialized histograms" << endl;

		Long64_t nentries = dataTree->GetEntries();

		for (int ientry=0; ientry<nentries; ++ientry){
			dataTree->GetEntry(ientry);
			if (counter==0){ isUnique = new_isUnique34B_1234B; }
			else { isUnique = new_isUnique12B_1234B; } 
			if ( new_t < tMin || new_m < mMin || new_t>tMax || new_m>mMax ) {continue; }
			idx_t = (int)( (new_t-tMin)/tStep ); 
			idx_m = (int)( (new_m-mMin)/mStep );
			//cout << idx_t << ", " << idx_m << endl;
			histIdx=num_tBins*idx_m+idx_t;
			hists_meta_mpi0[histIdx]->Fill(new_mpi0,new_meta,new_accWeight);
			hists_mpi0eta_t[histIdx]->Fill(new_m,new_t,new_accWeight);
			full_meta_mpi0->Fill(new_mpi0,new_meta,new_accWeight);
			full_meta->Fill(new_meta,new_accWeight);
			full_mpi0->Fill(new_mpi0,new_accWeight);
		}
		cout << "Loaded all the data" << endl;

		///////////////////////////////////////////
		// DOING 1D FIT TO USE AS INITIALIZATION FOR 2D
		///////////////////////////////////////////
		cout << "-----------------------\nSTARTING FULL ETA FIT!\n---------------------------- " << endl;
		double fitPars_eta[5];
		int fit_nentries=full_meta_mpi0->GetEntries();
		cout << "TOTAL ENTRIES: " << fit_nentries << endl;
		gStyle->SetOptFit();
		TF1 * feta = new TF1("meta_fit",gaus,yfitMin,yfitMax,numDOFsig1D); 
		feta->SetParameters(2000,0.55,0.03,0, 0);
		feta->SetParLimits(0,0,5000);
		feta->SetParLimits(1,0.52,0.585);
		feta->SetParLimits(2,0,0.1);
		full_meta->Fit(feta,"RLV");
		feta->GetParameters(fitPars_eta);
		allCanvases->SaveAs(("deckPlots/"+branchName+"/project_full_meta.png").c_str());

		cout << "-----------------------\nSTARTING FULL PI0 FIT!\n---------------------------- " << endl;
		double fitPars_pi0[5];
		allCanvases->Clear();
		TF1 * fpi0 = new TF1("mpi0_fit",gaus,xfitMin,xfitMax,numDOFsig1D); 
		fpi0->SetParameters(500,0.135,0.015,0, 0);
		fpi0->SetParLimits(0,0,1000);
		fpi0->SetParLimits(1,0.125,0.145);
		fpi0->SetParLimits(2,0.005,0.02);
		full_mpi0->Fit(fpi0,"RLV");
		fpi0->GetParameters(fitPars_pi0);
		allCanvases->SaveAs(("deckPlots/"+branchName+"/project_full_mpi0.png").c_str());

		cout << "-----------------------\nSTARTING FULL PI0ETA FIT!\n---------------------------- " << endl;
		TF2 * f2 = new TF2("f2",g2,xfitMin,xfitMax,yfitMin,yfitMax,numDOFsig); 
		TF2 * fflat2 = new TF2("fflat2",gflat2,xfitMin,xfitMax,yfitMin,yfitMax,numDOFsigFlat); 
		double peakBin=1200;
		// We will try to estimate the initialization of the 2D gaussian using the 1D gaussians
		// we get the peak from looking at the histogram. Take the peak*2*pi*sigmaX*sigmaY
		// take directly the mean, sigma from the pi0 and eta 1D fits since they dont scale
		// for the flat bkg we can take some avg between the two const parameters, and for the linear we can take the values from the individual fits. These 3 have to be divided by 100 since
		// 	it might make sense the scale grows with the number of bins we project onto it. In this case the number of bins is 100.
		double full_initParams[8]={peakBin*2*TMath::Pi()*fitPars_pi0[2]*fitPars_eta[2],fitPars_pi0[1], fitPars_pi0[2], fitPars_eta[1], fitPars_eta[2], (fitPars_eta[3]+fitPars_pi0[3])/2/100, fitPars_eta[4]/100, fitPars_pi0[4]/100};
		double full_convParams[8];
		double fitPars[8];
		double fitParErrors[8];
		const double *pointParError;
		f2->SetParameters(full_initParams[0],full_initParams[1],full_initParams[2],full_initParams[3],full_initParams[4],full_initParams[5],full_initParams[6],full_initParams[7]); 
		f2->SetParLimits(0,0,full_initParams[0]*10);
		f2->SetParLimits(1,0.13,0.14);
		f2->SetParLimits(2,0.005,0.015);
		f2->SetParLimits(3,0.52,0.585);
		f2->SetParLimits(4,0,0.1);
		allCanvases->Clear();
		f2->SetLineColorAlpha(kRed, 0.6);
		f2->SetLineWidth(1);
		full_meta_mpi0->Fit(f2,"RLV");
		f2->GetParameters(full_convParams);
		full_meta_mpi0->Draw("COLZ");
		//full_meta_mpi0->Draw("SURF2");
		//f2->Draw("SURF SAME");
		allCanvases->SaveAs(("deckPlots/"+branchName+"/full_meta_mpi0.png").c_str());	
		allCanvases->Clear();


		cout << "-----------------------\nSTARTING BINNED FITS!\n---------------------------- " << endl;


		///////////////////////////////////////////
		// STARTING THE INDIVIDUAL FITS IN EACH BIN
		///////////////////////////////////////////
		double scaled_initParams[8];
		double scaleFactor;
		double maxYield=DBL_MIN;
		TRandom randomGenerator;
		int maxResamplings=10;
		double scales[3]={0.1,0.2,0.3};
		int scaleThresholds[3]={2,4,6};
		for(int i=0; i<numHists; ++i){
        		scaleFactor = hists_meta_mpi0[i]->GetEntries()/fit_nentries;
			for (int j=0; j<8; ++j){
				if ( j==0 || j==5 || j==6 || j==7){
					// the only paramters that scale are the massHist of  the 2D gaus and the paramters for the planar bkg
					scaled_initParams[j] = full_convParams[j]*scaleFactor;
				}
				else {
					// We cant scale everything since the mean and variance shouldn't change
					scaled_initParams[j] = full_convParams[j];
				}
			}
			f2->SetParameters(scaled_initParams[0],scaled_initParams[1],scaled_initParams[2],scaled_initParams[3],scaled_initParams[4],scaled_initParams[5],scaled_initParams[6],scaled_initParams[7]); 
			f2->SetParLimits(0,0,hists_meta_mpi0[i]->GetEntries()); // just restrict to be positive
			f2->SetParLimits(1,0.13,0.14);
			f2->SetParLimits(2,0.005,0.015);
			f2->SetParLimits(3,0.52,0.585);
			f2->SetParLimits(4,0,0.1);
			//f2->SetParLimits(5,0,scaled_initParams[5]*10);
			//f2->SetParLimits(6,-scaled_initParams[5]*10,scaled_initParams[5]*10);
			//f2->SetParLimits(7,-scaled_initParams[6]*10,scaled_initParams[6]*10);
			cout << "\n\n*******************\nFITTING NEXT BIN************************"<<endl;
			Int_t fitStatus = hists_meta_mpi0[i]->Fit(f2,"RLE");

			if ( fitStatus != 0 ) {  
				f2->SetParameters(scaled_initParams[0],scaled_initParams[1],scaled_initParams[2],scaled_initParams[3],scaled_initParams[4],scaled_initParams[5],0,0); 
				f2->FixParameter(6,0);
				f2->FixParameter(7,0);
				fitStatus = hists_meta_mpi0[i]->Fit(f2,"RLE");
				if (fitStatus != 0) { cout << "FUDGE WE STILL NEED TO FIX THIS.... EXITING..." << endl; logFile << i << endl; }
			}
			//f2->SetParameter(0,randomGenerator.Uniform( scaled_initParams[0]*(1-scale),scaled_initParams[0]*(1+scale) ) );
			//f2->SetParameter(5,randomGenerator.Uniform( scaled_initParams[5]*(1-scale),scaled_initParams[5]*(1+scale) ) );
			


			// just resample the amplitudes and copy over the peak,widths
			//double scale;
			//for ( int maxIters=0; maxIters<maxResamplings; ++maxIters ){
			//	if (fitStatus == 0 ) { break; }
			//	if (maxIters==scaleThresholds[0]) { scale = scales[0]; } 
			//	if (maxIters==scaleThresholds[1]) { scale = scales[1]; } 
			//	if (maxIters==scaleThresholds[2]) { scale = scales[2]; } 
			//	cout << "                  *NEXT RESAMPLE*"<<endl;
			//	cout << randomGenerator.Uniform( scaled_initParams[0]*(1-scale),scaled_initParams[0]*(1+scale)) << endl;
			//	f2->SetParameter(0,randomGenerator.Uniform( scaled_initParams[0]*(1-scale),scaled_initParams[0]*(1+scale) ) );
			//	f2->SetParameter(5,randomGenerator.Uniform( scaled_initParams[5]*(1-scale),scaled_initParams[5]*(1+scale) ) );
			//	f2->SetParameter(6,randomGenerator.Uniform( scaled_initParams[6]*(1-scale),scaled_initParams[6]*(1+scale) ) );
			//	f2->SetParameter(7,randomGenerator.Uniform( scaled_initParams[7]*(1-scale),scaled_initParams[7]*(1+scale) ) );
			//	fitStatus = hists_meta_mpi0[i]->Fit(f2,"RLEV");
			//}
			//if (fitStatus != 0) { cout << "FUDGE WE STILL NEED TO FIX THIS.... EXITING..." << endl; exit(0); }
			
			dof=f2->GetNDF();
			chiSq=f2->GetChisquare();
			f2->GetParameters(fitPars);
			pointParError = f2->GetParErrors();

			allCanvases->Clear();
			chiSqPerDOF->Fill(chiSq/dof);
			hists_meta_mpi0[i]->SetTitle(("scaleFactor:"+std::to_string(scaleFactor)+"  ChiSqPerDOF: "+std::to_string(chiSq/dof)).c_str());
			hists_meta_mpi0[i]->Draw("COLZ");
			allCanvases->SaveAs(("deckPlots/"+branchName+"/meta_mpi0-"+to_string(i)+".png").c_str());	
			//allCanvases->Clear();
			//hists_mpi0eta_t[i]->Draw("COLZ");
			//allCanvases->SaveAs(("deckPlots/mpi0eta_t-"+to_string(i)+".png").c_str());	
			
			yields[i]=differentialArea*fitPars[0];
			yieldErrors[i]=differentialArea*pointParError[0];

			cout << "Yield, error = " << fitPars[0] << ", " << pointParError[0];
			cout << "   Scaled Yield, Scaled error = " << differentialArea*fitPars[0] << ", " << differentialArea*pointParError[0];
			if (maxYield<differentialArea*fitPars[0]){
				maxYield=differentialArea*fitPars[0];
			}
		}
		cout << "    maxYield = " << maxYield << endl;

		allCanvases->Clear();
		chiSqPerDOF->Draw();
		allCanvases->SaveAs(("deckPlots/"+branchName+"/chiSq.png").c_str());

		std::vector<string> massBinTitles; 
		for (Int_t i=0; i<num_massBins; ++i){
			massBinTitles.push_back( std::to_string( mMin+i*mStep )+" < M(#pi_{0}#eta) < " + std::to_string( mMin+(i+1)*mStep ) );
		}
		TH1F* massHist;
		gStyle->SetErrorX(0.000001); // remove the x-error bars
		for (Int_t iAmp=0; iAmp<massBinTitles.size(); ++iAmp){
			cout << "Plotting histogram for mass bin: " << iAmp << endl;
			massHist = new TH1F(("massHist_"+std::to_string(iAmp)+"_"+std::to_string(counter)).c_str(),massBinTitles[iAmp].c_str(),num_tBins,tMin,tMax);	
			massHist->SetAxisRange(0,maxYield*1.1,"Y");
			allCanvases_yields->cd( iAmp+1 );	
			cout << "Plotting on pad : " << iAmp+1 << endl;
			for (Int_t i=0; i<num_tBins; ++i){
				massHist->SetBinContent( i+1, yields[i+num_tBins*iAmp] );
				massHist->SetBinError( i+1, yieldErrors[i+num_tBins*iAmp]);
				//cout << "yield in bin: " << i+num_tBins*iAmp << " is " << yields[i+num_tBins*iAmp] << endl;
			}
			massHist->SetMarkerStyle(kFullCircle);
			if (counter==0){
				massHist->SetMarkerColor(kCyan-2);
				massHist->Draw("E1 PMC");
			}
			else { 
				massHist->SetMarkerColor(kRed-3);
				massHist->Draw("E1 PMC SAME");
			}
		}

		cout << "Initial parameter values: " << endl;
		for ( auto parValue: full_initParams ) {
			cout << parValue << endl;
		}

		cout << "ALL THE YIELDS: " <<endl;
		for (auto yield : yields ) {
			cout << yield << endl;
		}
	}
	gStyle->SetOptStat(0);
	allCanvases_yields->SaveAs("deckPlots/yields.png");



}
