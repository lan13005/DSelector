double LOW_MASS=0.7;
double UPP_MASS=2.0;
int NUMBER_BINS=60;
double lowMass=LOW_MASS;
double stepMass=(UPP_MASS-LOW_MASS)/NUMBER_BINS;
double upMass=LOW_MASS+stepMass;
char slowMass[5];
char supMass[5];

void overlaySingleBin(int binN,bool isLast){
        gStyle->SetOptStat(kFALSE);
        string fitDir = "divideRoot/bin_";
        string binNum = to_string(binN);
	string outputFile = "/twopi_plot.root";
        TFile* infile = TFile::Open((fitDir+binNum+outputFile).c_str());
        TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
        std::vector<std::string> names1D = {"M2pi","cosTheta","Phi","phi","psi","t"};

        int namesSize1D = static_cast<int>(names1D.size());
        TH1F *any1DHist_dat;
        TH1F *any1DHist_acc;
	TPaveText *pt = new TPaveText();
	sprintf(slowMass,"%.2lf",lowMass);
	sprintf(supMass,"%.2lf",upMass);
        for (int histIdx=0; histIdx<namesSize1D; ++histIdx){
   		pt->AddText(("BIN_"+to_string(binN)+" or Mass from ["+slowMass+","+supMass+"] GeV").c_str()); 
                infile->GetObject((names1D[histIdx]+"dat").c_str(),any1DHist_dat);
                infile->GetObject((names1D[histIdx]+"acc").c_str(),any1DHist_acc);

                any1DHist_dat->Draw();
		pt->Draw();
                allCanvases->Update();

                //any1DHist_acc->SetLineColor(kGreen);
    		any1DHist_acc->SetFillStyle( 3001);
    		any1DHist_acc->SetFillColor( kRed);
    		any1DHist_acc->SetLineColor( 0);
                any1DHist_acc->Draw("HIST SAME");
                
                //allCanvases->SaveAs(("overlayPlots/"+names1D[histIdx]+"/"+names1D[histIdx]+"-bin"+to_string(binN)+".png").c_str());
                if (!isLast){
                	allCanvases->Print(("overlayPlots/"+names1D[histIdx]+".pdf(").c_str(),"pdf");
		}
		else {
                	allCanvases->Print(("overlayPlots/"+names1D[histIdx]+".pdf)").c_str(),"pdf");
		}
		allCanvases->Clear();
		pt->Clear();
        }
	lowMass+=stepMass;
	upMass+=stepMass;
}

void overlayBins(){
	int maxBin=60;
	for (int iBin=0; iBin<maxBin;++iBin){
		if(iBin<maxBin-1){
			overlaySingleBin(iBin,false);
		}
		else {
			overlaySingleBin(iBin,true);
		}
	}
}
