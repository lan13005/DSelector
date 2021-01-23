double lowMass=0.7;
double uppMass=2.0;
int nBins=65;
double stepMass=(uppMass-lowMass)/nBins;
string fitName="EtaPi_fit";

char slowMass[5];
char suppMass[5];

vector<string> group={};
void overlaySingleBin(int binN,bool isLast){
        gStyle->SetOptStat(kFALSE);
        string fitLoc = fitName+"/bin_";
        string binNum = to_string(binN);
	string outputFile = "/etapi_plot.root";
        TFile* infile = TFile::Open((fitLoc+binNum+outputFile).c_str());
        TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
        std::vector<std::string> names1D = {"Metapi","cosTheta","Phi","phi","psi","t"};

        int namesSize1D = static_cast<int>(names1D.size());
        TH1F *any1DHist_dat;
        TH1F *any1DHist_acc;
	TPaveText *pt = new TPaveText();

        double dLowMass=lowMass+binN*stepMass;
        double dUppMass=lowMass+(binN+1)*stepMass;
	sprintf(slowMass,"%.2lf",dLowMass);
	sprintf(suppMass,"%.2lf",dUppMass);
        for (int histIdx=0; histIdx<namesSize1D; ++histIdx){
   		pt->AddText(("BIN_"+to_string(binN)+" or Mass from ["+slowMass+","+suppMass+"] GeV").c_str()); 
                infile->GetObject((names1D[histIdx]+"dat").c_str(),any1DHist_dat);
                infile->GetObject((names1D[histIdx]+"acc").c_str(),any1DHist_acc);

                any1DHist_dat->Draw();
		pt->Draw();
                allCanvases->Update();

                //any1DHist_acc->SetLineColor(kGreen);
                //any1DHist_acc->SetFillStyle( 3144);
    		any1DHist_acc->SetFillColorAlpha( kOrange,0.5);

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
}

void overlayBins(){
	for (int iBin=0; iBin<nBins;++iBin){
		if(iBin<nBins-1){
			overlaySingleBin(iBin,false);
		}
		else {
			overlaySingleBin(iBin,true);
		}
	}
}
