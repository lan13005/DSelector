double lowMass=0.7;
double uppMass=2.0;
int nBins=65;
double stepMass=(uppMass-lowMass)/nBins;
string fitName="EtaPi_fit";

char slowMass[5];
char suppMass[5];

vector<string> groups={"_S0+","_S0+_D0+","_S0+_D0+_D1+","_S0+_D0+_D1+_D2+","_P1+",""};
void overlaySingleBin(int iBin,int nBins, vector<string> names1D, vector<TCanvas*> allCanvases){
	TPaveText *pt = new TPaveText();
        double dLowMass=lowMass+iBin*stepMass;
        double dUppMass=lowMass+(iBin+1)*stepMass;
	sprintf(slowMass,"%.2lf",dLowMass);
	sprintf(suppMass,"%.2lf",dUppMass);
	pt->Clear();
   	pt->AddText(("BIN_"+to_string(iBin)+" or Mass from ["+slowMass+","+suppMass+"] GeV").c_str()); 

        gStyle->SetOptStat(kFALSE);

        string fitLoc = fitName+"/bin_";
        string binNum = to_string(iBin);
        cout << "Moved into " << fitLoc+binNum << endl;

        TH1F *any1DHist_dat;
        TH1F *any1DHist_acc;
        vector<TH1F> any1DHists_dat;
        vector<TH1F> any1DHists_acc;
        cout << "Defined some variables..." << endl;

        int igroup=1;
        for (auto group:groups){
	    string outputFile = "/etapi_plot"+group+".root";
            cout << "opening: " << fitLoc+binNum+outputFile << endl;
            TFile* infile = TFile::Open((fitLoc+binNum+outputFile).c_str());

            for (int histIdx=0; histIdx<(int)names1D.size(); ++histIdx){
                infile->GetObject((names1D[histIdx]+"dat").c_str(),any1DHist_dat);
                infile->GetObject((names1D[histIdx]+"acc").c_str(),any1DHist_acc);

                allCanvases[histIdx]->cd(igroup);
                any1DHist_dat->Draw();
                any1DHist_dat->SetTitle(group.c_str());
                any1DHist_dat->SetMinimum(0);
                allCanvases[histIdx]->Update();

                //any1DHist_acc->SetLineColor(kGreen);
                //any1DHist_acc->SetFillStyle( 3144);
    	    	any1DHist_acc->SetFillColorAlpha( kOrange,0.5);

    	    	any1DHist_acc->SetLineColor( 0);
                any1DHist_acc->Draw("HIST SAME");

                if (igroup==1){
                    // draw a pavetext showing the mass range for only the first pad
	            pt->Draw();
                }
            }
            ++igroup;
        }
        // we could have put this into the above loop but then we would have to open the same root file a lot more times
        for (int histIdx=0; histIdx<(int)names1D.size(); ++histIdx){
            if (iBin==0){
                allCanvases[histIdx]->Print(("overlayPlots/"+names1D[histIdx]+".pdf(").c_str(),"pdf");
            }
            if (iBin==(nBins-1)){
                allCanvases[histIdx]->Print(("overlayPlots/"+names1D[histIdx]+".pdf)").c_str(),"pdf");
            }
            else{
                allCanvases[histIdx]->Print(("overlayPlots/"+names1D[histIdx]+".pdf").c_str(),"pdf");
            }
        }
}

void overlayBins(){
        int ngroups=(int)groups.size();
        int nrows=(int)sqrt(ngroups);
        int ncols;
        if (nrows*nrows<ngroups)
            ncols=nrows+1;
        else
            ncols=nrows;

        TCanvas* anyCanvas;
        vector<TCanvas*> allCanvases;
        std::vector<std::string> names1D = {"Metapi","cosTheta","Phi000","Phi045","Phi090","Phi135","phi",
            "psi000","psi045","psi090","psi135","t"};
        for (auto name: names1D){
            anyCanvas = new TCanvas(("c"+name).c_str(),"",1440,900);
            anyCanvas->Divide(ncols,nrows);
            allCanvases.push_back(anyCanvas);
        }
        cout << "Defined all the canvases" << endl;

	for (int iBin=0; iBin<nBins;++iBin){
	    	overlaySingleBin(iBin,nBins,names1D,allCanvases);
	}
}
