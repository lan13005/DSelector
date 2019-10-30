double selectPi0Proton=1.4;

class overlayPlots{
	private:
		std::vector<TH1F*> overlayHists;
		std::vector<EColor> colors={kBlue,kMagenta,kGreen};
        	TLegend *leg = new TLegend(0.1,0.8,0.3,0.9);
		TLine* cutLine;
		int numberHists=0;
		std::set<string> _trackNames;

	public:
		overlayPlots( std::set<string> trackNames ){ 
			_trackNames=trackNames; 
		}
		void fillHist(TH1F* newHist){	
			if ( _trackNames.find(newHist->GetName()) != _trackNames.end() ){
				overlayHists.push_back(newHist);
				++numberHists;
				assert(numberHists<=3);
			}
		}
		
		void drawVLine( double xCut ){
			cout << "   Drawing VLine at " << xCut << " with maximum " << overlayHists[0]->GetMaximum() << endl;
			cutLine = new TLine(xCut,0,xCut,overlayHists[0]->GetMaximum());
			cutLine->SetLineWidth(3);
			cutLine->SetLineColor(kRed);	
			cutLine->SetLineStyle(9);
			cutLine->Draw("SAME");
		}

		void plot(string fileName, bool b_xCut, double xCut) {
			cout << "Creating hist " << overlayHists[0]->GetName() << endl;
			TCanvas* overlayCanvas = new TCanvas( ("canvas_"+fileName).c_str(),"",1440,900);
			overlayHists[0]->SetLineColor(colors[0]);
			leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
			overlayHists[0]->Draw();
			for (int i=1; i<overlayHists.size(); ++i){
				cout << "Overlaying hist " << overlayHists[i]->GetName() << endl;
				overlayHists[i]->SetLineColor(colors[i]);
				overlayHists[i]->Draw("same");
				leg->AddEntry(overlayHists[i],overlayHists[i]->GetTitle() , "l");
			}
			if (b_xCut){
				drawVLine(xCut);
				//cout << "   Drawing VLine at " << xCut << " with maximum " << overlayHists[0]->GetMaximum() << endl;
				//cutLine = new TLine(xCut,0,xCut,overlayHists[0]->GetMaximum());
				//cutLine->SetLineWidth(6);
				//cutLine->SetLineColor(kRed);	
				//cutLine->Draw("SAME");
			}			

			leg->Draw();
			//fileName.append(overlayHists[0]->GetName());
			//fileName.append("_overlaid.png");
			overlayCanvas->SaveAs( (fileName).c_str() );
		}

};


void makeGraphs(){
	TFile* file = TFile::Open("pi0eta_data_hists_DSelector.root");
	TIter keyList(file->GetListOfKeys());
	TKey *key;
		
	std::set<string> trackHists;
        trackHists = {"pi0eta1D_Cut_tAll","pi0eta1D_1_0_mR_0"};
	overlayPlots pi0eta1D_tAllVsSub(trackHists);
	trackHists = {"pi0proton1D_Cut_ASBS"};
	overlayPlots pi0proton1D_Cut_ASBS(trackHists);

	TCanvas *c1 = new TCanvas("c1","",1440,900);
	int i=0;
   	while ((key = (TKey*)keyList())) {
   	   	TClass *cl = gROOT->GetClass(key->GetClassName());
   	   	if (cl->InheritsFrom("TH1")){
			string fileName="newGraphs/";
   	   		TH1F *h = (TH1F*)key->ReadObj();
   	   		h->Draw("COLZ HIST");
			fileName.append(h->GetName());
			fileName.append(".png");
   	   		c1->SaveAs((fileName).c_str());

			// Wait until we have finally used up TH1 object first. Otherwise casting it into TH1F early creates some problems
			pi0eta1D_tAllVsSub.fillHist(h);
			pi0proton1D_Cut_ASBS.fillHist(h);
		}
   	}
	pi0eta1D_tAllVsSub.plot("newGraphs/pi0eta1D_tAllVsSub.png",false,0);
	pi0proton1D_Cut_ASBS.plot("newGraphs/pi0proton1D_Cut_ASBS.png",true,selectPi0Proton);


}
