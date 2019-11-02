double selectPi0Proton=1.4;

class overlayPlots{
	private:
		std::vector<TH1F*> overlayHists;
		std::vector<double> histWeights;
        	TLegend *leg = new TLegend(0.1,0.8,0.3,0.9);
		int colors[10] = {4, 6, 8, 7, 9, 30, 27, 46, 41};
		TLine* cutLine;
		int numberHists=0;
		std::map<string, double> _trackNames;

	public:
		overlayPlots( std::map<string, double> trackNames ){ 
			_trackNames=trackNames; 
		}
		void fillHist(TH1F* newHist){	
			if ( _trackNames.find(newHist->GetName()) != _trackNames.end() ){
				overlayHists.push_back(newHist);
				histWeights.push_back( _trackNames[newHist->GetName()] );
				++numberHists;
				//assert(numberHists<=3);
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

		void plot(string fileName, bool b_xCut, std::vector<double> xCut) {
			cout << "Creating hist " << overlayHists[0]->GetName() << endl;
			TCanvas* overlayCanvas = new TCanvas( ("canvas_"+fileName).c_str(),"",1440,900);
			overlayHists[0]->Scale(histWeights[0]);
			overlayHists[0]->SetLineColor( colors[0] );
			leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
			overlayHists[0]->Draw();
			for (int i=1; i<overlayHists.size(); ++i){
				cout << "Overlaying hist " << overlayHists[i]->GetName() << endl;
				overlayHists[i]->Scale(histWeights[i]);
				overlayHists[i]->SetLineColor( colors[i] );
				overlayHists[i]->Draw("same");
				leg->AddEntry(overlayHists[i],overlayHists[i]->GetTitle() , "l");
			}
			if (b_xCut){
				for (auto i : xCut){
					drawVLine(i);
				}
			}			

			leg->Draw();
			//fileName.append(overlayHists[0]->GetName());
			//fileName.append("_overlaid.png");
			overlayCanvas->SaveAs( (fileName).c_str() );
		}

};


 //  Side by side comparison
class sideBySide2D{
	private:
		std::vector<TH2F*> overlayHists;
		std::vector<double> histWeights;
        	TLegend *leg = new TLegend(0.1,0.8,0.3,0.9);
		TEllipse* cutEllipse;
		int numberHists=0;
		std::map<string, double> _trackNames;

	public:
		sideBySide2D( std::map<string, double> trackNames ){ 
			_trackNames=trackNames; 
		}
		void fillHist(TH2F* newHist){	
			if ( _trackNames.find(newHist->GetName()) != _trackNames.end() ){
				overlayHists.push_back(newHist);
				histWeights.push_back( _trackNames[newHist->GetName()] );
				++numberHists;
			}
		}
		
		void drawEllipse( double x, double y, double xr, double yr ){
			cutEllipse = new TEllipse(x,y,xr,yr);
			cutEllipse->SetLineWidth(3);
			cutEllipse->SetLineColor(kRed);	
        		cutEllipse->SetFillStyle(0);
			cutEllipse->SetLineStyle(9);
			cutEllipse->Draw("SAME");
		}

		void plot(string fileName, string cutShape, std::vector<double [4]> xCut) {
			cout << "Creating hist " << overlayHists[0]->GetName() << endl;
			TCanvas* overlayCanvas = new TCanvas( ("canvas_"+fileName).c_str(),"",1440,900);
			overlayHists[0]->Scale(histWeights[0]);
			leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
			overlayHists[0]->Draw("COLZ");
			//for (int i=1; i<overlayHists.size(); ++i){
			//	cout << "Overlaying hist " << overlayHists[i]->GetName() << endl;
			//	overlayHists[i]->Scale(histWeights[i]);
			//	overlayHists[i]->SetLineColor( colors[i] );
			//	overlayHists[i]->Draw("same");
			//	leg->AddEntry(overlayHists[i],overlayHists[i]->GetTitle() , "l");
			//}
			if (cutShape=="ellipse"){
				for (auto i : xCut){
					drawEllipse(i[0],i[1],i[2],i[3]);
				}
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
		
	std::map<string, double> trackHists;
        trackHists = { {"pi0eta1D_Cut_tAll", 1}, {"pi0eta1D_1_0_mR_0",1} };
	overlayPlots pi0eta1D_tAllVsSub(trackHists);
	trackHists = { {"pi0proton1D_Cut_ASBS",1} };
	overlayPlots pi0proton1D_Cut_ASBS(trackHists);
	trackHists = { {"pi0eta1D_RectSBSubRegion4", 1},{"pi0eta1D_RectSBSubRegion0268", 0.25},{"pi0eta1D_RectSBSubRegion17", 0.5},{"pi0eta1D_RectSBSubRegion35", 0.5} };
	overlayPlots pi0eta1D_RectSBSubRegion(trackHists);
	trackHists = { {"pi0eta_Meas_mEllipsePre",1 } };
	sideBySide2D pi0eta_Meas_mEllipsePre(trackHists);	

	TCanvas *c1 = new TCanvas("c1","",1440,900);
	int i=0;
   	while ((key = (TKey*)keyList())) {
   	   	TClass *cl = gROOT->GetClass(key->GetClassName());
   	   	if (cl->InheritsFrom("TH2")){
			string fileName="newGraphs/";
   	   		TH2F *h = (TH2F*)key->ReadObj();
   	   		h->Draw("COLZ HIST");
			fileName.append(h->GetName());
			fileName.append(".png");
   	   		c1->SaveAs((fileName).c_str());
			pi0eta_Meas_mEllipsePre.fillHist(h);
		}
		else if (cl->InheritsFrom("TH1")){
			string fileName="newGraphs/";
   	   		TH1F *h = (TH1F*)key->ReadObj();
   	   		h->Draw("COLZ HIST");
			fileName.append(h->GetName());
			fileName.append(".png");
   	   		c1->SaveAs((fileName).c_str());
			// Wait until we have finally used up TH1 object first. Otherwise casting it into TH1F early creates some problems
			pi0eta1D_tAllVsSub.fillHist(h);
			pi0proton1D_Cut_ASBS.fillHist(h);
			pi0eta1D_RectSBSubRegion.fillHist(h);
		}
   	}

	std::vector<double> lineCutThresholds;
	std::vector<double [4]> cutThreshold2D; 
	lineCutThresholds={0};
	pi0eta1D_tAllVsSub.plot("newGraphs/pi0eta1D_tAllVsSub.png",false,lineCutThresholds);
	pi0eta1D_RectSBSubRegion.plot("newGraphs/pi0eta1D_RectSBSubRegions.png",false,lineCutThresholds);
	lineCutThresholds={selectPi0Proton};
	pi0proton1D_Cut_ASBS.plot("newGraphs/pi0proton1D_Cut_ASBS.png",true,lineCutThresholds);
	
	cutThreshold2D = { {0.134,0.538,0.013,0.04 }, {0.134,0.538,0.0155,0.05 }, {0.134,0.538, 0.0205,0.07} };
	pi0eta_Meas_mEllipsePre.plot("newGraphs/pi0eta_Meas_mEllipsePre_withCut.png","ellipse",cutThreshold2D);
}
