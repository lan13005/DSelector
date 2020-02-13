#include "makeGraphs.h"
#include "/d/grid15/ln16/pi0eta/q-values/makeDiagnosticHists.h"

double selectPi0Proton=1.4;
double selectEtaProton=2;

class overlayPlots{
	private:
		std::vector<TH1F*> overlayHists;
		std::vector<double> histWeights;
		std::vector<int> histFillColors;
        	TLegend *leg = new TLegend(0.7,0.75,0.9,0.9);
		int colors[10] = {4, 6, 8, 7, 9, 30, 27, 46, 41};
		TLine* cutLine;
		int numberHists=0;
		std::map<string, double> _trackNames;
		double maximum1D = DBL_MIN;
		double minimum1D = DBL_MAX;

	public:
		overlayPlots( std::map<string, double> trackNames ){ 
			_trackNames=trackNames; 
		}
		void fillHist(TH1F* newHist){	
			if ( _trackNames.find(newHist->GetName()) != _trackNames.end() ){
				overlayHists.push_back(newHist);
				histWeights.push_back( _trackNames[newHist->GetName()] );
				++numberHists;
				if ( maximum1D < newHist->GetMaximum() ) {
					maximum1D = newHist->GetMaximum();
				}
				if ( minimum1D > newHist->GetMinimum() ) {
					minimum1D = newHist->GetMinimum();
				}
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
			overlayHists[0]->SetLineWidth( 2 ) ;
			leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
			overlayHists[0]->SetTitle("");
			overlayHists[0]->Draw("HIST");
			overlayHists[0]->SetAxisRange(minimum1D,maximum1D*1.05,"Y");
			for (int i=1; i<overlayHists.size(); ++i){
				cout << "Overlaying hist " << overlayHists[i]->GetName() << endl;
				leg->AddEntry(overlayHists[i],overlayHists[i]->GetTitle() , "l");
				overlayHists[i]->Scale(histWeights[i]);
				overlayHists[i]->SetLineColor( colors[i] );
				overlayHists[i]->SetLineWidth( 2 ) ;
				overlayHists[i]->Draw("SAME HIST");
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
			//leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
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

			//leg->Draw();
			//fileName.append(overlayHists[0]->GetName());
			//fileName.append("_overlaid.png");
			overlayCanvas->SaveAs( (fileName).c_str() );
		}

};



void makeGraphs(){
	gStyle->SetOptStat(0);
	TFile* file = TFile::Open("/d/grid15/ln16/pi0eta/092419/pi0eta_all_tLT1_2018_1_hists_DSelector.root");
	//TFile* file = TFile::Open("/d/grid15/ln16/pi0eta/092419/eta3pi/pi0eta_seanResoution_reco_3pi0_hists_DSelector.root");
	TIter keyList(file->GetListOfKeys());
	TKey *key;

	TH1F* totalHist;
	TH1F* bkgHist;
	TH1F* sigHist;
	file->GetObject("pi0eta1D_0_1_1pR_1",bkgHist);
	file->GetObject("pi0eta1D_1_1_1_1", totalHist);
	file->GetObject("pi0eta1D_1_0_mR_0",sigHist);
	makeStackedHist(totalHist, sigHist, bkgHist, "pi0eta1D_totalBkgSig","newGraphs");


		
	std::map<string, double> trackHists;
	trackHists = { {"pi0proton1D_mMandelstamT_mdelta", 1 } };
	overlayPlots pi0proton1D_mMandelstamT_mdelta(trackHists);
	trackHists = { {"etaproton1D_mMandelstamT_mdelta", 1 } };
	overlayPlots etaproton1D_mMandelstamT_mdelta(trackHists);
	trackHists = { {"pi0eta1D_RectSBSubRegion4",  1 },{"pi0eta1D_RectSBSubRegion0268",  0.25 },{"pi0eta1D_RectSBSubRegion17",  0.5 },{"pi0eta1D_RectSBSubRegion35",  0.5 } };
	overlayPlots pi0eta1D_RectSBSubRegion(trackHists);
	trackHists = { {"pi0eta1D_RectSBSubRegion4_fixed",  1 },{"pi0eta1D_RectSBSubRegion0268_fixed",  0.25 },{"pi0eta1D_RectSBSubRegion17_fixed",  0.5 },{"pi0eta1D_RectSBSubRegion35_fixed",  0.5 } };
	overlayPlots pi0eta1D_RectSBSubRegion_fixed(trackHists);
	trackHists = { {"pi0eta_mEllipsePre",1 } };
	sideBySide2D pi0eta_mEllipsePre(trackHists);	
	trackHists = { {"pi0eta_Meas_mEllipsePre",1 } };
	sideBySide2D pi0eta_Meas_mEllipsePre_showEllipse(trackHists);	
        trackHists = { {"pi0Mass_Kin_mEllipsePre",  1 }, {"pi0MassFCAL_Kin_mEllipsePre", 1 }, {"pi0MassBCAL_Kin_mEllipsePre", 1 }, {"pi0MassSPLIT_Kin_mEllipsePre", 1 } };
	overlayPlots pi0MassDiffSubDetectors(trackHists);
        trackHists = { {"etaMass_Kin_mEllipsePre",  1 }, {"etaMassFCAL_Kin_mEllipsePre", 1 }, {"etaMassBCAL_Kin_mEllipsePre", 1 }, {"etaMassSPLIT_Kin_mEllipsePre", 1 } };
	overlayPlots etaMassDiffSubDetectors(trackHists);
	trackHists = { {"pi0eta1D_mMandelstamT", 1}, {"pi0eta1D_Cut", 1} };
	overlayPlots pi0eta1DtAlltCut(trackHists);
	trackHists = { {"pi0proton1D_mMandelstamT_mdelta", 1}, {"pi0proton1D_mDelta", 1} };
	overlayPlots pi0proton1D_beforeAfterT(trackHists);

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
			if ( strcmp(h->GetName(),"pi0eta_mEllipsePre")==0 ){
				cout << "DRAWING RECTANGULAR SIDEBANDS" << endl;
				drawRectSB(0.135881-2*0.0076, 0.135881+2*0.0076, 0.548625-2*0.0191, 0.548625+2*0.0191, 0.01, 0.02);
				fileName="newGraphs/";
				fileName.append(h->GetName());
				fileName.append("_withRectCut.png");
   	   			c1->SaveAs((fileName).c_str());
			}	
			pi0eta_mEllipsePre.fillHist(h);
			pi0eta_Meas_mEllipsePre_showEllipse.fillHist(h);
		}
		else if (cl->InheritsFrom("TH1")){
			string fileName="newGraphs/";
   	   		TH1F *h = (TH1F*)key->ReadObj();
   	   		h->Draw("COLZ HIST");
			fileName.append(h->GetName());
			fileName.append(".png");
   	   		c1->SaveAs((fileName).c_str());
			// Wait until we have finally used up TH1 object first. Otherwise casting it into TH1F early creates some problems
			pi0proton1D_mMandelstamT_mdelta.fillHist(h);
			etaproton1D_mMandelstamT_mdelta.fillHist(h);
			pi0eta1D_RectSBSubRegion.fillHist(h);
			pi0eta1D_RectSBSubRegion_fixed.fillHist(h);
			pi0MassDiffSubDetectors.fillHist(h);
			etaMassDiffSubDetectors.fillHist(h);
			pi0eta1DtAlltCut.fillHist(h);
			trackHists = { {"pi0proton1D_mMandelstamT_mdelta", 1}, {"pi0proton1D_mDelta", 1} };
			if ( strcmp(h->GetName(),"pi0proton1D_mMandelstamT_mdelta")==0 ){ h->SetTitle("all t'");}
			if ( strcmp(h->GetName(),"pi0proton1D_mDelta")==0 ){ h->SetTitle("|t'|<1GeV^2");}
			pi0proton1D_beforeAfterT.fillHist(h);
		}
   	}

	std::vector<double> lineCutThresholds;
	std::vector<double [4]> cutThreshold2D; 
	lineCutThresholds={0};
	pi0eta1D_RectSBSubRegion.plot("newGraphs/pi0eta1D_RectSBSubRegions.png",false,lineCutThresholds);
	pi0eta1D_RectSBSubRegion_fixed.plot("newGraphs/pi0eta1D_RectSBSubRegions_fixed.png",false,lineCutThresholds);
	pi0MassDiffSubDetectors.plot("newGraphs/pi0MassDiffSubDetectors.png",false,lineCutThresholds);
	etaMassDiffSubDetectors.plot("newGraphs/etaMassDiffSubDetectors.png",false,lineCutThresholds);
	pi0eta1DtAlltCut.plot("newGraphs/pi0eta1DtAlltCut.png",false, lineCutThresholds);
	pi0proton1D_beforeAfterT.plot("newGraphs/pi0proton1D_beforeAfterT.png",false, lineCutThresholds);

	lineCutThresholds={selectPi0Proton};
	pi0proton1D_mMandelstamT_mdelta.plot("newGraphs/pi0proton1D_mMandelstamT_mdelta_showCut.png",true,lineCutThresholds);
	lineCutThresholds={selectEtaProton};
	etaproton1D_mMandelstamT_mdelta.plot("newGraphs/etaproton1D_mMandelstamT_mdelta_showCut.png",true,lineCutThresholds);
	
	//cutThreshold2D = { {0.134,0.538,0.013,0.04 }, {0.134,0.538,0.0155,0.05 }, {0.134,0.538, 0.0205,0.07} };
	//pi0eta_Meas_mEllipsePre_showEllipse.plot("newGraphs/pi0eta_Meas_mEllipsePre_showEllipse.png","ellipse",cutThreshold2D);
	
	cutThreshold2D = { {0.135881, 0.548625, 2*0.0076, 2*0.0191 } }; // kinFit
	//cutThreshold2D = { {0.13381, 0.5388, 3*0.006, 3*0.0264 } };//eta3pi
	pi0eta_mEllipsePre.plot("newGraphs/pi0eta_mEllipsePre_withCut.png","ellipse",cutThreshold2D);
}
