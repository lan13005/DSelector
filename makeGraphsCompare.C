#include "makeGraphs.h"
#include "/d/grid15/ln16/pi0eta/q-values/makeDiagnosticHists.h"

double selectPi0Proton=1.4;
double selectEtaProton=2;

bool overlay1D=false;

class overlayPlots{
	private:
		std::vector<TH1F*> overlayHists;
		std::vector<TH1*> overlayHists2;
		std::vector<double> histWeights;
		std::vector<int> histFillColors;
        	TLegend *leg = new TLegend(0.7,0.75,0.9,0.9);
		int colors[10] = {4, 6, 8, 7, 9, 30, 27, 46, 41};
		TLine* cutLine;
		int numberHists=0;
		std::map<string, double> _trackNames;
		double maximum1D = DBL_MIN;
		double minimum1D = DBL_MAX;
		double maximum1D2 = DBL_MIN;
		double minimum1D2 = DBL_MAX;

	public:
		overlayPlots( std::map<string, double> trackNames ){ 
			_trackNames=trackNames; 
		}
		void fillHist(TH1F* newHist, TH1* newHist2){	
			if ( _trackNames.find(newHist->GetName()) != _trackNames.end() ){
				overlayHists.push_back(newHist);
				overlayHists2.push_back(newHist2);
				double weight = _trackNames[newHist->GetName()];
				histWeights.push_back( weight );
				++numberHists;
				if ( maximum1D < newHist->GetMaximum()*weight ) {
					maximum1D = newHist->GetMaximum()*weight;
				}
				if ( minimum1D > newHist->GetMinimum()*weight ) {
					minimum1D = newHist->GetMinimum()*weight;
				}
				if ( maximum1D2 < newHist2->GetMaximum()*weight ) {
					maximum1D2 = newHist2->GetMaximum()*weight;
				}
				if ( minimum1D2 > newHist2->GetMinimum()*weight ) {
					minimum1D2 = newHist2->GetMinimum()*weight;
				}
				cout << "Current maximum: " << maximum1D2 << endl; 
				cout << "Current minimum: " << minimum1D2 << endl; 
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
			overlayCanvas->Divide(2,1);
			overlayCanvas->cd(1);
			overlayHists[0]->Scale(histWeights[0]);
			cout << "Scaling original hist by " << histWeights[0] << endl;
			overlayHists[0]->SetLineColor( colors[0] );
			overlayHists[0]->SetLineWidth( 2 ) ;
			leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
			overlayHists[0]->SetTitle("");
			overlayHists[0]->Draw("HIST");
			overlayHists[0]->GetXaxis()->SetTitleSize(0.05);
			overlayHists[0]->GetYaxis()->SetTitleSize(0.05);
			overlayHists[0]->SetAxisRange(minimum1D,maximum1D*1.05,"Y");
			for (int i=1; i<overlayHists.size(); ++i){
				cout << "Overlaying hist " << overlayHists[i]->GetName() << endl;
				leg->AddEntry(overlayHists[i],overlayHists[i]->GetTitle() , "l");
				overlayHists[i]->Scale(histWeights[i]);
				cout << " Scaled by " << histWeights[i] << endl;
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
			//

			leg->Clear();
			overlayCanvas->cd(2);
			overlayHists2[0]->Scale(histWeights[0]);
			cout << "Scaling original hist by " << histWeights[0] << endl;
			overlayHists2[0]->SetLineColor( colors[0] );
			overlayHists2[0]->SetLineWidth( 2 ) ;
			leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
			overlayHists2[0]->SetTitle("");
			overlayHists2[0]->Draw("HIST");
			overlayHists2[0]->GetXaxis()->SetTitleSize(0.05);
			overlayHists2[0]->GetYaxis()->SetTitleSize(0.05);
			overlayHists2[0]->SetAxisRange(minimum1D2,maximum1D2*1.05,"Y");
			cout << "Setting axis range for right panel: " << minimum1D2 << ", " << maximum1D2 << endl;
			for (int i=1; i<overlayHists2.size(); ++i){
				cout << "Overlaying hist " << overlayHists2[i]->GetName() << endl;
				leg->AddEntry(overlayHists2[i],overlayHists2[i]->GetTitle() , "l");
				overlayHists2[i]->Scale(histWeights[i]);
				cout << " Scaled by " << histWeights[i] << endl;
				overlayHists2[i]->SetLineColor( colors[i] );
				overlayHists2[i]->SetLineWidth( 2 ) ;
				overlayHists2[i]->Draw("SAME HIST");
			}
			if (b_xCut){
				for (auto i : xCut){
					drawVLine(i);
				}
			}			
			leg->Draw();

			overlayCanvas->SaveAs( (fileName).c_str() );
		}

};

 //  Side by side comparison
class sideBySide2D{
	private:
		std::vector<TH2F*> overlayHists;
		std::vector<TH2*> overlayHists2;
		std::vector<double> histWeights;
        	TLegend *leg = new TLegend(0.1,0.8,0.3,0.9);
		TEllipse* cutEllipse;
		int numberHists=0;
		std::map<string, double> _trackNames;

	public:
		sideBySide2D( std::map<string, double> trackNames ){ 
			_trackNames=trackNames; 
		}
		void fillHist(TH2F* newHist, TH2* newHist2){	
			if ( _trackNames.find(newHist->GetName()) != _trackNames.end() ){
				overlayHists.push_back(newHist);
				overlayHists2.push_back(newHist2);
				histWeights.push_back( _trackNames[newHist->GetName()] );
				++numberHists;
			}
		}
		
		void drawEllipse( double x, double y, double xr, double yr ){
			cout << "makiing ellipse" << endl;
			cutEllipse = new TEllipse(x,y,xr,yr);
			cutEllipse->SetLineWidth(3);
			cutEllipse->SetLineColor(kRed);	
        		cutEllipse->SetFillStyle(0);
			cutEllipse->SetLineStyle(9);
			cout << "Drawing ellipse" << endl;
			cutEllipse->Draw("SAME");
			cout << "Drew ellipse" << endl;
		}

		void plot(string fileName, string cutShape, std::vector<double [4]> xCut) {
			cout << "Creating hist " << overlayHists[0]->GetName() << endl;
			TCanvas* overlayCanvas = new TCanvas( ("canvas_"+fileName).c_str(),"",1440,900);
			overlayCanvas->Divide(2,1);
			overlayCanvas->cd(1);
			overlayHists[0]->Scale(histWeights[0]);
			overlayHists[0]->Draw("COLZ");
			//leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
			//overlayHists[0]->GetXaxis()->SetTitleSize(0.05);
			//overlayHists[0]->GetYaxis()->SetTitleSize(0.05);
			//for (int i=1; i<overlayHists.size(); ++i){
			//	cout << "Overlaying hist " << overlayHists[i]->GetName() << endl;
			//	overlayHists[i]->Scale(histWeights[i]);
			//	overlayHists[i]->SetLineColor( colors[i] );
			//	overlayHists[i]->Draw("same");
			//	leg->AddEntry(overlayHists[i],overlayHists[i]->GetTitle() , "l");
			//}
			if (cutShape=="ellipse"){
				cout << "Drawing first ellipse" << endl;
				for (auto i : xCut){
					cout << "i = {" <<i[0] <<","<<i[1] << "," << i[2] << "," << i[3] << "}" << endl;
					drawEllipse(i[0],i[1],i[2],i[3]);
				}
			}			
			overlayCanvas->cd(2);
			cout << "going to second pad" << endl;
			overlayHists2[0]->Scale(histWeights[0]);
			cout << "Scaled the histogram" << endl;
			overlayHists2[0]->Draw("COLZ");
			cout << " Drew base histogram" << endl;
			if (cutShape=="ellipse"){
				cout << "Drawing second ellipse" << endl;
				for (auto i : xCut){
					cout << "i2 = {" <<i[0] <<","<<i[1] << "," << i[2] << "," << i[3] <<  "}" << endl;
					drawEllipse(i[0],i[1],i[2],i[3]);
				}
			}			

			//leg->Draw();
			//fileName.append(overlayHists[0]->GetName());
			//fileName.append("_overlaid.png");
			overlayCanvas->SaveAs( (fileName).c_str() );
		}

};



void makeGraphsCompare(){	
	TFile* file = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_test_hists_DSelector.root");
	TFile* file2 = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_test_defUT_hists_DSelector.root");

	gStyle->SetOptStat(0);
	TIter keyList(file->GetListOfKeys());
	TKey *key;

	TH1F* totalHist;
	TH1F* bkgHist;
	TH1F* sigHist;

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
	trackHists = { {"pi0proton1D_baseAsymCut", 1}, {"pi0proton1D_baseAsymCut_mDelta_fastEta", 1} };
	overlayPlots pi0proton1D_baseAsym(trackHists);
	trackHists = { {"etaproton1D_baseAsymCut", 1}, {"etaproton1D_baseAsymCut_mDelta_fastPi0", 1} };
	overlayPlots etaproton1D_baseAsym(trackHists);

	TCanvas *c1D = new TCanvas("c1D","",1440,900);
	if ( overlay1D ) {
		c1D->Divide(1,1);
	} else {
		c1D->Divide(2,1);
	}
	TCanvas *c2D = new TCanvas("c2D","",1440,900);
	c2D->Divide(2,1);
	int i=0;
   	while ((key = (TKey*)keyList())) {
   	   	TClass *cl = gROOT->GetClass(key->GetClassName());
   	   	if (cl->InheritsFrom("TH2")){
			string fileName="newGraphs/";
   	   		TH2F *h = (TH2F*)key->ReadObj();
			h->GetXaxis()->SetTitleSize(0.04);
			h->GetYaxis()->SetTitleSize(0.04);
			c2D->cd(1);
   	   		h->Draw("COLZ HIST");
			c2D->cd(2);
			TH2 *h2;
			file2->GetObject(h->GetName(),h2);
			h2->Draw("COLZ HIST");
			fileName.append(h->GetName());
			fileName.append(".png");
   	   		c2D->SaveAs((fileName).c_str());
			if ( strcmp(h->GetName(),"pi0eta_mEllipsePre")==0 ){
				cout << "DRAWING RECTANGULAR SIDEBANDS" << endl;
				drawRectSB(0.135881-2*0.0076, 0.135881+2*0.0076, 0.548625-2*0.0191, 0.548625+2*0.0191, 0.01, 0.02);
				fileName="newGraphs/";
				fileName.append(h->GetName());
				fileName.append("_withRectCut.png");
   	   			c2D->SaveAs((fileName).c_str());
			}	
			pi0eta_mEllipsePre.fillHist(h,h2);
			pi0eta_Meas_mEllipsePre_showEllipse.fillHist(h,h2);
		}
		else if (cl->InheritsFrom("TH1")){
			string fileName="newGraphs/";
   	   		TH1F *h = (TH1F*)key->ReadObj();
			h->GetXaxis()->SetTitleSize(0.04);
			h->GetYaxis()->SetTitleSize(0.04);
			c1D->cd(1);
   	   		h->Draw("HIST");
			TH1* h2;
			file2->GetObject(h->GetName(),h2);
			if (overlay1D) {
				h2->SetLineColor(kGreen+1);
				h2->Scale(0.4);
				h2->Draw("SAME");
				h2->SetLineColor(kBlue);
			} else { 
				c1D->cd(2);
				h2->Draw("HIST");
			}
			fileName.append(h->GetName());
			fileName.append(".png");
   	   		c1D->SaveAs((fileName).c_str());
			// Wait until we have finally used up TH1 object first. Otherwise casting it into TH1F early creates some problems
			pi0proton1D_mMandelstamT_mdelta.fillHist(h,h2);
			etaproton1D_mMandelstamT_mdelta.fillHist(h,h2);
			pi0proton1D_baseAsym.fillHist(h,h2);
			etaproton1D_baseAsym.fillHist(h,h2);
			pi0eta1D_RectSBSubRegion.fillHist(h,h2);
			pi0eta1D_RectSBSubRegion_fixed.fillHist(h,h2);
			pi0MassDiffSubDetectors.fillHist(h,h2);
			etaMassDiffSubDetectors.fillHist(h,h2);
			if ( strcmp(h->GetName(),"pi0proton1D_mMandelstamT_mdelta")==0 ){ h->SetTitle("all t'");}
			if ( strcmp(h->GetName(),"pi0proton1D_mDelta")==0 ){ h->SetTitle("|t'|<1GeV^2");}
			if ( strcmp(h->GetName(),"pi0eta1D_mMandelstamT")==0 ){ h->SetTitle("all t'");}
			if ( strcmp(h->GetName(),"pi0eta1D_Cut")==0 ){ h->SetTitle("|t'|<1GeV^2");}
			pi0eta1DtAlltCut.fillHist(h,h2);
			pi0proton1D_beforeAfterT.fillHist(h,h2);
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
	pi0proton1D_baseAsym.plot("newGraphs/pi0proton1D_baseAsym.png",false,lineCutThresholds);
	etaproton1D_baseAsym.plot("newGraphs/etaproton1D_baseAsym.png",false,lineCutThresholds);

	lineCutThresholds={selectPi0Proton};
	pi0proton1D_mMandelstamT_mdelta.plot("newGraphs/pi0proton1D_mMandelstamT_mdelta_showCut.png",true,lineCutThresholds);
	lineCutThresholds={selectEtaProton};
	etaproton1D_mMandelstamT_mdelta.plot("newGraphs/etaproton1D_mMandelstamT_mdelta_showCut.png",true,lineCutThresholds);
	
	//cutThreshold2D = { {0.134,0.538,0.013,0.04 }, {0.134,0.538,0.0155,0.05 }, {0.134,0.538, 0.0205,0.07} };
	//pi0eta_Meas_mEllipsePre_showEllipse.plot("newGraphs/pi0eta_Meas_mEllipsePre_showEllipse.png","ellipse",cutThreshold2D);
	
	cutThreshold2D = { {0.135784, 0.548036, 2*0.00753584, 2*0.0170809 } }; // kinFit
	//cutThreshold2D = { {0.13381, 0.5388, 3*0.006, 3*0.0264 } };//eta3pi
	pi0eta_mEllipsePre.plot("newGraphs/pi0eta_mEllipsePre_withCut.png","ellipse",cutThreshold2D);
}
