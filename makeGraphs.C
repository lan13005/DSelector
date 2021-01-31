#include "makeGraphs.h"
#include "/d/grid15/ln16/pi0eta/q-values/makeDiagnosticHists.h"
#include "DSelector_helperFuncs.h"

//#include "/d/grid13/gluex/gluex_top/gluex_style.C"

double selectPi0Proton=1.4;
double selectEtaProton=2;

class overlayPlots{
	private:
		std::vector<TH1F*> overlayHists;
		std::vector<double> histWeights;
		std::vector<int> histFillColors;
        	TLegend *leg = new TLegend(0.6,0.75,0.9,0.95);
		int colors[10] = {4, 6, kOrange, 7, 9, 30, 27, 46, 41};
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
                        //gluex_style();
			cout << "Creating overlay hist " << overlayHists[0]->GetName() << endl;
			TCanvas* overlayCanvas = new TCanvas( ("canvas_"+fileName).c_str(),"",1440,900);
			overlayHists[0]->Scale(histWeights[0]);
			overlayHists[0]->SetLineColor( colors[0] );
			overlayHists[0]->SetLineWidth( 3 ) ;
			leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
			//overlayHists[0]->SetTitle("");
			overlayHists[0]->Draw("HIST");
			//overlayHists[0]->GetXaxis()->SetTitleSize(0.05);
			//overlayHists[0]->GetYaxis()->SetTitleSize(0.05);
			overlayHists[0]->SetAxisRange(minimum1D,maximum1D*1.05,"Y");
			for (int i=1; i<overlayHists.size(); ++i){
				cout << "Overlaying hist " << overlayHists[i]->GetName() << endl;
				leg->AddEntry(overlayHists[i],overlayHists[i]->GetTitle() , "l");
				overlayHists[i]->Scale(histWeights[i]);
				overlayHists[i]->SetLineColor( colors[i] );
				overlayHists[i]->SetLineWidth( 3 ) ;
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
                        //gluex_style();
			cout << "Creating side by side hist " << overlayHists[0]->GetName() << endl;
			TCanvas* overlayCanvas = new TCanvas( ("canvas_"+fileName).c_str(),"",1440,900);
			overlayHists[0]->Scale(histWeights[0]);
			//leg->AddEntry(overlayHists[0],overlayHists[0]->GetTitle() , "l");
			overlayHists[0]->Draw("COLZ");
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



void makeGraphs(string fileLoc){
        //gluex_style();


        TFile* file = TFile::Open(fileLoc.c_str());
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
        trackHists = { {"pi0Mass_Kin_mEllipsePre",  1 }, {"pi0MassFCAL_Kin_mEllipsePre", 1 }, {"pi0MassBCAL_Kin_mEllipsePre", 1 }, {"pi0MassSPLIT_Kin_mEllipsePre", 1 } };
	overlayPlots pi0MassDiffSubDetectors(trackHists);
        trackHists = { {"etaMass_Kin_mEllipsePre",  1 }, {"etaMassFCAL_Kin_mEllipsePre", 1 }, {"etaMassBCAL_Kin_mEllipsePre", 1 }, {"etaMassSPLIT_Kin_mEllipsePre", 1 } };
	overlayPlots etaMassDiffSubDetectors(trackHists);
	trackHists = { {"pi0eta1D_mMandelstamT", 1}, {"pi0eta1D_Cut", 1} };
	overlayPlots pi0eta1DtAlltCut(trackHists);
	trackHists = { {"pi0proton1D_mMandelstamT_mdelta", 1}, {"pi0proton1D_mDelta", 1} };
	overlayPlots pi0proton1D_beforeAfterT(trackHists);
        trackHists = { {"pi0Mass_Kin_mEllipsePre",  1 } };
        overlayPlots pi0Mass_sbRegions(trackHists);
        trackHists = { {"etaMass_Kin_mEllipsePre",  1 } };
        overlayPlots etaMass_sbRegions(trackHists);
	//trackHists = { {"pi0proton1D_baseAsymCut", 1}, {"pi0proton1D_baseAsymCut_mDelta_fastEta", 1} };
	//overlayPlots pi0proton1D_baseAsym(trackHists);
	//trackHists = { {"etaproton1D_baseAsymCut", 1}, {"etaproton1D_baseAsymCut_mDelta_fastPi0", 1} };
	//overlayPlots etaproton1D_baseAsym(trackHists);

	//trackHists = { {"pi0eta1D_baseAsymCut_res", 1}, {"pi0eta1D_pFastEtaBin0", 1}, {"pi0eta1D_pFastEtaBin1", 1}, {"pi0eta1D_pFastEtaBin2", 1}, {"pi0eta1D_pFastEtaBin3", 1}, {"pi0eta1D_pFastEtaBin4", 1} , {"pi0eta1D_pFastEtaBin5", 1} };
	//overlayPlots pi0eta1D_pFastEta(trackHists);
	//trackHists = { {"pi0eta1D_baseAsymCut_res", 1}, {"pi0eta1D_pFastPi0Bin0", 1}, {"pi0eta1D_pFastPi0Bin1", 1}, {"pi0eta1D_pFastPi0Bin2", 1}, {"pi0eta1D_pFastPi0Bin3", 1}, {"pi0eta1D_pFastPi0Bin4", 1}, {"pi0eta1D_pFastPi0Bin5", 1} };
	//overlayPlots pi0eta1D_pFastPi0(trackHists);
	//trackHists = { {"pi0eta1D_pFastEtaBin0", 1} };
	//overlayPlots pi0eta1D_pFastEtaBin0(trackHists);
	//trackHists = { {"pi0eta1D_pFastEtaBin1", 1} };
	//overlayPlots pi0eta1D_pFastEtaBin1(trackHists);
	//trackHists = { {"pi0eta1D_pFastEtaBin2", 1} };
	//overlayPlots pi0eta1D_pFastEtaBin2(trackHists);
	//trackHists = { {"pi0eta1D_pFastEtaBin3", 1} };
	//overlayPlots pi0eta1D_pFastEtaBin3(trackHists);
	//trackHists = { {"pi0eta1D_pFastEtaBin4", 1} };
	//overlayPlots pi0eta1D_pFastEtaBin4(trackHists);
	//trackHists = { {"pi0eta1D_pFastEtaBin5", 1} };
	//overlayPlots pi0eta1D_pFastEtaBin5(trackHists);
	//trackHists = { {"pi0eta1D_pFastPi0Bin0", 1} };
	//overlayPlots pi0eta1D_pFastPi0Bin0(trackHists);
	//trackHists = { {"pi0eta1D_pFastPi0Bin1", 1} };
	//overlayPlots pi0eta1D_pFastPi0Bin1(trackHists);
	//trackHists = { {"pi0eta1D_pFastPi0Bin2", 1} };
	//overlayPlots pi0eta1D_pFastPi0Bin2(trackHists);
	//trackHists = { {"pi0eta1D_pFastPi0Bin3", 1} };
	//overlayPlots pi0eta1D_pFastPi0Bin3(trackHists);
	//trackHists = { {"pi0eta1D_pFastPi0Bin4", 1} };
	//overlayPlots pi0eta1D_pFastPi0Bin4(trackHists);
	//trackHists = { {"pi0eta1D_pFastPi0Bin5", 1} };
	//overlayPlots pi0eta1D_pFastPi0Bin5(trackHists);
	//	

        //trackHists = { {"pi0proton1D_baseAsymCut_mDelta_fastEta", 1} };
        //overlayPlots pi0proton1D_baseAsymCut_mDelta_fastEta(trackHists);
        //trackHists = { {"etaproton1D_baseAsymCut_mDelta_fastPi0", 1} };
        //overlayPlots etaproton1D_baseAsymCut_mDelta_fastPi0(trackHists);

        //trackHists = { {"//labThetaSPLIT_eta_mEllipsePre", 1}, {"labThetaFCAL_eta_mEllipsePre", 1}, {"labThetaBCAL_eta_mEllipsePre", 1} }; 
	//overlayPlots labTheta_eta(trackHists);
        //trackHists = { {"//labThetaSPLIT_pi0_mEllipsePre", 1}, {"labThetaFCAL_pi0_mEllipsePre", 1}, {"labThetaBCAL_pi0_mEllipsePre", 1} }; 
	//overlayPlots labTheta_pi0(trackHists);

	TCanvas *c1 = new TCanvas("c1","",1440,900);
	int i=0;
   	while ((key = (TKey*)keyList())) {
   	   	TClass *cl = gROOT->GetClass(key->GetClassName());
   	   	if (cl->InheritsFrom("TH2")){
			string fileName="newGraphs/";
   	   		TH2F *h = (TH2F*)key->ReadObj();
			//h->GetXaxis()->SetTitleSize(0.04);
			//h->GetYaxis()->SetTitleSize(0.04);
   	   		h->Draw("COLZ HIST");
			fileName.append(h->GetName());
			fileName.append(".png");
                        gluex_style->cd();
                        gROOT->ForceStyle();
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
		}
		else if (cl->InheritsFrom("TH1")){
			string fileName="newGraphs/";
   	   		TH1F *h = (TH1F*)key->ReadObj();
			if (strncmp(h->GetName(),"prodPlanePSphi", strlen("prodPlanePSphi"))==0 ) { continue; } 
			//h->GetXaxis()->SetTitleSize(0.04);
			//h->GetYaxis()->SetTitleSize(0.04);
   	   		h->Draw("COLZ HIST");
			fileName.append(h->GetName());
			fileName.append(".png");
                        gluex_style->cd();
                        gROOT->ForceStyle();
   	   		c1->SaveAs((fileName).c_str());
			// Wait until we have finally used up TH1 object first. Otherwise casting it into TH1F early creates some problems
			pi0proton1D_mMandelstamT_mdelta.fillHist(h);
			etaproton1D_mMandelstamT_mdelta.fillHist(h);
			//pi0proton1D_baseAsym.fillHist(h);
			//etaproton1D_baseAsym.fillHist(h);
			pi0eta1D_RectSBSubRegion.fillHist(h);
			pi0eta1D_RectSBSubRegion_fixed.fillHist(h);
			pi0MassDiffSubDetectors.fillHist(h);
			etaMassDiffSubDetectors.fillHist(h);
			//pi0eta1D_pFastEtaBin0.fillHist(h); 
			//pi0eta1D_pFastEtaBin1.fillHist(h); 
			//pi0eta1D_pFastEtaBin2.fillHist(h); 
			//pi0eta1D_pFastEtaBin3.fillHist(h); 
			//pi0eta1D_pFastEtaBin4.fillHist(h); 
			//pi0eta1D_pFastEtaBin5.fillHist(h); 
			//pi0eta1D_pFastPi0Bin0.fillHist(h); 
			//pi0eta1D_pFastPi0Bin1.fillHist(h); 
			//pi0eta1D_pFastPi0Bin2.fillHist(h); 
			//pi0eta1D_pFastPi0Bin3.fillHist(h); 
			//pi0eta1D_pFastPi0Bin4.fillHist(h); 
			//pi0eta1D_pFastPi0Bin5.fillHist(h); 
			//pi0eta1D_pFastEta.fillHist(h);
			//pi0eta1D_pFastPi0.fillHist(h);
                        //etaproton1D_baseAsymCut_mDelta_fastPi0.fillHist(h);
                        //pi0proton1D_baseAsymCut_mDelta_fastEta.fillHist(h);
                        if ( strcmp(h->GetName(),"pi0proton1D_baseAsymCut_mDelta_fastEta")==0 ){ h->SetTitle("M(pi0proton) fastEta"); }
                        if ( strcmp(h->GetName(),"etaproton1D_baseAsymCut_mDelta_fastPi0")==0 ){ h->SetTitle("M(etaproton) fastPi0"); }
			if ( strcmp(h->GetName(),"pi0proton1D_mMandelstamT_mdelta")==0 ){ h->SetTitle("all t'");}
			if ( strcmp(h->GetName(),"pi0proton1D_mDelta")==0 ){ h->SetTitle("|t'|<1GeV^2");}
			if ( strcmp(h->GetName(),"pi0eta1D_mMandelstamT")==0 ){ h->SetTitle("all t'");}
			if ( strcmp(h->GetName(),"pi0eta1D_Cut")==0 ){ h->SetTitle("|t'|<1GeV^2");}
                        // title changes for the legend for Mpi0eta in different t1 bins
                        if ( strcmp(h->GetName(),"pi0eta1D_baseAsymCut_res")==0 ) { h->SetTitle("Mpi0eta"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastEtaBin0")==0 ) { h->SetTitle("Mpi0eta: 0.0<t_{#eta}<0.2"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastEtaBin1")==0 ) { h->SetTitle("Mpi0eta: 0.2<t_{#eta}<0.4"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastEtaBin2")==0 ) { h->SetTitle("Mpi0eta: 0.4<t_{#eta}<0.6"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastEtaBin3")==0 ) { h->SetTitle("Mpi0eta: 0.6<t_{#eta}<0.8"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastEtaBin4")==0 ) { h->SetTitle("Mpi0eta: 0.8<t_{#eta}<1"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastEtaBin5")==0 ) { h->SetTitle("Mpi0eta: 1<t_{#eta}"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastPi0Bin0")==0 ) { h->SetTitle("Mpi0eta: 0.0<t_{#pi}<0.2"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastPi0Bin1")==0 ) { h->SetTitle("Mpi0eta: 0.2<t_{#pi}<0.4"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastPi0Bin2")==0 ) { h->SetTitle("Mpi0eta: 0.4<t_{#pi}<0.6"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastPi0Bin3")==0 ) { h->SetTitle("Mpi0eta: 0.6<t_{#pi}<0.8"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastPi0Bin4")==0 ) { h->SetTitle("Mpi0eta: 0.8<t_{#pi}<1"); }
                        if ( strcmp(h->GetName(),"pi0eta1D_pFastPi0Bin5")==0 ) { h->SetTitle("Mpi0eta: 1<t_{#pi}"); }
			pi0eta1DtAlltCut.fillHist(h);
			pi0proton1D_beforeAfterT.fillHist(h);
                        //labTheta_pi0.fillHist(h);
                        //labTheta_eta.fillHist(h);
                        etaMass_sbRegions.fillHist(h);
                        pi0Mass_sbRegions.fillHist(h);
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
	//pi0proton1D_baseAsym.plot("newGraphs/pi0proton1D_baseAsym.png",false,lineCutThresholds);
	//etaproton1D_baseAsym.plot("newGraphs/etaproton1D_baseAsym.png",false,lineCutThresholds);


        //lineCutThresholds = {4,10,17};
        //labTheta_pi0.plot("newGraphs/labTheta_pi0.png",true,lineCutThresholds);
        //lineCutThresholds = {3,6,11};
        //labTheta_eta.plot("newGraphs/labTheta_eta.png",true,lineCutThresholds);


	//lineCutThresholds = {0.9, 1.060, 1.24, 1.4, 1.65, 1.9, 2.15, 2.4, 2.65, 2.9};
	//pi0eta1D_pFastEta.plot("newGraphs/pi0eta1D_pFastEta.png",true,lineCutThresholds);
	//pi0eta1D_pFastPi0.plot("newGraphs/pi0eta1D_pFastPi0.png",true,lineCutThresholds);
	//pi0eta1D_pFastEtaBin0.plot("newGraphs/pi0eta1D_pFastEtaBin0.png", true, lineCutThresholds); 
	//pi0eta1D_pFastEtaBin1.plot("newGraphs/pi0eta1D_pFastEtaBin1.png", true, lineCutThresholds); 
	//pi0eta1D_pFastEtaBin2.plot("newGraphs/pi0eta1D_pFastEtaBin2.png", true, lineCutThresholds); 
	//pi0eta1D_pFastEtaBin3.plot("newGraphs/pi0eta1D_pFastEtaBin3.png", true, lineCutThresholds); 
	//pi0eta1D_pFastEtaBin4.plot("newGraphs/pi0eta1D_pFastEtaBin4.png", true, lineCutThresholds); 
	//pi0eta1D_pFastEtaBin5.plot("newGraphs/pi0eta1D_pFastEtaBin5.png", true, lineCutThresholds); 
	//pi0eta1D_pFastPi0Bin0.plot("newGraphs/pi0eta1D_pFastPi0Bin0.png", true, lineCutThresholds); 
	//pi0eta1D_pFastPi0Bin1.plot("newGraphs/pi0eta1D_pFastPi0Bin1.png", true, lineCutThresholds); 
	//pi0eta1D_pFastPi0Bin2.plot("newGraphs/pi0eta1D_pFastPi0Bin2.png", true, lineCutThresholds); 
	//pi0eta1D_pFastPi0Bin3.plot("newGraphs/pi0eta1D_pFastPi0Bin3.png", true, lineCutThresholds); 
	//pi0eta1D_pFastPi0Bin4.plot("newGraphs/pi0eta1D_pFastPi0Bin4.png", true, lineCutThresholds); 
	//pi0eta1D_pFastPi0Bin5.plot("newGraphs/pi0eta1D_pFastPi0Bin5.png", true, lineCutThresholds); 

	lineCutThresholds={selectPi0Proton};
	pi0proton1D_mMandelstamT_mdelta.plot("newGraphs/pi0proton1D_mMandelstamT_mdelta_showCut.png",true,lineCutThresholds);
        //lineCutThresholds={selectEtaProton};
	//etaproton1D_mMandelstamT_mdelta.plot("newGraphs/etaproton1D_mMandelstamT_mdelta_showCut.png",true,lineCutThresholds);

        //lineCutThresholds={1.5,1.65,1.9,2.2,2.4,2.8};
        //etaproton1D_baseAsymCut_mDelta_fastPi0.plot("newGraphs/etaproton1D_baseAsymCut_mDelta_fastPi0_showBins.png",true,lineCutThresholds);
        //lineCutThresholds={1.15,1.4,1.6,1.75,1.95,2.4};
        //pi0proton1D_baseAsymCut_mDelta_fastEta.plot("newGraphs/pi0proton1D_baseAsymCut_mDelta_fastEta_showBins.png",true,lineCutThresholds);

	cutThreshold2D = { {0.135784, 0.548036, 2*0.00753584, 2*0.0170809 } }; // kinFit
	//cutThreshold2D = { {0.13381, 0.5388, 3*0.006, 3*0.0264 } };//eta3pi
	pi0eta_mEllipsePre.plot("newGraphs/pi0eta_mEllipsePre_withCut.png","ellipse",cutThreshold2D);

        lineCutThresholds = { {
            pi0Mean-pi0Std*(pi0Sig+pi0Skip+pi0SB),
            pi0Mean-pi0Std*(pi0Sig+pi0Skip),
            pi0Mean-pi0Std*(pi0Sig),
            pi0Mean+pi0Std*(pi0Sig),
            pi0Mean+pi0Std*(pi0Sig+pi0Skip),
            pi0Mean+pi0Std*(pi0Sig+pi0Skip+pi0SB),
        } };
        pi0Mass_sbRegions.plot("newGraphs/pi0Mass_sbRegions.png",true, lineCutThresholds);
        lineCutThresholds = { {
            etaMean-etaStd*(etaSig+etaSkip+etaSB),
            etaMean-etaStd*(etaSig+etaSkip),
            etaMean-etaStd*(etaSig),
            etaMean+etaStd*(etaSig),
            etaMean+etaStd*(etaSig+etaSkip),
            etaMean+etaStd*(etaSig+etaSkip+etaSB),
        } };
        etaMass_sbRegions.plot("newGraphs/etaMass_sbRegions.png",true, lineCutThresholds);



}
