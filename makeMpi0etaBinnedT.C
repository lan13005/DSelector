#include "/d/grid13/gluex/gluex_top/gluex_style.C"

void makeMpi0etaBinnedT(string fileLoc){
        gluex_style();
	TFile* dataFile = TFile::Open(fileLoc.c_str());
	TCanvas *allCanvases = new TCanvas("","",1440,900);
	gPad->SetLeftMargin(0.2);
	gPad->SetBottomMargin(0.2);
	allCanvases->Divide(3,3,0,0);
	gStyle->SetOptStat(0);


	TH1F *anyHist;
	string baseName = "pi0eta1D_Cut_tBin";
	string masses[10] = {"0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"};

	//gStyle->SetTitleSize(0.08,"t");


	double maxY = DBL_MIN;
	double minY = DBL_MAX;
	for (int iBin=1; iBin<10; ++iBin){
		dataFile->GetObject((baseName+to_string(iBin)).c_str(), anyHist);
		if (maxY < anyHist->GetMaximum() ){
			maxY = anyHist->GetMaximum();
		}
		if (minY > anyHist->GetMinimum() ){
			minY = anyHist->GetMinimum();
		}
	}


	for (int iBin=1; iBin<10; ++iBin){
		allCanvases->cd(iBin);
		dataFile->GetObject((baseName+to_string(iBin)).c_str(), anyHist);
		anyHist->SetAxisRange(minY, maxY*1.05,"Y");
		anyHist->SetTitle((masses[iBin-1]+" < |t| < "+masses[iBin]).c_str() );
		//anyHist->GetXaxis()->SetLabelSize(0.075);
		//anyHist->GetYaxis()->SetLabelSize(0.075);
		//anyHist->GetXaxis()->SetTitleSize(0.075);
		//anyHist->GetYaxis()->SetTitleSize(0.075);
		anyHist->Draw("HIST");
		anyHist->GetYaxis()->SetTitleSize(0.05);
		anyHist->GetXaxis()->SetTitleSize(0.05);
                gStyle->SetTitleSize(0.1,"t");
	}

	allCanvases->SaveAs("newGraphs/Mpi0etaBinnedT.png");
}
