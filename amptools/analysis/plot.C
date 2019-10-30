#include "plot.h"

void plot(){
	TCanvas *canvas = new TCanvas("canvas","canvas",1440, 900);
	gStyle->SetOptStat(0);
	canvas->Divide(3,3,0,0);
	int numBins=65;
	int numAmps=8; // 7 amplitudes with 1 total amplitude


	loadAmplitudes* ampData = new loadAmplitudes("amplitudes.txt", numAmps, numBins);
	ampData->load();

	TH1F* amplitude;
	std::vector<string> ampNames = { "S0m", "P0m", "P1m", "D0m", "D1m", "P1p", "D1p" }; 

	amplitude = new TH1F("","",numBins,0.7,2.0);	
	cout << "Plotting total" << endl;
	double maxAmp = 0;
	canvas->cd(1);
	for (Int_t i=0; i<numBins; ++i){
		amplitude->SetBinContent( i, ampData->getAmp("total")[2*i] );
		amplitude->SetBinError( i, 0);//ampData->getAmp("total")[2*i+1] );
	}
	amplitude->SetTitle("total");
	amplitude->Draw("E");
	maxAmp = amplitude->GetMaximum();
	cout << "max amplitude value in total summed amplitudes is " << maxAmp << " will use to set y axis of future histograms" << endl;

	for (Int_t iAmp=0; iAmp<ampNames.size(); ++iAmp){
		amplitude = new TH1F("","",numBins,0.7,2.0);	
		amplitude->SetAxisRange(0,maxAmp,"Y");
		cout << "Plotting " << ampNames[iAmp] << endl;
		canvas->cd(iAmp+2);	
		for (Int_t i=0; i<numBins; ++i){
			amplitude->SetBinContent( i, ampData->getAmp(ampNames[iAmp])[2*i] );
			amplitude->SetBinError( i, ampData->getAmp(ampNames[iAmp])[2*i+1] );
		}
		amplitude->SetTitle(ampNames[iAmp].c_str());
		amplitude->Draw("E");

	}
	canvas->SaveAs("amplitudes.pdf");
}
