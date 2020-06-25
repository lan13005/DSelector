// Apparently it is interesting to look at the asymmetry between the fast eta and the fast pion in the double 
// regge region. Here we are loading the histograms and making the asymmetry. The input histograms have different
// upper thresholds for teta or tpi which is set at 1 down to 0.25 in 0.25 increments
//

void drawHistAndAcc(TH1F* dataHist, TH1F* accHist, TCanvas* c1){
    dataHist->SetMinimum(0);
    dataHist->Draw();
    c1->Update(); // need to update otherwise the gPad will not grab the correct value
    //scale accHist to the pad coordinates
    Float_t rightmax = 1.1*accHist->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    accHist->SetLineColor(kRed);
    accHist->Scale(scale);

    //draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    accHist->Draw("SAME");
    axis->Draw("SAME");
}

void asymmetryFastEtaPi0(){
    gStyle->SetOptStat(0);

    TFile* dataFile = new TFile("degALL_data_2017_mEllipsePre_hists_DSelector.root");
    TFile* genFile = new TFile("degALL_gen_2017_hists_DSelector.root");
    TFile* accFile = new TFile("degALL_flat_2017_mEllipsePre_hists_DSelector.root");

    TH1F* dHist_fastEta_data;
    TH1F* dHist_fastPi0_data;
    TH1F* dHist_fastEta_acc;
    TH1F* dHist_fastPi0_acc;
    TH1F* dHist_Mpi0eta_gen;

    string baseFastEtaString = "pi0eta1D_pFastEtaThresh";
    string baseFastPi0String = "pi0eta1D_pFastPi0Thresh";
    TCanvas* allCanvases = new TCanvas("","",1440,900);
    for (int it=0; it<4; ++ it){
        allCanvases->Clear();
        allCanvases->Divide(3,2);

        cout << "Getting histograms with fast t threshold " << it  << endl;
        dataFile->GetObject((baseFastEtaString+to_string(it)).c_str(),dHist_fastEta_data);
        dataFile->GetObject((baseFastPi0String+to_string(it)).c_str(),dHist_fastPi0_data);
        accFile->GetObject((baseFastEtaString+to_string(it)).c_str(),dHist_fastEta_acc);
        accFile->GetObject((baseFastPi0String+to_string(it)).c_str(),dHist_fastPi0_acc);
        genFile->GetObject("Mpi0eta",dHist_Mpi0eta_gen);

        dHist_fastEta_data->SetLineColor(kBlue-3);
        dHist_fastPi0_data->SetLineColor(kBlue-3);

        // get the acceptances
        TH1F* dHist_fastEta_acceptance = (TH1F*)dHist_fastEta_acc->Clone();
        TH1F* dHist_fastPi0_acceptance = (TH1F*)dHist_fastPi0_acc->Clone();
        dHist_fastEta_acceptance->Divide(dHist_Mpi0eta_gen);
        dHist_fastPi0_acceptance->Divide(dHist_Mpi0eta_gen);
        // need to save a copy of  the acceptances since when drawing the histogram on a double axis we have to scale
        TH1F* dHist_fastEta_acceptance_unscaled = (TH1F*)dHist_fastEta_acceptance->Clone();
        TH1F* dHist_fastPi0_acceptance_unscaled = (TH1F*)dHist_fastPi0_acceptance->Clone();

        // Draw the data and acceptance histograms overlaid 
        allCanvases->cd(1);
        dHist_fastEta_data->SetTitle("M(#pi#eta) yield - fast #eta");
        drawHistAndAcc(dHist_fastEta_data,dHist_fastEta_acceptance, allCanvases);
        allCanvases->cd(2);
        dHist_fastPi0_data->SetTitle("M(#pi#eta) yield - fast #pi");
        drawHistAndAcc(dHist_fastPi0_data,dHist_fastPi0_acceptance, allCanvases);

        // Acceptance correct
        TH1F* dHist_fastEta_corrected = (TH1F*)dHist_fastEta_data->Clone();
        TH1F* dHist_fastPi0_corrected = (TH1F*)dHist_fastPi0_data->Clone();
        dHist_fastEta_corrected->SetTitle("M(#pi#eta) Yield fast #eta  acceptance corrected");
        dHist_fastPi0_corrected->SetTitle("M(#pi#eta) Yield fast #pi acceptance corrected");
        allCanvases->cd(3);
        dHist_fastEta_corrected->Divide(dHist_fastEta_acceptance_unscaled);
        dHist_fastEta_corrected->Draw();
        allCanvases->cd(4);
        dHist_fastPi0_corrected->Divide(dHist_fastPi0_acceptance_unscaled);
        dHist_fastPi0_corrected->Draw();

        allCanvases->cd(5); 
        cout << "Drawing the asymmetries" << endl;
        TH1F* dHist_asymFastEtaFastPi = (TH1F*)dHist_fastEta_corrected->GetAsymmetry(dHist_fastPi0_corrected);
        dHist_asymFastEtaFastPi->SetTitle("#frac{N_{fast#eta}-N_{fast#pi}}{N_{fast#eta}+N_{fast#pi}}");
        dHist_asymFastEtaFastPi->Draw();

        allCanvases->cd(6);
	TPaveText *histID_pt = new TPaveText(0.3,0.3,0.6,0.6);
        double threshold = 0.25*(1+it);
	histID_pt->AddText("By Fast I Mean:");
        histID_pt->AddText(("t_{#eta} / t_{#pi} < "+std::to_string(threshold)).c_str());
	histID_pt->Draw();

        //allCanvases->SaveAs(((string)"asymmetryNfastEtafastPi/asymmetry"+to_string(it)+".png").c_str());
        if(it==0){
            allCanvases->Print("asymmetryNfastEtafastPi/asymmetry.pdf(","pdf");
        }
        else if (it==3){
            allCanvases->Print("asymmetryNfastEtafastPi/asymmetry.pdf)","pdf");
        }
        else{
            allCanvases->Print("asymmetryNfastEtafastPi/asymmetry.pdf","pdf");
        }
    }
}
