void getMacceptance(){
        gStyle->SetOptStat(kFALSE);
        //CutCLOrderNeg
        TFile* infile_acc = TFile::Open("v20_pi0eta_flat_hists_DSelector_pi0eta.root");
        TFile* infile_gen = TFile::Open("v20_pi0eta_gen_hists_DSelector_pi0eta.root");
        TCanvas *allCanvases = new TCanvas("anyHists","",1920,1080);
        std::vector<std::string> names1D_acc;
        std::vector<std::string> names1D_gen;

        // 1D hists
        names1D_acc.push_back("pi0eta1DNoCut_AccSubEBeam_AccSubBaryBkg1");
        names1D_gen.push_back("pi0eta1D8GeVPlus");

        int namesSize1D = static_cast<int>(names1D_acc.size());

        TH1F *any1DHist_acc;
        TH1F *any1DHist_gen;
        for (int histIdx=0; histIdx<namesSize1D; ++histIdx){
                allCanvases->Clear();
                infile_acc->GetObject(names1D_acc[histIdx].c_str(),any1DHist_acc);
                infile_gen->GetObject(names1D_gen[histIdx].c_str(),any1DHist_gen);
                
                any1DHist_acc->SetTitle("acc-Left | gen-Right");
                any1DHist_acc->Draw();
                allCanvases->Update();

                //scale hint1 to the pad coordinates
                Float_t rightmax = 1.1*any1DHist_gen->GetMaximum();
                Float_t scale    = gPad->GetUymax()/rightmax;
                any1DHist_gen->SetLineColor(kRed);
                any1DHist_gen->Scale(scale);
                any1DHist_gen->Draw("HISTSAME");
                //draw an axis on the right side
                TGaxis*axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                                         gPad->GetUxmax(),gPad->GetUymax(),
                                         0,rightmax,510,"+L");
                
                axis->SetLineColor(kRed);
                axis->SetLabelColor(kRed);
                axis->Draw();

                allCanvases->SaveAs("mass_superimposed.png");
                allCanvases->Clear();

		any1DHist_acc->Divide(any1DHist_gen);	
                any1DHist_acc->SetTitle("Acceptance");
                any1DHist_acc->Draw("HIST");
                allCanvases->SaveAs("mass_acceptance.png");
        }
}

