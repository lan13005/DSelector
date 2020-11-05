#include "/d/grid13/gluex/gluex_top/gluex_style.C"
string fileType="pdf";

void makeAdd20172018(){
    gSystem->mkdir("newGraphs/Mpi0eta_phase1");
    
    gluex_style();
    TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
    TFile* data2017 = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_data_2017_mEllipsePre_hists_DSelector.root");
    TFile* data2018_1 = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_data_2018_1_mEllipsePre_hists_DSelector.root");
    TFile* data2018_8 = TFile::Open("/d/grid15/ln16/pi0eta/092419/degALL_data_2018_8_mEllipsePre_hists_DSelector.root");
    TH1F *data2017_hist;
    TH1F *data2018_1_hist;
    TH1F *data2018_8_hist;

    gStyle->SetOptStat(0);

    data2017->GetObject("pi0eta1D_mMandelstamT",data2017_hist);
    data2018_1->GetObject("pi0eta1D_mMandelstamT",data2018_1_hist);
    data2018_8->GetObject("pi0eta1D_mMandelstamT",data2018_8_hist);
    TH1F* pi0eta_tAll = (TH1F*)data2017_hist->Clone();
    cout << "nentries 2017: " << data2017_hist->GetEntries() << endl;
    cout << "nentries 2018_1: " << data2018_1_hist->GetEntries() << endl;
    cout << "nentries 2018_8: " << data2018_8_hist->GetEntries() << endl;
    cout << "total entries: " << data2017_hist->GetEntries()+data2018_1_hist->GetEntries()+data2018_8_hist->GetEntries() << endl;
    pi0eta_tAll->Add(data2018_1_hist);
    pi0eta_tAll->Add(data2018_8_hist);
    int integral = (int)pi0eta_tAll->Integral();
    cout << "Integral: " << integral << endl;
    pi0eta_tAll->SetLineColor(kMagenta);
    pi0eta_tAll->Draw("HIST");
    pi0eta_tAll->SetTitle(("Number of Events: "+to_string(integral)).c_str());

    data2017->GetObject("pi0eta1D_Cut",data2017_hist);
    data2018_1->GetObject("pi0eta1D_Cut",data2018_1_hist);
    data2018_8->GetObject("pi0eta1D_Cut",data2018_8_hist);
    TH1F* pi0eta_tLT1 = (TH1F*)data2017_hist->Clone();
    pi0eta_tLT1->Add(data2018_1_hist);
    pi0eta_tLT1->Add(data2018_8_hist);
    pi0eta_tLT1->SetLineColor(kBlue);
    pi0eta_tLT1->Draw("SAME HIST");

    allCanvases->SaveAs(("newGraphs/Mpi0eta_phase1/Mpi0eta."+fileType).c_str());

    allCanvases->Clear();
    pi0eta_tLT1->Draw("HIST");
    pi0eta_tLT1->SetTitle("");
    allCanvases->SaveAs(("newGraphs/Mpi0eta_phase1/Mpi0eta_tLT1."+fileType).c_str());

}
