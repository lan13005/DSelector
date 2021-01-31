//#include "/d/grid13/gluex/gluex_top/gluex_style.C"
string fileType="pdf";

void makeSpecialtyPlots(string fileLoc){
    //gluex_style();
    TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
    TFile* data2017 = TFile::Open(fileLoc.c_str());
    TH1F *data2017_hist;
    gStyle->SetOptStat(0);
    allCanvases->SetLogy(1);
    data2017->GetObject("mandelstam_tp_mMandelstamT",data2017_hist);
    data2017_hist->SetTitle("");
    data2017_hist->Draw("HIST");
    allCanvases->SaveAs(("newGraphs/mandelstam_tp_log."+fileType).c_str());
    allCanvases->SetLogy(0);
}
