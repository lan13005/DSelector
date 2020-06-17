// Was used as a check to see if the isNotRepeated uniqueness tracking scheme was doing what it should be
// I think the result of this was yes

void graphCheckUniqueContrib(){	
        TFile* dataFile=new TFile("/d/grid15/ln16/pi0eta/092419/pi0eta_a0_recotreeFlat_DSelector.root");
	TTree *dataTree;
	dataFile->GetObject("pi0eta_a0_recotree_flat",dataTree);
    	TCanvas *allCanvases = new TCanvas("anyHists","",1440,900);
	THStack* stackedHists = new THStack("stackedHists","");

	TH1F* Meta_tot = new TH1F( "Meta", "M(#eta) ; M(#eta) GeV; Events/0.002 GeV", 300, 0.25, 0.85 );
	TH1F* Meta_bkg = new TH1F( "Meta_bkg", "M(#eta) ; M(#eta) GeV; Events/0.002 GeV", 300, 0.25, 0.85 );
	TH1F* Mpi0_tot = new TH1F( "Mpi0", "M(#pi_{0}) ; M(#pi_{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
	TH1F* Mpi0_bkg = new TH1F( "Mpi0_bkg", "M(#pi_{0}) ; M(#pi_{0}) GeV; Events/0.001 GeV", 200, 0.05, 0.25 );
	TH1F* Mpi0eta_tot = new TH1F( "Mpi0eta", "M(#pi_{0}#eta) ; M(#pi_{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );
	TH1F* Mpi0eta_bkg = new TH1F( "Mpi0eta_bkg", "M(#pi_{0}#eta) ; M(#pi_{0}#eta) GeV; Events/0.01 GeV", 350, 0, 3.5 );

	gStyle->SetOptFit(111);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);

	double Meta;
	double Mpi0;
	double Mpi0eta;
        bool isUniqueEtaB;
        bool isUniquePi0B;
        bool isUniquePi0EtaB;
	dataTree->SetBranchAddress("Meta",&Meta);
	dataTree->SetBranchAddress("Mpi0",&Mpi0);
	dataTree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        dataTree->SetBranchAddress("isNotRepeated_eta",&isUniqueEtaB);
        dataTree->SetBranchAddress("isNotRepeated_pi0",&isUniquePi0B);
        dataTree->SetBranchAddress("isNotRepeated_pi0eta",&isUniquePi0EtaB);
	Long64_t nentries=dataTree->GetEntries();
	for (int ientry=0; ientry<nentries; ientry++)
	{
	    dataTree->GetEntry(ientry);
            //cout << isUniqueEtaB << ", " << isUniquePi0B << "," << isUniquePi0EtaB << endl;

            if(!isUniqueEtaB){
               Meta_bkg->Fill(Meta); 
            }
            Meta_tot->Fill(Meta);
            if(!isUniquePi0B){
               Mpi0_bkg->Fill(Mpi0); 
            }
            Mpi0_tot->Fill(Mpi0);
            if(!isUniquePi0EtaB){
               Mpi0eta_bkg->Fill(Mpi0eta); 
            }
            Mpi0eta_tot->Fill(Mpi0eta);
        }

        //Meta_tot->SetFillColorAlpha(kBlue,0.3);
        Meta_bkg->SetFillColorAlpha(kMagenta,0.3);
        stackedHists->Add(Meta_tot,"HIST");
        stackedHists->Add(Meta_bkg,"HIST");
        stackedHists->Draw("nostack");
	stackedHists->GetXaxis()->SetTitle(Meta_tot->GetXaxis()->GetTitle());
	stackedHists->SetTitle("M(#eta)");
        allCanvases->SaveAs("Meta.png");

        allCanvases->Clear();
        stackedHists = new THStack("stackedHists","");
        //Mpi0_tot->SetFillColorAlpha(kBlue,0.3);
        Mpi0_bkg->SetFillColorAlpha(kMagenta,0.3);
        stackedHists->Add(Mpi0_tot,"HIST");
        stackedHists->Add(Mpi0_bkg,"HIST");
        stackedHists->Draw("nostack");
	stackedHists->GetXaxis()->SetTitle(Mpi0_tot->GetXaxis()->GetTitle());
	stackedHists->SetTitle("M(#pi0)");
        allCanvases->SaveAs("Mpi0.png");

        allCanvases->Clear();
        stackedHists = new THStack("stackedHists","");
        //Mpi0eta_tot->SetFillColorAlpha(kBlue,0.3);
        Mpi0eta_bkg->SetFillColorAlpha(kMagenta,0.3);
        stackedHists->Add(Mpi0eta_tot,"HIST");
        stackedHists->Add(Mpi0eta_bkg,"HIST");
        stackedHists->Draw("nostack");
	stackedHists->GetXaxis()->SetTitle(Mpi0eta_tot->GetXaxis()->GetTitle());
	stackedHists->SetTitle("M(#pi0eta)");
        allCanvases->SaveAs("Mpi0eta.png");
}
