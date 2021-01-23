void countEventsInBins(){
	TString nameAffix = "amptools";
	TString treeName = "kin";
	//TString fileTag = "data";
	for (int i=0; i<60;++i){
	//	cout << "divideGEANT4/bin_"+std::to_string(i)+"/"+fileTag+"_pi0eta_"+nameAffix+"_"+std::to_string(i)+".root" << endl;
		TFile *myFile = TFile::Open(("EtaPi_fit_onlypi1_gen/bin_"+std::to_string(i)+"/amptools_onlypi1_gen_"+std::to_string(i)+".root").c_str());
		TTree *t = (TTree *)myFile->Get(treeName);
		cout << "Data: Entries in bin_"+std::to_string(i)+": "+t->GetEntries() << endl;;
		myFile = TFile::Open(("EtaPi_fit_onlypi1_gen/bin_"+std::to_string(i)+"/amptools_gen_"+std::to_string(i)+".root").c_str());
		t = (TTree *)myFile->Get(treeName);
		cout << "Reco: Entries in bin_"+std::to_string(i)+": "+t->GetEntries() << endl;;
		myFile = TFile::Open(("EtaPi_fit_onlypi1_gen/bin_"+std::to_string(i)+"/amptools_gen_"+std::to_string(i)+".root").c_str());
		t = (TTree *)myFile->Get(treeName);
		cout << "Gen: Entries in bin_"+std::to_string(i)+": "+t->GetEntries() << endl;;
		cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	}
	
}
