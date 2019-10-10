void countEventsInBins(){
	TString nameAffix = "tLT1";
	TString treeName = "kin";
	//TString fileTag = "data";
	for (int i=0; i<60;++i){
	//	cout << "divideGEANT4/bin_"+std::to_string(i)+"/"+fileTag+"_pi0eta_"+nameAffix+"_"+std::to_string(i)+".root" << endl;
		TFile *myFile = TFile::Open("divideRoot/bin_"+std::to_string(i)+"/"+"data_pi0eta_"+nameAffix+"_"+std::to_string(i)+".root");
		TTree *t = (TTree *)myFile->Get(treeName);
		cout << "Data: Entries in bin_"+std::to_string(i)+": "+t->GetEntries() << endl;;
		myFile = TFile::Open("divideRoot/bin_"+std::to_string(i)+"/"+"reco_pi0eta_"+nameAffix+"_"+std::to_string(i)+".root");
		t = (TTree *)myFile->Get(treeName);
		cout << "Reco: Entries in bin_"+std::to_string(i)+": "+t->GetEntries() << endl;;
		myFile = TFile::Open("divideRoot/bin_"+std::to_string(i)+"/"+"thrown_pi0eta_"+nameAffix+"_"+std::to_string(i)+".root");
		t = (TTree *)myFile->Get(treeName);
		cout << "Gen: Entries in bin_"+std::to_string(i)+": "+t->GetEntries() << endl;;
		cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	}
	
}
