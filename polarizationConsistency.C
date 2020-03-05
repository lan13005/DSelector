void polarizationConsistency(){
	string beamAngles[6] = {"deg000","deg045","deg090","deg135","degAMO","degALL"};
	string dataTypes[3] {"data","acc","gen"};
	TFile* newFile;
	TTree* newTree;
	int sumEvents=0;
	int totalEvents=0;
	for (auto dataType : dataTypes){
		for (auto beamAngle : beamAngles){
			newFile = TFile::Open(("selectedPolarizations/"beamAngle+"_"+dataType+"_treeFlat_DSelector.root").c_str());
			newFile->GetObject(("selectedPolarizations/"beamAngle+"_"+dataType+"_tree_flat").c_str(),newTree);
			cout << "Entries in " << beamAngle << ": " << newTree->GetEntries() << endl;
			if (beamAngle != "degALL"){ sumEvents += newTree->GetEntries(); }  
			else { totalEvents = newTree->GetEntries(); } 
		}
	}
	cout << "Summed entries = " << sumEvents << " where totalEntries = " << totalEvents << endl;  
	assert(sumEvents==totalEvents); 
	cout << "Everything looks good!" << endl;
}
