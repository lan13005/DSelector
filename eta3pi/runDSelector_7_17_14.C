// macro to process analysis TTree with DSelector
// We cannot just run this macro, the library doesnt load properly. We can run the following two lines of code
//.x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C
//.x runDSelector.C

#include <iostream> 
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"



//R__LOAD_LIBRARY(/d/grid15/ln16/gluex_top/gluex_root_analysis/gluex_root_analysis-0.5/Linux_CentOS7-x86_64-gcc4.8.5/lib/libDSelector.so) 
R__LOAD_LIBRARY(/d/home/ln16/gluex_top/gluex_root_analysis/gluex_root_analysis_1.7.0/Linux_CentOS7-x86_64-gcc4.8.5/lib/libDSelector.so)
R__LOAD_LIBRARY(libDSelector.so)
   
void runDSelector_7_17_14(bool useproof = 1, string path = "") 
{
	cout << "Loaded using R__LOAD_LIBRARY" << endl;
	// Load DSelector library
	gROOT->ProcessLine(".x $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");
	// change the directory that proof saves the data to
	//gEnv->SetValue("ProofLite.Sandbox", "/d/grid15/ln16/.proof");
	int proof_Nthreads = 6;
	//int proof_Nthreads = 50;

	// open ROOT files and TTree
	TString nameOfTree = "pi0eta__eta_pi0pi0pi0__B4_F1_M7_M17_Tree"; // pi0eta__B4_Tree is the old one
	//TString nameOfTree = "pi0eta__eta_pi0pi0pi0__B4_M17_M7_Tree";
	TChain *chain = new TChain(nameOfTree);
	
	// **********************************************************************************	
	// ************************** ------ PI0ETA BELOW ---------**************************	
	// **********************************************************************************	
	
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a2_10M/eta3pi/tree_pi0eta__eta_pi0pi0pi0__B4_M17_M7.root");
	//
	// sean resolutions "normal"
	chain->Add("/d/home/sdobbs/grid13/MCtest/pi0eta_a2/hddm/tree_pi0eta__eta_pi0pi0pi0__B4_F1_M7_M17_00*");
	// sean resolutions "new"
	//chain->Add("/d/home/sdobbs/grid13/MCtest/71001_0/tree_pi0eta__eta_pi0pi0pi0__B4_F1_M7_M17.root");
	//chain->Add("/d/home/sdobbs/grid13/MCtest/71001_3.save/tree_pi0eta__eta_pi0pi0pi0__B4_F1_M7_M17.root");

	TString degAngle = "deg000";
	// should change the name below from data to reco when running over MC
	degAngle="pi0eta_seanResoution_reco_3pi0";
	//degAngle="pi0pi0_f2_reco";
	//degAngle = "pi0pi0_May2_";

	//  ===== Section is for pulling in data by polarization for asymmetry ===== /////
	// dividing allData into gorups by polarization. Actually the files we import have
	// the root files separated by polarizations. 
	//ifstream myReadFile;
	//TString degFile = "separateRuns/"+degAngle+".txt";
	//myReadFile.open(degFile);
	//char output[120];
	//if (myReadFile.is_open()) {
	//	// read the degAngle files until the end of it... and Add the file to the chain
	//	 while (!myReadFile.eof()) {
	//		myReadFile >> output;
	//		chain->Add(output);
	//		cout << output << endl;
	//	}
	//}
	/// ========================================================================///////


	//char oldHome[80]="HOME=/d/home/ln16";
	//char newHome[80]="HOME=/d/grid15/ln16";
	//putenv(newHome);
	//system("echo New Home At:");
	//system("echo $HOME");

	//string options = sample;
	string options = "";
	if(useproof) { // add TTree to chain and use PROOFLiteManager
		string outputHistFileName;
                string outputTreeFileName;
		outputHistFileName = degAngle+"_hists_DSelector.root";
		outputTreeFileName = degAngle+"_tree_DSelector.root"; 
	
		DPROOFLiteManager::Process_Chain(chain, "DSelector_eta3pi0.C+",  proof_Nthreads, outputHistFileName, outputTreeFileName, options);
	}
	else { // get TTree and use standard TTree::Process
		//TFile *f = TFile::Open(fileName);
		//TTree *tree = (TTree*)f->Get("omega_skim_Tree");
		chain->Process("DSelector_eta3pi0.C+", options.data());
		
	}
	
	//putenv(oldHome);
	//system("echo Returning to  Home At:");
	//system("echo $HOME");
	cout << "\033[1;31mMAKE SURE IF YOU ARE RUNNING ON PI0PI0 DATA WE USE mergePi0.root WITH THE CORRECT FILE!\n1)Change file name in mergePi0.root\n2)root -l -b -q mergePi0.root \033[0m\n";
	cout << "\033[1;31mRemember that Mpi0, Meta, Mpi0eta fundamental branches have a cut applied on them! This is for SPlotting to define a better region \033[0m\n";

	return;
}
