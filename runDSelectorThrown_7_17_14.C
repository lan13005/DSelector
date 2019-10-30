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

R__LOAD_LIBRARY(libDSelector) 
   
void runDSelectorThrown_7_17_14(bool proof = 1, string path = "") 
{
	// Load DSelector library
	gROOT->ProcessLine(".x $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");
	int proof_Nthreads = 12;
	//int proof_Nthreads = 50;

	// open ROOT files and TTree
	//TChain *chain = new TChain("ksks_Tree");
	//Check the name of the tree in the root files
	TString nameOfTree = "Thrown_Tree"; // pi0eta__B4_Tree is the old one
	TChain *chain = new TChain(nameOfTree);

//	ifstream inf;
//	inf.open("root_trees");
//	string line;
//	int nfile = 0;
//	while(inf >> line) {
//		//cout << line << endl;
//		//if(nfile++ > 20)
//		//	break;
//		//if(nfile < 10)
//		//	continue;
//
//		cout << line << endl;
//		chain->Add(line.c_str());
//	}
//	//For a single chain we do this, we can uncomment the above loop to loop over files in a text
//	Old Data when I first started
	//chain->Add("/d/home/sdobbs/grid13/gluex_data/RunPeriod-2017-01/analysis-ver08/tree_pi0eta__B4/merged/tree_pi0eta__B4_030298.root"); //large file in sean's directory
	
	//chain->Add("/d/home/ln16/pi0etaUpdated/tree_thrown_gen_amp_030730_7_17_14_GEANT4.root");//_GEANT4.root"); //the code I generated from ifarm
	//chain->Add("../rootFiles/pi0eta/rootFiles_Geant4_7M_1.8.0/tree_thrown_gen_amp_030730_**");
	//
	//a0
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0_091519/tree_thrown.root");

	// a0a2
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2pi1_noPlugin_Geant4_30730/tree_thrown.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2pi1_noPlugin_Geant4_30730/tree_thrown.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2_noPlugin_Geant4_30730/tree_thrown.root");

	// Flat
	chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_noPlugin_Geant4_30730/tree_thrown.root");

	// The following two files work separately but not together...
	//chain->Add("rootFiles_Geant4_7M/tree_thrown_gen_amp_030730_514.root");//_GEANT4.root"); //the code I generated from ifarm
	//chain->Add("rootFiles_Geant4_7M/tree_thrown_gen_amp_030730_513.root");//_GEANT4.root"); //the code I generated from ifarm

	//string options = sample;
	string options = "";
	if(proof) { // add TTree to chain and use PROOFLiteManager
		string outputHistFileName = Form("v20_flat_gen_hists_DSelector_pi0eta.root");//_GEANT4.root");
                string outputTreeFileName = Form("v20_flat_gen_trees_DSelector_pi0eta.root");//_GEANT4.root");
		DPROOFLiteManager::Process_Chain(chain, "DSelector_thrown_7_17_14.C+",  proof_Nthreads, outputHistFileName, outputTreeFileName, options);
	}
	else { // get TTree and use standard TTree::Process
		//TFile *f = TFile::Open(fileName);
		//TTree *tree = (TTree*)f->Get("omega_skim_Tree");
		chain->Process("DSelector_thrown_7_17_14.C+", options.data());
		
	}
	
	return;
}
