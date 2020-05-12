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
	int proof_Nthreads = 24;
	//int proof_Nthreads = 50;

	// open ROOT files and TTree
	//TChain *chain = new TChain("ksks_Tree");
	//Check the name of the tree in the root files
	TString nameOfTree = "Thrown_Tree"; // pi0eta__B4_Tree is the old one
	TChain *chain = new TChain(nameOfTree);

	
	//a0
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0_091519/tree_thrown.root");

	// a0a2
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2pi1_noPlugin_Geant4_30730/tree_thrown.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2pi1_noPlugin_Geant4_30730/tree_thrown.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2_noPlugin_Geant4_30730/tree_thrown.root");

	// Flat
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_noPlugin_Geant4_30730_8to9GeV/tree_thrown.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_8GeVPlus_lustre_upTo3GeVResMass/tree_thrown.root");
	chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_2.1t/tree_thrown.root");

	// BA studies. Hddm filtered to make phi have some cos dependence
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_21t_hddmFiltered_8288_1628_uniquePolar/tree_thrown.root");

	// test
	//chain->Add("/d/grid13/ln16/MC/pi0eta_flat_2.3t/hddm/tree_thrown.root");
	string degAngle="degALL";
	string tag="gen";

	// The following two files work separately but not together...

	//string options = sample;
	string options = "";
	if(proof) { // add TTree to chain and use PROOFLiteManager
		//string outputHistFileName = Form("flatUpTo3GeVResMass_2_gen_hists_DSelector_pi0eta.root");//_GEANT4.root");
                //string outputTreeFileName = Form("flatUpTo3GeVResMass_2_gen_trees_DSelector_pi0eta.root");//_GEANT4.root");
		string outputHistFileName = degAngle+"_"+tag+"_2017_hists_DSelector.root";//_GEANT4.root");
                string outputTreeFileName = degAngle+"_"+tag+"_2017_trees_DSelector.root";//_GEANT4.root");
		DPROOFLiteManager::Process_Chain(chain, "DSelector_thrown.C+",  proof_Nthreads, outputHistFileName, outputTreeFileName, options);
	}
	else { // get TTree and use standard TTree::Process
		//TFile *f = TFile::Open(fileName);
		//TTree *tree = (TTree*)f->Get("omega_skim_Tree");
		chain->Process("DSelector_thrown.C+", options.data());
		
	}
	
	return;
}
