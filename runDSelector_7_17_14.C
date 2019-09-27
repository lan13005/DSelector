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
R__LOAD_LIBRARY(/d/grid13/gluex/gluex_top/gluex_root_analysis/gluex_root_analysis-1.3.0^rec1701v03/Linux_CentOS7-x86_64-gcc4.8.5/lib/libDSelector.so)
   
void runDSelector_7_17_14(bool useproof = 1, string path = "") 
{
	cout << "Loaded using R__LOAD_LIBRARY" << endl;
	// Load DSelector library
	gROOT->ProcessLine(".x $ROOT_ANALYSIS_HOME/scripts/Load_DSelector.C");
	// change the directory that proof saves the data to
	//gEnv->SetValue("ProofLite.Sandbox", "/d/grid15/ln16/.proof");
	int proof_Nthreads = 18;
	//int proof_Nthreads = 50;

	// open ROOT files and TTree
	//TString nameOfTree = "pi0eta__B3_F1_M7_M17_Tree"; // pi0eta__B4_Tree is the old one
	TString nameOfTree = "pi0eta__B4_M17_M7_Tree"; // pi0eta__B4_Tree is the old one
	//TString nameOfTree = "pi0pi0__B3_F1_U1_M7_Tree";
	TChain *chain = new TChain(nameOfTree);
	
	// **********************************************************************************	
	// ************************** ------ PI0ETA BELOW ---------**************************	
	// **********************************************************************************	
	
	// MC flat
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_noPlugin_Geant4_30730/tree_pi0eta__B4_M17_M7.root");
	//
	// a0 new hdgeant to test if mass shift still there
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/new_hdgeant/tree_pi0eta__B4_M17_M7.root");
	chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0_091519/tree_pi0eta__B4_M17_M7.root");
	
	// a0a2 recon_2017
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2_noPlugin_Geant4_30730/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2_noPlugin_Geant4_30730/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2pi1_noPlugin_Geant4_30730/tree_pi0eta__B4_M17_M7.root");
	//
	// vincent
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/vincent_noPlugin_Geant4_30461/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/vincent_noPlugin_Geant4_30730/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/vincent_noPlugin_Geant4_30461_withP/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/vincent_noPlugin_Geant4_30730_withP/tree_pi0eta__B4_M17_M7.root");



	//Real data for the 7_17_14 case MUCH LARGER
	//Trying to merge them now
	// allData
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/analysis-ver20/tree_pi0eta__B3_F1_M7_M17/merged/tree_pi0eta__B3_F1_M7_M17_03*.root");
	// This one contains the showerQuality variables
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/analysis-ver27/tree_pi0eta__B4_M17_M7/merged/tree_pi0eta__B4_M17_M7_03028*");
//	chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/analysis-ver27/tree_pi0eta__B4_M17_M7/merged/tree_pi0eta__B4_M17_M7*");
	//chain->Add("/d/grid15/ln16/pi0eta/121818/z_pi0eta_data/pi0eta_data_tree_DSelector.root");

	// **********************************************************************************	
	// ************************** ------ PI0PI0 BELOW ---------**************************	
	// **********************************************************************************	
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/analysis-ver20/tree_pi0pi0__B3_F1_M7/merged/tree_pi0pi0__B3_F1_M7_*");
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/tree_pi0pi0__B3_F1_U1_M7/merged/tree_pi0pi0__B3_F1_U1_M7_03*");
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/tree_pi0pi0__B3_F1_U1_M7/merged/tree_pi0pi0__B3_F1_U1_M7_03033*");
	//
	// MC with f2 resonance
	//chain->Add("/d/grid15/ln16/rootFiles/pi0pi0/Aug27Sep09/tree_pi0pi0__B3_F1_U1_M7.root");
	//chain->Add("/d/grid13/ln16/MC/pi0pi0_f2_interactive/hddm/combined/tree_pi0pi0__B3_F1_U1_M7.root");


	TString degAngle = "deg000";
	// should change the name below from data to reco when running over MC
	degAngle="pi0eta_a0_reco";
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
	
		DPROOFLiteManager::Process_Chain(chain, "DSelector_ver20.C+",  proof_Nthreads, outputHistFileName, outputTreeFileName, options);
	}
	else { // get TTree and use standard TTree::Process
		//TFile *f = TFile::Open(fileName);
		//TTree *tree = (TTree*)f->Get("omega_skim_Tree");
		chain->Process("DSelector_ver20.C+", options.data());
		
	}
	
	//putenv(oldHome);
	//system("echo Returning to  Home At:");
	//system("echo $HOME");
	cout << "\033[1;31mMAKE SURE IF YOU ARE RUNNING ON PI0PI0 DATA WE USE mergePi0.root WITH THE CORRECT FILE!\n1)Change file name in mergePi0.root\n2)root -l -b -q mergePi0.root \033[0m\n";
	cout << "\033[1;31mRemember that Mpi0, Meta, Mpi0eta fundamental branches have a cut applied on them! This is for SPlotting to define a better region \033[0m\n";

	return;
}
