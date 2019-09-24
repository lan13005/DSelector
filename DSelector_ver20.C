#include "DSelector_ver20.h"
bool NoCut=0;
// degXXX where XXX = {000,045,090,135,All} where All is polarization independent. Actually anything other than the first 4 cases work but
// MUST BE ATLEAST 3 CHARACTERS LONG.
//string degAngle = "a0a2a2pi1_";
string degAngle="pi0eta_a0_reco";
bool showOutput = false;
bool showMassCalc = false;
bool onlyNamesPi0_1 = true; // true if we want to show only the histograms with _1 in their names so we can merge them with _2

// Here is an attempt to shorten the run time. We will only output the essential histograms
bool outputGJ_HEL=false;
bool outputThetaRegion=false;
bool outputTimingHists=false;
bool outputVanHoveCutEffect=false;
bool outputProdPlanePS_BeamAsym=false;
bool outputMandel_tBinnedInMass=false;
bool outputUERegion=false;
bool outputChiSqRegion=false;
bool outputPi0Resolutions=true;
bool outputMassBinnedE=false;
bool output_mEllipse_UE_ChiSq=true;
bool outputBkgSubOnBkgAndSignalRegions=false;
bool outputCorrelationBetweenMasses=false;
bool outputPhotonXandShowerLoc=false;
bool outputMassShift=true;

int itersToRun = 0;
int finalStateComboID=0;
void DSelector_ver20::Init(TTree *locTree)
{
	targetCenter = {0,0,65};

        if (is_pi0eta) {
                lowMass = 0.7;
                upMass = 2.0;
                etaProtonBaryonCut = 1.65;
                pi0ProtonBaryonCut = 2;
		binScale = (upMass-lowMass)/numBinsMass;
                //ellipseX = 0.134547; ellipseY = 0.541950; ellipseXr = 0.025449; ellipseYr = 0.069267;
                //using the kin data
                ellipseX = 0.135881; ellipseY = 0.548625; ellipseXr = 0.0160375; ellipseYr = 0.025671;
                ellipseXBS1 = 0.135881; ellipseYBS1 = 0.548625; ellipseXrBS1 = 0.022; ellipseYrBS1 = 0.06;
                ellipseXBS2 = 0.135881; ellipseYBS2 = 0.548625; ellipseXrBS2 = 0.045; ellipseYrBS2 = 0.165;
		ellipseXr_loose=0.0391; ellipseYr_loose=0.131;
		areaRatio = 0.067; //double checked these values;
        }
        else {
                lowMass = 0.3;
                upMass = 1.6;
                pi0ProtonBaryonCut = 1.8;
                etaProtonBaryonCut = 1.8;
                //ellipseX = 0.134547; ellipseY = 0.134547; ellipseXr = 0.025449; ellipseYr = 0.025449;
                //using the kin data
                ellipseX = 0.135881; ellipseY = 0.135881; ellipseXr = 0.0160375; ellipseYr = 0.0160375;
                ellipseXBS1 = 0.135881; ellipseYBS1 = 0.135881; ellipseXrBS1 = 0.022; ellipseYrBS1 = 0.022;
                ellipseXBS2 = 0.135881; ellipseYBS2 = 0.135881; ellipseXrBS2 = 0.045; ellipseYrBS2 = 0.045;
		ellipseXr_loose=0.0391; ellipseYr_loose=0.0391;
		areaRatio = 0.167; //double checked these values
        }


	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

        //USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
        dOutputFileName = degAngle+"DSelector_output.root"; //"" for none
        dOutputTreeFileName = degAngle+"tree_DSelector.root"; //"" for none
        dFlatTreeFileName = degAngle+"treeFlat_DSelector.root"; //output flat tree (one combo per tree entry), "" for none
        dFlatTreeName = degAngle+"tree_flat"; //if blank, default name will be chosen

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.4, 9.05));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	// dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	// dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/		

	// cutString will be used by alot of histograms when defining names 
	// Just defining alot of stuff before hand which will be used in the labeling of the histograms 
	std::string cutString;
	std::string cutsApplied="";
	std::string cutsBase="";
	std::vector<std::string> cutVariations;
	std::string massBin;

        dHist_BeamAngle = new TH1F("BeamAngle", "Beam Angle with no cuts applied;Beam Angle (GeV)", 180,0,180);
        dHist_BeamAngle->SetYTitle("Events / Degree");
	dHist_Cuts = new TH1F("CutsPassed", "Number of times a cut has been passed", 15,0,15);
	for (int i =0; i<3; ++i){
		if (is_pi0eta){
			dHist_checkEllipseBS[i] = new TH2F(("checkEllipseBS"+std::to_string(i)).c_str(), ";#pi_{0} Mass (GeV) with Events / 0.001 GeV;#eta Mass (GeV) with Events / 0.0025 GeV", atof(pi0BinRange[0].c_str()), atof(pi0BinRange[1].c_str()), atof(pi0BinRange[2].c_str()), atof(etaBinRange[0].c_str()), atof(etaBinRange[1].c_str()), atof(etaBinRange[2].c_str()));
		}
		else {
			dHist_checkEllipseBS[i] = new TH2F(("checkEllipseBS"+std::to_string(i)).c_str(), ";#pi_{0} Mass (GeV) with Events / 0.001 GeV;#pi_{0} Mass (GeV) with Events / 0.001 GeV", atof(pi0BinRange[0].c_str()), atof(pi0BinRange[1].c_str()), atof(pi0BinRange[2].c_str()), atof(pi0BinRange[0].c_str()), atof(pi0BinRange[1].c_str()), atof(pi0BinRange[2].c_str()));
		}
	}

// **************************************** START INITIALZING VARIABLES TO USE WITH HISTO BUILDING ********************************************//
// **************************************** START INITIALZING VARIABLES TO USE WITH HISTO BUILDING ********************************************//

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////// ********* DEFINING ALL HISTOGRAMS ****************/////////////////////////////////////////////

	// initialize id to -1. We have to have some initialization of histVals and histCuts or else it wont run properly.
	id = -1; // have to reinitialize id
	for (int i = 0; i < 2; ++i){
	       if (i<1) { cutString = ""; cutsToApply = 1; cutsBase="None";  applyAccSub=noAccSub;}
	       else { cutString = "Cut"; cutsToApply = allGeneralCutsPassed; cutsBase="GeneralCuts";  applyAccSub=weight;}
		cutsApplied=cutsBase;

               ++id; histVals[id] = {8, applyAccSub, locMissingMassSquared};
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mMMSq; cutsApplied="mMMsq";}
	       histList[id] = {("MissingMassSquared"+cutString).c_str(), ("Cuts="+cutsApplied+";Missing Mass Squared (GeV/c^{2})^{2}").c_str(), "200", "-0.2", "0.2", "Events / 0.002 GeV/c^{2}"};
		cutsApplied=cutsBase;
	
	       ++id; histVals[id] = {1, applyAccSub, locBeamE};
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mBeamE; cutsApplied="mBeamE";}
               histList[id] = {("BeamEnergy"+cutString).c_str(), ("Cuts="+cutsApplied+";Beam Energy (GeV)").c_str(), "120", "0", "12", "Events / 0.1 GeV" };
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {8, applyAccSub, locYDotZ_GJ};
	       histList[id] = {("yDotZ_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";Cos(theta) angle between Z and Y in GJ").c_str(), "100", "0", "1", "Events / 0.01"};
	       histCuts[id] = cutsToApply;
	       // Invariant Mass Hists
	       if (is_pi0eta) {
	       		++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
			if (cutString=="") { histCuts[id] = cutsToApply; }
			else { histCuts[id] = mMPi0P14; cutsApplied="mMPi0P14";}
	       		histList[id] = {("pi0proton1D"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
	       		++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin};
			cutsApplied=cutsBase;
	       		histList[id] = {("etaproton1D"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#etaproton) (GeV)").c_str(), "450", "0", "4.5", "Events / 0.01 GeV"};
	       		histCuts[id] = cutsToApply;
	       		++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locPi0Proton_Kin};
	       		histList[id] = {("pi0etaPi0Proton"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "90", "0.", "4.5", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#pi_{0}Proton) (GeV) with Events / 0.05 GeV"};
	       		histCuts[id] = cutsToApply;
	       		++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locEtaProton_Kin};
	       		histList[id] = {("pi0etaEtaProton"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "100", "0.", "5", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#etaProton) (GeV) with Events / 0.05 GeV"};
	       		histCuts[id] = cutsToApply;
		}
		else {
			// ********** NOT SURE IF WE SHOULD INCLUDE mMPi0P14 IN PI0PI0 REACION
	       		++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
			if (cutString=="") { histCuts[id] = cutsToApply; }
			else { histCuts[id] = mMPi0P14; cutsApplied="mMPi0P14";}
	       		histList[id] = {("pi0proton1D_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
	       		//histCuts[id] = cutsToApply;
			cutsApplied=cutsBase;
	       		++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin};
			if (cutString=="") { histCuts[id] = cutsToApply; }
			else { histCuts[id] = mMPi0P14; cutsApplied="mMPi0P14";}
	       		histList[id] = {("pi0proton1D_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
	       		//histCuts[id] = cutsToApply;
			cutsApplied=cutsBase;
	       		++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locPi0Proton_Kin};
	       		histList[id] = {("pi0pi0Pi0Proton_1"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "90", "0.", "4.5", "M(#pi_{0}#pi_{0}) (GeV) with Events / 0.025 GeV", "M(#pi_{0}Proton) (GeV) with Events / 0.05 GeV"};
	       		histCuts[id] = cutsToApply;
	       		++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locEtaProton_Kin};
	       		histList[id] = {("pi0pi0Pi0Proton_2"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "90", "0.", "4.5", "M(#pi_{0}#pi_{0}) (GeV) with Events / 0.025 GeV", "M(#pi_{0}Proton) (GeV) with Events / 0.05 GeV"};
	       		histCuts[id] = cutsToApply;
			
	       		++id; histVals[id] = {5, applyAccSub, locPi0E_Kin};
	       		histList[id] = {("pi0E_1"+cutString).c_str(), ("Cuts="+cutsApplied+"*pSelectf2;E (GeV)").c_str(), "100", "0", "10", "Events / 0.1 GeV"};
	       		histCuts[id] = cutsToApply;
	       		++id; histVals[id] = {6, applyAccSub, locEtaE_Kin};
	       		histList[id] = {("pi0E_2"+cutString).c_str(), ("Cuts="+cutsApplied+"*pSelectf2;E (GeV)").c_str(), "100", "0", "10", "Events / 0.1 GeV"};
	       		histCuts[id] = cutsToApply;
		}

	       // Kinematic Hists
	       ++id; histVals[id] = {8, applyAccSub, locPhi}; 
	       histList[id] = {("prodPlanePSphi"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi").c_str(), "180", "-540", "540", "Events / 6 degrees"};
	       histCuts[id] = cutsToApply;
	       //++id; histVals[id] = {8, applyAccSub, cosTheta_decayPlane_hel};
	       //histList[id] = {("decayPlane_cosTheta_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of decay plane").c_str(), "100","-1","1", "Events / 0.02"};
	       //histCuts[id] = cutsToApply;
	       //++id; histVals[id] = {8, applyAccSub, phi_decayPlane_hel};
	       //histList[id] = {("decayPlane_phi_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of decay plane").c_str(), "160","-1.6","1.6", "Events / 0.02"};
	       //histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {8, applyAccSub, omega};
	       histList[id] = {("vanHove_omega"+cutString).c_str(), ("Cuts="+cutsApplied+";Van Hove Omega Plot").c_str(), "120","-360","360", "60", "Events / 6 degrees"};
	       histCuts[id] = cutsToApply;


	       if(outputGJ_HEL){
		       if (is_pi0eta) {
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_pi0_hel};
		       		histList[id] = {("pi0_cosTheta_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_eta_hel};
		       		histList[id] = {("eta_cosTheta_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #eta").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_pi0eta_hel};
		       		histList[id] = {("pi0eta_cosTheta_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#eta").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_eta_hel};
		       		histList[id] = {("eta_phi_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #eta").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_pi0_hel};
		       		histList[id] = {("pi0_phi_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_pi0eta_hel};
		       		histList[id] = {("pi0eta_phi_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}+#eta").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, theta_pi0_GJ};
		       		histList[id] = {("pi0_theta_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, theta_eta_GJ};
		       		histList[id] = {("eta_theta_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta of #eta").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_eta_GJ};
		       		histList[id] = {("eta_cosTheta_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #eta").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_pi0_GJ};
		       		histList[id] = {("pi0_cosTheta_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_eta_GJ};
		       		histList[id] = {("eta_cosTheta_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta)Events / 0.01 GeV", "Cos(#theta) of #eta Events / 0.2"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_pi0_GJ};
		       		histList[id] = {("pi0_cosTheta_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #pi_{0} vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta) Events / 0.01 GeV", "Cos(#theta) of #pi_0 Events / 0.2"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, phi_pi0_GJ};
		       		histList[id] = {("pi0_phi_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+" phi of #pi_{0} vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "90","-180","180", "M(#pi_{0}#eta) Events / 0.01 GeV", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, phi_eta_GJ};
		       		histList[id] = {("eta_phi_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+" phi of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "90","-180","180", "M(#pi_{0}#eta)Events / 0.01 GeV", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_eta_GJ};
		       		histList[id] = {("eta_phi_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #eta").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_pi0_GJ};
		       		histList[id] = {("pi0_phi_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
			}
			else {
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_pi0_hel};
		       		histList[id] = {("pi0_cosTheta_hel_1"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_eta_hel};
		       		histList[id] = {("pi0_cosTheta_hel_2"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_pi0eta_hel};
		       		histList[id] = {("pi0pi0_cosTheta_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_eta_hel};
		       		histList[id] = {("pi0_phi_hel_1"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_pi0_hel};
		       		histList[id] = {("pi0_phi_hel_2"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_pi0eta_hel};
		       		histList[id] = {("pi0pi0_phi_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}+#pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, theta_pi0_GJ};
		       		histList[id] = {("pi0_theta_GJ_1"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, theta_eta_GJ};
		       		histList[id] = {("pi0_theta_GJ_2"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_eta_GJ};
		       		histList[id] = {("pi0_cosTheta_GJ_1"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, cosTheta_pi0_GJ};
		       		histList[id] = {("pi0_cosTheta_GJ_2"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_eta_GJ};
		       		histList[id] = {("pi0_cosTheta_GJvsM_1"+cutString).c_str(), ("Cuts="+cutsApplied+" cosTheta vs M(#pi_{0}#pi_{0})").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta) Events / 0.01 GeV", "Cos(#theta) Events / 0.2"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_pi0_GJ};
		       		histList[id] = {("pi0_cosTheta_GJvsM_2"+cutString).c_str(), ("Cuts="+cutsApplied+" cosTheta vs M(#pi_{0}#pi_{0})").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta) Events / 0.01 GeV", "Cos(#theta) Events / 0.2"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, phi_pi0_GJ};
		       		histList[id] = {("pi0_phi_GJvsM_1"+cutString).c_str(), ("Cuts="+cutsApplied+" phi of #pi_{0} vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "90","-180","180", "M(#pi_{0}#eta) Events / 0.01 GeV", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, phi_eta_GJ};
		       		histList[id] = {("pi0_phi_GJvsM_2"+cutString).c_str(), ("Cuts="+cutsApplied+" phi of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "90","-180","180", "M(#pi_{0}#eta)Events / 0.01 GeV", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_eta_GJ};
		       		histList[id] = {("pi0_phi_GJ_1"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, phi_pi0_GJ};
		       		histList[id] = {("pi0_phi_GJ_2"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
		       		histCuts[id] = cutsToApply;
			}
	       		if (is_pi0eta) { 
	       			++id; histVals[id] = {8, applyAccSub, cosTheta_pi0eta_CM};
	       			histList[id] = {("pi0eta_cosTheta_CM"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#eta").c_str(), "100","-1.00","1.00", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {3, applyAccSub, cosTheta_pi0_CM};
	       			histList[id] = {("pi0_cosTheta_CM"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {7, applyAccSub, cosTheta_eta_CM};
	       			histList[id] = {("eta_cosTheta_CM"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #eta").c_str(), "100","-1.00","1.00", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, phi_pi0eta_CM};
	       			histList[id] = {("pi0eta_phi_CM"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}+#eta").c_str(), "90","-180","180", "Events / 4 degrees"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {3, applyAccSub, phi_pi0_CM};
	       			histList[id] = {("pi0_phi_CM"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {7, applyAccSub, phi_eta_CM};
	       			histList[id] = {("eta_phi_CM"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #eta").c_str(), "90","-180","180", "Events / 4 degrees"};
	       			histCuts[id] = cutsToApply;
			}
			else {
	       			++id; histVals[id] = {8, applyAccSub, cosTheta_pi0eta_CM};
	       			histList[id] = {("pi0pi0_cosTheta_CM"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {3, applyAccSub, cosTheta_pi0_CM};
	       			histList[id] = {("pi0_cosTheta_CM_1"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {7, applyAccSub, cosTheta_eta_CM};
	       			histList[id] = {("pi0_cosTheta_CM_2"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, phi_pi0eta_CM};
	       			histList[id] = {("pi0pi0_phi_CM"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}+#pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {3, applyAccSub, phi_pi0_CM};
	       			histList[id] = {("pi0_phi_CM_1"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {7, applyAccSub, phi_eta_CM};
	       			histList[id] = {("pi0_phi_CM_2"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}").c_str(), "90","-180","180", "Events / 4 degrees"};
	       			histCuts[id] = cutsToApply;
			}
		}
	
	       // Invariant Mass Hists EBeam > 8 GeV
	//       cutsApplied=cutsBase+"*pBeamE8GeVPlus";
	//       ++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
	//       histList[id] = {("pi0proton1D8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
	//       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	//       ++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin}; 
	//       histList[id] = {("etaproton1D8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#etaproton) (GeV)").c_str(), "450", "0", "4.5", "Events / 0.01 GeV"};
	//       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	//       ++id; histVals[id] = {11, applyAccSub, locEtaProton_Kin,locPi0Eta_Kin}; 
	//       histList[id] = {("pi0etaProton8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "90", "0.", "4.5", "160", "0.", "4", "M(#etaProton) (GeV) with Events / 0.05 GeV", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV"};
	//       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	//       ++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locPi0Proton_Kin};
	//       histList[id] = {("pi0etaPi0Proton8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "90", "0.", "4.5", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#pi_{0}Proton) (GeV) with Events / 0.05 GeV"};
	//       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	//       ++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locEtaProton_Kin};
	//       histList[id] = {("pi0etaEtaProton8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "100", "0.", "5", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#etaProton) (GeV) with Events / 0.05 GeV"};
	//       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	//	cutsApplied=cutsBase;
	
	       // Charged Track Hists
	       ++id; histVals[id] = {2, applyAccSub, locRProton}; 
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mRProton; cutsApplied="mRProton";}
	       histList[id] = {("RadiusProton"+cutString).c_str(), ("Cuts="+cutsApplied+";Radius(proton) (cm)").c_str(), "200", "0" , "10", "Events / 0.05 cm"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {2, applyAccSub, locdzProton};
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mZMin; cutsApplied="mZMin";}
	       histList[id] = {("dzProton"+cutString).c_str(), ("Cuts="+cutsApplied+";z(Proton) (cm)").c_str(), "160", "0" , "160", "Events / 1 cm"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {2, applyAccSub, locdEdxCDCProton};
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mdEdxCDC; cutsApplied="mdEdxCDC";}
	       histList[id] = {("dEdxProtonCDC"+cutString).c_str(), ("Cuts="+cutsApplied+";dEdx(proton) GeV/cm").c_str(), "200", "0." , "0.00003", "Events / 1.5E-7 GeV/cm"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {2, applyAccSub, locPzProton}; 
	       histList[id] = {("PzProton"+cutString).c_str(), ("Cuts="+cutsApplied+";Pz GeV").c_str(), "200", "0" , "5", "Events / 0.025 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locXProton,locYProton}; 
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mRProton; cutsApplied="mRProton";}
	       histList[id] = {("XYplaneProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "100", "-4", "4", "100", "-4", "4", "x(proton) (cm) with Events / 0.08 cm", "y(proton) (cm) with Events / 0.08 cm"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {2, applyAccSub, locRProton,locdzProton};
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mRProtonZMin; cutsApplied="mRProtonZMin";}
	       histList[id] = {("RZplaneProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "100", "0", "4", "100", "0", "100", "R(proton) (cm) wtih Events / 0.04 cm","z(proton) (cm) with Events / 1 cm"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {2, applyAccSub, locMagP3Proton,locdEdxCDCProton};
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mdEdxCDC; cutsApplied="dzRP";} // recall this is no where the final value of histCuts is set!
	       histList[id] = {("P3dEdxCDCProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "100", "0", "4", "100", "0", "0.00003", "Momentum(proton) CDC (GeV/c) with Events / 0.04 GeV/c", "dEdx(proton) CDC (GeV/cm) with Events / 3E-7"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {2, applyAccSub, locMagP3Proton,locdEdxFDCProton};
	       histList[id] = {("P3dEdxFDCProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "100", "0", "4", "100", "0", "0.00003", "Momentum(proton) FDC (GeV/c) with Events / 0.04 GeV/c", "dEdx(proton) FDC (GeV/cm) with Events / 3E-7"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locPzProton,locPtProton}; 
	       histList[id] = {("PzPtProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "125", "0", "5", "125", "0", "2.5", "Pz(proton) (GeV/c) with Events / 0.04 GeV/c", "Pt(proton) (GeV/c) with Events / 0.02 GeV/c"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locPzProton,locPolarAngleProton};
	       histList[id] = {("PzThetaProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "125", "0", "5", "250", "0", "100", "Pz(proton) (GeV/c) with Events / 0.04 GeV/c", "Theta(proton) (degrees) with Events / 0.4 degrees"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locPtProton,locPolarAngleProton};
	       histList[id] = {("PtThetaProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "125", "0", "5", "250", "0", "100", "Pz(proton) (GeV/c) with Events / 0.04 GeV/c", "Theta(proton) (degrees) with Events / 0.4 degrees"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	       histList[id] = {("PolarAngleProton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(Proton) degrees").c_str()," 200", "0", "100", "Events / 0.5 degrees"};
	       histCuts[id] = cutsToApply;
		if (outputThetaRegion){
	       		++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	       		histList[id] = {("ThetaRegion1Proton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(degrees)").c_str(), "125", "0", "100", "Events / 0.08 degrees"};
	       		histCuts[id] = cutsToApply*pReg1; 
	       		++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	       		histList[id] = {("ThetaRegion2Proton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(degrees)").c_str(), "125", "0", "100", "Events / 0.08 degrees"};
	       		histCuts[id] = cutsToApply*pReg2; 
	       		++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	       		histList[id] = {("ThetaRegion3Proton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(degrees)").c_str(), "125", "0", "100", "Events / 0.08 degrees"};
	       		histCuts[id] = cutsToApply*pReg3; 
	       		++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	       		histList[id] = {("ThetaRegion4Proton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(degrees)").c_str(), "125", "0", "100", "Events / 0.08 degrees"};
	       		histCuts[id] = cutsToApply*pReg4; 
		}
	
	
	       // Neutral Track Hists
	       ++id; histVals[id] = {13, applyAccSub, photonEnergies[0]};
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mPhotonE; cutsApplied="mPhotonE";}
	       histList[id] = {("PhotonShowerE1"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy(#gamma) GeV").c_str(), "200", "0" , "10", "Events / 0.05 GeV"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {13, applyAccSub, photonThetas[0]}; 
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mPhotonTheta; cutsApplied="mPhotonTheta";}
	       histList[id] = {("PhotonShowerTheta1"+cutString).c_str(), ("Cuts="+cutsApplied+";PolarAngle(#gamma) degrees").c_str(), "300", "0" , "150", "Events / 0.5 degrees"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {13, applyAccSub, photonEnergies[0], photonThetas[0]};
	       histList[id] = {("thetaEPhotonShower"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "10", "400", "0", "150", "E(#gamma) (GeV) wtih Events / 0.05", "#theta(#gamma) (radians) with Events / 0.25"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, photonXs_Shower[0], photonYs_Shower[0]}; 
	       histList[id] = {("XYPhotonShower"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "260", "-130", "130", "260", "-130", "130", "x(cm) with Events / 1 cm", "y(cm) with Events / 1 cm"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, photonXs_Shower[0], photonYs_Shower[0]}; 
	       histList[id] = {("XYPhotonShowerBCAL"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "260", "-130", "130", "260", "-130", "130", "x(cm) with Events / 1 cm", "y(cm) with Events / 1 cm"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, photonPhis[0], photonThetas[0]};
	       histList[id] = {("ThetaPhiPhotonShower"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "360", "-180", "180", "130", "0", "130", "#phi(degrees) with Events / 1 cm", "#theta(degrees) with Events / 1 cm"};
	       histCuts[id] = cutsToApply;
	       // These below are tracked differently than the charged tracks and neutral tracks so they have their own unique identifiers.
	       ++id; histVals[id] = {12, applyAccSub, locPhotonDijFCAL}; 
	       if (cutString=="") { histCuts[id] = cutsToApply; id_noCutDij3=id;}
	       else { histCuts[id] = mdij3; cutsApplied="mdij3";} 
	       histList[id] = {("3DistanceBetweenPhotonsFCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";FCAL - Distance between Pairs of Photons (cm)").c_str(), "500", "0", "250", "Events / 0.5 cm"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {14, applyAccSub, locPhotonDijBCAL}; 
	       histList[id] = {("distanceOnCylinderBCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";BCAL - Distance on cylinder between Pairs of Photons (cm)").c_str(), "400", "0", "400", "Events / 1 cm"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {14, applyAccSub, locPhotonAij}; 
	       histList[id] = {("angleBetweenPhotonsBCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";Angle between photons in BCAL (degrees)").c_str(), "360", "0", "180", "Entries / 0.5cm"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {14, applyAccSub, locPhotonZij, locPhotonPij}; 
	       histList[id] = {("deltaZvsPhiPhotonBCAL"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "400", "100", "0", "200", "#DeltaZ  (cm) with Events / 2 cm",  "#Delta#phi between showers in BCAL (degrees) with Events / 2 degrees"};
	       histCuts[id] = cutsToApply;
	
	       // Timing
	       if(outputTimingHists){
	       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton,RFtimeProton}; 
	       		histList[id] = {("timeProtonRFvsTheta"+cutString).c_str(),("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees","T_{Proton} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply;
	       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	       		histList[id] = {("timeProtonRFvsP3"+cutString).c_str(),("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c", "T_{Proton} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply;
	       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton,RFtimeProton}; 
	       		histList[id] = {("timeBCALRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees","T_{BCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	       		histList[id] = {("timeBCALRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120","-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{BCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
	       		histList[id] = {("timeTOFRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees", "T_{TOF} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	       		histList[id] = {("timeTOFRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{TOF} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
	       		histList[id] = {("timeFCALRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200"," 0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 cm", "T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	       		histList[id] = {("timeFCALRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c", "T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
	       		histList[id] = {("timeSTARTRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees", "T_{START} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	       		histList[id] = {("timeSTARTRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{START} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
	       		histList[id] = {("timeSYSNULLRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees","T_{SYSNULL} - T_{RF} (ns) with Events / 0.1 ns "};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	       		histList[id] = {("timeSYSNULLRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6" , "6","Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{SYSNULL} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 



	       		++id; histVals[id] = {8, applyAccSub, photonThetas[0], photonDeltaTs[0]}; 
	       		histList[id] = {("timePhotonBCALFCALRFvsTheta"+cutString).c_str(),("Cuts="+cutsApplied).c_str(),"300","0","120","120","-6","6", "#theta(Proton) degrees with Events / 0.4 degrees","T_{BCAL/FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	       		++id; histVals[id] = {8, applyAccSub, photonDeltaTs[0]}; 
	       		histList[id] = {("timePhotonBCALFCALRF"+cutString).c_str(),("Cuts="+cutsApplied+";T_{BCAL/FCAL} - T_{RF} (ns)").c_str(),"120","-6","6", "Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		histCuts[id] = cutsToApply;// recall this is not where we actually fill the histCuts, its just for checking sizes of things at this point... 
	       		++id; histVals[id] = {8, applyAccSub, photonThetas[0], photonDeltaTs[0]}; 
	       		histList[id] = {("timePhotonFCALRFvsTheta"+cutString).c_str(),("Cuts="+cutsApplied).c_str(),"60","0","12","120","-6","6", "#theta(Proton) degrees with Events / 0.2 degrees","T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, photonEnergies[0], photonDeltaTs[0]};
	       		histList[id] = {("timePhotonFCALRFvsEnergy"+cutString).c_str(),("Cuts="+cutsApplied).c_str(),"200","0","12","120","-6","6", "Photon momentum (GeV/c) with Events / 0.06 GeV/c", "T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply*pPhotonInFCAL[0]; 
	       		++id; histVals[id] = {8, applyAccSub, photonDeltaTs[0]};                                                                                    
	       		histList[id] = {("timePhotonFCALRF"+cutString).c_str(),("Cuts="+cutsApplied+";T_{FCAL} - T_{RF} (ns)").c_str(),"120","-6","6", "Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, photonEnergies[0], photonDeltaTs[0]};
	       		histList[id] = {("timePhotonBCALRFvsEnergy"+cutString).c_str(),("Cuts="+cutsApplied).c_str(),"200","0","6","120","-6","6", "Photon momentum (GeV/c) with Events / 0.03 GeV/c","T_{BCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply; 
	       		++id; histVals[id] = {8, applyAccSub, photonDeltaTs[0]}; 
	       		histList[id] = {("timePhotonBCALRF"+cutString).c_str(),("Cuts="+cutsApplied+";T_{BCAL} - T_{RF} (ns)").c_str(),"120","-6","6", "Events / 0.1 ns"};
	       		histCuts[id] = cutsToApply;                                                                                                  
		}
	

	       // Non track related hists
	       ++id; histVals[id] = {8, applyAccSub, locCLKinFit}; 
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mChiSq; cutsApplied="mChiSq";} 
	       histList[id] = {("P4CLKinFit"+cutString).c_str(), ("Cuts="+cutsApplied+";Confidence Level").c_str() , "200", "0", ".30", "Entries / 0.0015"};
	       cutsApplied=cutsBase;
	       ++id; histVals[id] = {8, applyAccSub, locUnusedEnergy};
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mUE; cutsApplied="mUE";} 
	       histList[id] = {("UnusedEnergy"+cutString).c_str(), ("Cuts="+cutsApplied+";Unused Energy (GeV)").c_str(), "100", "0", "1", "Entries / 0.01 GeV"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {8, applyAccSub, locChiSqKinFit}; 
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mChiSq; cutsApplied="mChiSq";} 
	       histList[id] = {("P4ChiSqKinFit"+cutString).c_str(), ("Cuts="+cutsApplied+";Chi Squared").c_str(), "150", "0", "300", "Entries / 1"};
		cutsApplied=cutsBase;
	       ++id; histVals[id] = {8, applyAccSub, locDOFKinFit}; 
	       if (cutString=="") { histCuts[id] = cutsToApply; }
	       else { histCuts[id] = mChiSq; cutsApplied="mChiSq";} 
	       histList[id] = {("P4DOFKinFit"+cutString).c_str(), ("Cuts="+cutsApplied+";DOF").c_str(), "30", "0", "30", "Entries / 1"};
		cutsApplied=cutsBase;

		//Shower shape variables
		// neutral showers
	       ++id; histVals[id] = {13, applyAccSub, E1E9_FCAL[0]}; 
	       histList[id] = {("E1E9_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";E1E9_FCAL").c_str(), "100", "-1", "2", "Entries / 0.03"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, E9E25_FCAL[0]}; 
	       histList[id] = {("E9E25_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";E9E25_FCAL").c_str(), "100", "-1", "2", "Entries / 0.03"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, SumU_FCAL[0]}; 
	       histList[id] = {("SumU_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SumU_FCAL").c_str(), "200", "-1", "9", "Entries / 0.5"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, SumV_FCAL[0]}; 
	       histList[id] = {("SumV_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SumV_FCAL").c_str(), "200", "-1", "9", "Entries / 0.5"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, DOCA_FCAL[0]}; 
	       histList[id] = {("DOCA_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";DOCA_FCAL").c_str(), "160", "0", "16", "Entries / 0.1 cm"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, showerQuality_FCAL[0]}; 
	       histList[id] = {("showerQuality_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";showerQuality_FCAL").c_str(), "50", "0", "1", "Entries / 0.02"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, Energy_BCALPreshower[0]}; 
	       histList[id] = {("Energy_BCALPreshower"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy_BCALPreshower").c_str(), "200", "-1", "9", "Entries / 0.5 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, Energy_BCAL[0]}; 
	       histList[id] = {("Energy_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, SigLong_BCAL[0]}; 
	       histList[id] = {("sigLong_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SigLong_BCAL").c_str(), "200", "0", "10", "Entries / 0.05 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, SigTrans_BCAL[0]}; 
	       histList[id] = {("sigTrans_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTrans_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, SigTheta_BCAL[0]}; 
	       histList[id] = {("sigTheta_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTheta_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, DeltaPhi_BCAL[0]}; 
	       histList[id] = {("DeltaPhi_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";DeltaPhi_BCAL").c_str(), "100", "-5", "5", "Entries / 0.1 degree"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {13, applyAccSub, DeltaZ_BCAL[0]}; 
	       histList[id] = {("DeltaZ_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";DeltaZ_BCAL").c_str(), "200", "-100", "100", "Entries / 1 cm"};
	       histCuts[id] = cutsToApply;
		
		// charged showers
	       ++id; histVals[id] = {2, applyAccSub, locE1E9_FCAL_proton}; 
	       histList[id] = {("E1E9_FCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";E1E9_FCAL").c_str(), "100", "-1", "2", "Entries / 0.03"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locE9E25_FCAL_proton}; 
	       histList[id] = {("E9E25_FCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";E9E25_FCAL").c_str(), "100", "-1", "2", "Entries / 0.03"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locSumU_FCAL_proton}; 
	       histList[id] = {("SumU_FCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";SumU_FCAL").c_str(), "200", "-1", "9", "Entries / 0.05"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locSumV_FCAL_proton}; 
	       histList[id] = {("SumV_FCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";SumV_FCAL").c_str(), "200", "-1", "9", "Entries / 0.05"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locEnergy_BCALPreshower_proton}; 
	       histList[id] = {("Energy_BCALPreshowerproton"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy_BCALPreshower").c_str(), "100", "-1", "9", "Entries / 0.05 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locEnergy_BCAL_proton}; 
	       histList[id] = {("Energy_BCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locSigLong_BCAL_proton}; 
	       histList[id] = {("sigLong_BCAL_proton"+cutString).c_str(), ("Cuts="+cutsApplied+";SigLong_BCAL").c_str(), "200", "0", "10", "Entries / 0.05 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locSigTrans_BCAL_proton}; 
	       histList[id] = {("sigTrans_BCAL_proton"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTrans_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	       histCuts[id] = cutsToApply;
	       ++id; histVals[id] = {2, applyAccSub, locSigTheta_BCAL_proton}; 
	       histList[id] = {("sigTheta_BCAL_proton"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTheta_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	       histCuts[id] = cutsToApply;
	}


	if (outputVanHoveCutEffect){
		cutVariations = {"", "Cut", "pVanHove"};
		for (int j=0; j<3; ++j){
		       if (j==0) { cutsToApply = 1; cutsApplied="None"; applyAccSub=noAccSub;}
		       if (j==1) { cutsToApply = allGeneralCutsPassed; cutsApplied="GeneralCuts"; applyAccSub=weight;}
		       if (j==2) { cutsToApply = allGeneralCutsPassed*pVanHove; cutsApplied="GeneralCuts+pVanHove"; applyAccSub=weight;}
			cutString = cutVariations[j];
			if (is_pi0eta) {
				++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, mandelstam_tp};
				histList[id] = {("pi0eta_tVspi0etaMass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "100", "0.8", "1.8", "100" ,"0", "6", "M(#pi_{0}#eta) with Events / 0.01  GeV", "-t' momentum transfer of #pi_{0}+#eta with Events / 0.01 GeV"};
				histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {11, applyAccSub, mandelstam_tp};
		       		histList[id] = {("pi0eta_t"+cutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using recoil+target").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
		       		histCuts[id] = cutsToApply;
			}
			else {
				++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, mandelstam_tp};
				histList[id] = {("pi0pi0_tVspi0pi0Mass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "100", "0.3", "2", "100" ,"0", "6", "M(#pi_{0}#pi_{0}) with Events / 0.017  GeV", "-t' momentum transfer of #pi_{0}+#pi_{0} with Events / 0.06 GeV"};
				histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {11, applyAccSub, mandelstam_tp};
		       		histList[id] = {("pi0pi0_t"+cutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using recoil+target").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
		       		histCuts[id] = cutsToApply;
			}
		       ++id; histVals[id] = {8, applyAccSub, vanHove_x, vanHove_y};
		       histList[id] = {("vanHove"+cutString).c_str(), ("Cuts="+cutsApplied+";Van Hove Plot").c_str(), "60","-3","3", "60", "-3", "3", "Events / 0.1 degrees", "Events / 0.1 degrees"};
		       histCuts[id] = cutsToApply;
		}
	}


	cutVariations = {"", "Cut", "mUE", "mUEChiSq100"};
	for (int j=0; j<4; ++j){
	       if (j==0) { cutsToApply = 1; cutsApplied="None"; applyAccSub=noAccSub;}
	       if (j==1) { cutsToApply = allGeneralCutsPassed; cutsApplied="GeneralCuts"; applyAccSub=weight;}
	       if (j==2) { cutsToApply = mUE; cutsApplied="mUE"; applyAccSub=weight;}
	       if (j==3) { cutsToApply = mUEChiSq*chiSq100; cutsApplied="mUEChiSq+ChiSq100"; applyAccSub=weight;}
	       cutString = cutVariations[j];
	       ++id; histVals[id] = {8, applyAccSub, locNumExtraNeutralShowers}; 
	       histList[id] = {("numExtraNeutralShowers"+cutString).c_str(),("Cuts="+cutsApplied+";Number of extra showers").c_str(),"10","0","10", "Events / 1"};
	       histCuts[id] = cutsToApply;                                                                                                  
	}


	if(outputProdPlanePS_BeamAsym){
		for (int i = 0; i < numBinsMass; ++i) {
			cutString="Cut";
		       applyAccSub = weight;
		       // these will always have the general cuts applied, we already have the noCuts version of this (not binned though) 
		       cutsToApply = mBeamE*pBeamAsymE;
			cutsApplied="mBeamE*pBeamAsymE*p_phiMassBinned[i]";
		       massBin = "Bin"+std::to_string(i);
		       ++id; histVals[id] = {8, applyAccSub, locPhi}; 
		       histList[id] = {("prodPlanePSphi"+massBin).c_str(), ("Cuts="+cutsApplied+"+phi binned M(#pi_{0}#eta), General Cuts + 8 GeV < EBeam < 8.7 GeV").c_str(), "180", "-540","540", "Events / 6 degrees"};
		       histCuts[id] = cutsToApply*p_phiMassBinned[i];
			if (is_pi0eta){
		       ++id; histVals[id] = {8, applyAccSub, locPhi, locPi0Eta_Kin}; 
		       histList[id] = {("prodPlanePSphiVsMass"+massBin).c_str(), ("Cuts="+cutsApplied+"+phi vs mass binned in M(#pi_{0}#eta), General Cuts + 8 GeV < EBeam < 8.7 GeV").c_str(), "180", "-540","540", "300", "0", "3", "Events / 6 degrees", "Events / 0.1 GeV"};
		       histCuts[id] = cutsToApply*p_phiMassBinned[i];
			}
			else {
		       ++id; histVals[id] = {8, applyAccSub, locPhi, locPi0Eta_Kin}; 
		       histList[id] = {("prodPlanePSphiVsMass"+massBin).c_str(), ("Cuts="+cutsApplied+"+phi vs mass binned in M(#pi_{0}#pi_{0}), General Cuts + 8 GeV < EBeam < 8.7 GeV").c_str(), "180", "-540","540", "300", "0", "3", "Events / 6 degrees", "Events / 0.1 GeV"};
		       histCuts[id] = cutsToApply*p_phiMassBinned[i];
			}

			for (int iDelta=0; iDelta<2; ++iDelta){
		       		cutsToApply = allGeneralCutsPassed;
				if (is_pi0eta){
		       			++id; histVals[id] = {8, applyAccSub, 1};
		       			histList[id] = {("numPassConeDelta"+std::to_string(delta[iDelta])+"Mass"+massBin+cutString).c_str(), ("Cuts="+cutsApplied+";pi0InCone,etaInCone,largeAngle,withinCone,notWithinCone").c_str(),"5","0","5","Events / 1"};
		       			histCuts[id] = cutsToApply;
				}
				else {
		       			++id; histVals[id] = {8, applyAccSub, 1};
		       			histList[id] = {("numPassConeDelta"+std::to_string(delta[iDelta])+"Mass"+massBin+cutString).c_str(), ("Cuts="+cutsApplied+";pi0InCone_1,piInCone_2,largeAngle,withinCone,notWithinCone").c_str(),"5","0","5","Events / 1"};
		       			histCuts[id] = cutsToApply;
				}
			}
		}
	}
	
	if(outputMandel_tBinnedInMass){
		for (int i = 0; i < numBinsMass_t; ++i){
		       // these will always have the general cuts applied, we already have the noCuts version of this (not binned though) 
		       applyAccSub = weight;
		       cutsToApply = allGeneralCutsPassed;
			cutsApplied="GeneralCuts*p_tMassBinned[i]";
		       massBin = "Bin"+std::to_string(i);
			if (is_pi0eta){
		       		++id; histVals[id] = {11, applyAccSub, mandelstam_tp};  
		       		histList[id] = {("pi0eta_t"+massBin).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer binned in M(#pi_{0}#eta)").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
		       		histCuts[id] = cutsToApply*p_tMassBinned[i];
		       		++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, mandelstam_tp};
		       		histList[id] = {("pi0eta_tVspi0etaMass"+massBin).c_str(), ("Cuts="+cutsApplied).c_str(), "100", "0.8", "1.8", "100" ,"0", "2", "M(#pi_{0}#eta) with Events / 0.01  GeV","-t' momentum transfer of #pi_{0}+#eta with Events / 0.02 GeV"};
		       		histCuts[id] = cutsToApply*p_tMassBinned[i];
			}
			else {
		       		++id; histVals[id] = {11, applyAccSub, mandelstam_tp};  
		       		histList[id] = {("pi0pi0_t"+massBin).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer binned in M(#pi_{0}#pi_{0})").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
		       		histCuts[id] = cutsToApply*p_tMassBinned[i];
		       		++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, mandelstam_tp};
		       		histList[id] = {("pi0pi0_tVspi0pi0Mass"+massBin).c_str(), ("Cuts="+cutsApplied).c_str(), "100", "0.3", "1.3", "100" ,"0", "2", "M(#pi_{0}#pi_{0}) with Events / 0.01  GeV","-t' momentum transfer of #pi_{0}+#pi_{0} with Events / 0.02 GeV"};
		       		histCuts[id] = cutsToApply*p_tMassBinned[i];
			}
		}
	}
	
	//std::string rfCutString;
	// we will use i=0 to be the noCut graph. i>0 would be differnet CL cuts applied where i is the order.
	//for (int i = 0; i < 7; i++){
	//       if (i == 0){ cutsToApply = 1; cutsApplied="None"; applyAccSub=noAccSub;}
	//       if (i == 1){ cutsToApply = pCLKinFit1*pDiffCL; cutsApplied="pCLKinFit1*pDiffCL"; applyAccSub=weight;}
	//       if (i == 2){ cutsToApply = pCLKinFit*pDiffCL; cutsApplied="pCLKinFit*pDiffCL"; applyAccSub=weight;}
	//       if (i == 3){ cutsToApply = pCLKinFit3*pDiffCL; cutsApplied="pCLKinFit3*pDiffCL"; applyAccSub=weight;}
	//       if (i == 4){ cutsToApply = pCLKinFit4*pDiffCL; cutsApplied="pCLKinFit4*pDiffCL"; applyAccSub=weight;}
	//       if (i == 5){ cutsToApply = pCLKinFit5*pDiffCL; cutsApplied="pCLKinFit5*pDiffCL"; applyAccSub=weight;}
	//       if (i == 6){ cutsToApply = pCLKinFit6*pDiffCL; cutsApplied="pCLKinFit6*pDiffCL"; applyAccSub=weight;}
	//       rfCutString = "CutCLOrderNeg"+std::to_string(i);
	//       ++id; histVals[id] = {8, applyAccSub, locDeltaTRF}; 
	//       histList[id] = {("RFTime"+rfCutString).c_str(), ("Cuts="+cutsApplied+"; RF Time (ns)").c_str(), "400", "-20", "20", "Entries / 0.03 ns"};
	//	histCuts[id] = cutsToApply;
	//}
	
	if (outputUERegion){
		std::string ueCutString;
		// we will use i=0 to be the noCut graph. i>0 would be differnet CL cuts applied where i is the order.
		iUpUE=0.1;
		for (int i = 0; i < numRegions_UE; i++){
        	        ueCutString = std::to_string(iUpUE);
        	        // remove trailing 0's
        	        ueCutString.erase ( ueCutString.find_last_not_of('0') + 1, std::string::npos );
        	        iUpUE+=0.1;
			//reuse ueCutString to set the region name
		       cutsToApply = mEllipseUE_pre*chiSq100; cutsApplied=("mEllipseUE_pre & ChiSq100 & UE<"+ueCutString).c_str(); applyAccSub=weight;
			if (is_pi0eta){
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassUEregion"+std::to_string(i)).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("etaMassUEregion"+std::to_string(i)).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply;
			}
			else {
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassUEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("pi0MassUEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
			}
		}
	}

	if (outputPi0Resolutions){
		//pi0 resolutions
		std::string lowECutString;
		std::string upECutString;
		for (int i = 0; i < numRegions_E; i++){
        	        lowECutString = std::to_string(iLowE[i]);
        	        upECutString = std::to_string(iUpE[i]);
        	        // remove trailing 0's
        	        lowECutString.erase ( lowECutString.find_last_not_of('0') + 1, std::string::npos );
        	        upECutString.erase ( upECutString.find_last_not_of('0') + 1, std::string::npos );
			//reuse ueCutString to set the region name
			// ************* NOTE THAT mEllipse_pre SHOULD ONLY BE USED WHEN WEIGHT INCLUDES BKG_SUB RATHER mEllipse SHOULD BE USED 
		        cutsToApply = mEllipse_pre; cutsApplied=("mEllipse_pre & "+lowECutString+"<E(pi0)<"+upECutString+" & pSelectf2").c_str(); applyAccSub=weight;
			if (!is_pi0eta){
		        	++id; histVals[id] = {5, applyAccSub, theta_pi0_lab};
		        	histList[id] = {("thetaLabEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";#theta lab of #pi_{0}").c_str(), "100","0","70", "Events / 0.07"};
				// these are not the final cuts we apply! Its in the second section
		        	histCuts[id] = cutsToApply;
		        	++id; histVals[id] = {6, applyAccSub, theta_eta_lab};
		        	histList[id] = {("thetaLabEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";#theta lab of #pi_{0}").c_str(), "100","0","70", "Events / 0.07"};
				// these are not the final cuts we apply! Its in the second section
		        	histCuts[id] = cutsToApply;


		        	++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		        	histList[id] = {("pi0MassKinEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				// these are not the final cuts we apply! Its in the second section
		        	histCuts[id] = cutsToApply;
		        	++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		        	histList[id] = {("pi0MassEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				// these are not the final cuts we apply! Its in the second section
		        	histCuts[id] = cutsToApply;
		        	++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		        	histList[id] = {("pi0MassKinEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		        	histCuts[id] = cutsToApply;
		        	++id; histVals[id] = {6, applyAccSub, locEtaMass};
		        	histList[id] = {("pi0MassEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		        	histCuts[id] = cutsToApply;

		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassFCALKinEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InFCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassBCALKinEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InBCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("pi0MassFCALKinEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(),pi0BinRange[0], pi0BinRange[1], pi0BinRange[2] , "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInFCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("pi0MassBCALKinEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInBCAL; 

		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0MassFCALEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InFCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0MassBCALEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InBCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("pi0MassFCALEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(),pi0BinRange[0], pi0BinRange[1], pi0BinRange[2] , "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInFCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("pi0MassBCALEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInBCAL; 
			}
		}
	}

	if(outputChiSqRegion){	
		std::string chiSqCutString;
		iUpChiSq=5;
		// we will use i=0 to be the noCut graph. i>0 would be differnet CL cuts applied where i is the order.
		for (int i = 0; i < numRegions_ChiSq; i++){
        	        chiSqCutString = std::to_string(iUpChiSq);
        	        // remove trailing 0's
        	        chiSqCutString.erase ( chiSqCutString.find_last_not_of('0') + 1, std::string::npos );
        	        iUpChiSq+=5;
		       	cutsToApply = mEllipseUEChiSq_pre; cutsApplied=("mEllipseUEChiSq_pre & ChiSq<"+chiSqCutString).c_str(); applyAccSub=weight;
			if (is_pi0eta){
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassChiSqregion"+std::to_string(i)).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("etaMassChiSqregion"+std::to_string(i)).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply;
			}
			else {
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassChiSqregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("pi0MassChiSqregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
			}
		}
	}

	std::string pi0etaCutString;
	for (int i = 0; i < 5; i++){
	       // 0 = noCut; 1 = Gen Cuts + EBeam > 6; 2 = Gen Cuts + BaryBkg
	       if (i==0) { cutsToApply = 1; cutsApplied="None"; applyAccSub=noAccSub;}
	       if (i==1) { cutsToApply = allGeneralCutsPassed; cutsApplied="GeneralCuts"; applyAccSub=weight;}
	       if (i==2) { cutsToApply = allGeneralCutsPassed*pEtaProtonBaryonCut*ppi0ProtonBaryonCut; cutsApplied="GeneralCuts+BaryCuts"; applyAccSub=weight;}
	       if (i==3) { cutsToApply = allGeneralCutsPassed*pVanHove; cutsApplied="GeneralCuts+pVanHove"; applyAccSub=weight;}
		// If we were to compare allGeneralCutsPassed graph with the one just with bkg subtraction we should use mEllipse instead of mEllipse_pre. This is because the 
		// uniqueness tracking can do more bad things since the region with weight=0 is still contributing to the removal of some combinations. Some graphs to look at is in
		// ~/Desktop/Work/gluex/
	       if (i==4) { cutsToApply = mEllipse; cutsApplied="mEllipse"; applyAccSub=weight;}
	       pi0etaCutString = "NoCut_AccSubEBeam_AccSubBaryBkg"+std::to_string(i);
		if (is_pi0eta){
	       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
	       		histList[id] = {("pi0eta1D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
	       		histCuts[id] = cutsToApply;
		}
		else {
	       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
	       		histList[id] = {("pi0pi01D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}#pi_{0}) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
	       		histCuts[id] = cutsToApply;
		}
	}


	if(outputMassBinnedE){
		for (int i = 0; i < 3; i++){
		       // 0 = noCut; 1 = Gen Cuts + EBeam > 6; 2 = Gen Cuts + BaryBkg
		       if (i==0) { cutsToApply = 1; cutsApplied="None"; applyAccSub=noAccSub;}
		       // this section does not use pBeamE since this section bins the M(pi0eta) in beam energy. 
		       if (i==1) { cutsToApply = mBeamE; cutsApplied="mBeamE"; applyAccSub=weight;}
		       if (i==2) { cutsToApply = mBeamE*pEtaProtonBaryonCut*ppi0ProtonBaryonCut; cutsApplied="mBeamE*BaryCuts"; applyAccSub=weight;}
		       pi0etaCutString = std::to_string(i);
			if (is_pi0eta){
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0eta1D30to46ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*3.0<beamE<4.6 GeV;M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE30to46;
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0eta1D46to62ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*4.6<beamE<6.2 GeV;M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE46to62;
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0eta1D62to78ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*6.2<beamE<7.8 GeV;M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE62to78;
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0eta1D78to94ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*7.8<beamE<9.4 GeV;M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE78to94;
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0eta1D94to11ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*9.4<beamE<11 GeV;M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE94to11;
			}
			else {
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0pi01D30to46ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*3.0<beamE<4.6 GeV;M(#pi_{0}#pi_{0}) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE30to46;
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0pi01D46to62ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*4.6<beamE<6.2 GeV;M(#pi_{0}#pi_{0}) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE46to62;
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0pi01D62to78ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*6.2<beamE<7.8 GeV;M(#pi_{0}#pi_{0}) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE62to78;
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0pi01D78to94ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*7.8<beamE<9.4 GeV;M(#pi_{0}#pi_{0}) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE78to94;
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0pi01D94to11ID"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+"*9.4<beamE<11 GeV;M(#pi_{0}#pi_{0}) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply*pBeamE94to11;
			}
		}
	}
	
	if(output_mEllipse_UE_ChiSq){
		for (int i = 0; i<5; i++){
		       if (i == 0){ cutString = ""; cutsToApply = 1; cutsApplied="None"; applyAccSub=noAccSub;}
		       if (i == 1){ cutString = "Cut"; cutsToApply = mEllipse_pre; cutsApplied="mEllipse_pre"; applyAccSub=weight;}
		       if (i == 2){ cutString = "noUECut"; cutsToApply = mEllipseUE_pre; cutsApplied="mEllipseUE_pre"; applyAccSub=weight;}
		       if (i == 3){ cutString = "noUEChiSqCut"; cutsToApply = mEllipseUEChiSq_pre; cutsApplied="mEllipseUEChiSq_pre"; applyAccSub=weight;}
			// This one should not have any accsub since it would always basicalyl filly with a negative weight due to being in the BKG regional ellipse
		       if (i == 4){ cutString = "CutBKG"; cutsToApply = mEllipse*pYellowBKG; cutsApplied="mEllipsepYellowBKG"; applyAccSub=weightBS;}
			if(is_pi0eta){
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0Mass_Kin"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassFCAL_Kin"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InFCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassBCAL_Kin"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InBCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassSplit_Kin"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InSplit; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("etaMass_Kin"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("etaMassFCAL_Kin"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(),etaBinRange[0], etaBinRange[1], etaBinRange[2] , "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInFCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("etaMassBCAL_Kin"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInBCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("etaMassSplit_Kin"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInSplit; 
				// masses with the measured values
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0Mass"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0MassFCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InFCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0MassBCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InBCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0MassSplit"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InSplit; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("etaMass"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("etaMassFCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(),etaBinRange[0], etaBinRange[1], etaBinRange[2] , "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInFCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("etaMassBCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInBCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("etaMassSplit"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInSplit; 

				// 2d mass plot
		       		++id; histVals[id] = {4, applyAccSub, locPi0Mass_Kin, locEtaMass_Kin};
		       		histList[id] = {("pi0eta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], etaBinRange[0], etaBinRange[1], etaBinRange[2], "#pi_{0} Mass (GeV) with Events / 0.001 GeV","#eta Mass (GeV) with Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply;
			}
			else{
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0Mass_Kin_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassFCAL_Kin_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InFCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassBCAL_Kin_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InBCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
		       		histList[id] = {("pi0MassSplit_Kin_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InSplit; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("pi0Mass_Kin_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("pi0MassFCAL_Kin_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(),pi0BinRange[0], pi0BinRange[1], pi0BinRange[2] , "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInFCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("pi0MassBCAL_Kin_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInBCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
		       		histList[id] = {("pi0MassSplit_Kin_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInSplit; 
				// using kin fit to see differences in mass plots
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0Mass_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0MassFCAL_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InFCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0MassBCAL_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InBCAL; 
		       		++id; histVals[id] = {5, applyAccSub, locPi0Mass};
		       		histList[id] = {("pi0MassSplit_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InSplit; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("pi0Mass_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("pi0MassFCAL_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(),pi0BinRange[0], pi0BinRange[1], pi0BinRange[2] , "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInFCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("pi0MassBCAL_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInBCAL; 
		       		++id; histVals[id] = {6, applyAccSub, locEtaMass};
		       		histList[id] = {("pi0MassSplit_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInSplit; 
				// 2D mass plot
		       		++id; histVals[id] = {4, applyAccSub, locPi0Mass_Kin, locEtaMass_Kin};
		       		histList[id] = {("pi0pi0"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "#pi_{0} Mass (GeV) with Events / 0.001 GeV","#pi_{0} Mass (GeV) with Events / 0.001 GeV"};
		       		histCuts[id] = cutsToApply;
				// JUST FOR THE PIPI SYSTEM
				// LOOKING AT WRONG COMBINATIONS OF PAIRED PHOTONS WITHTHE KIN DATA
		       		++id; histVals[id] = {15, applyAccSub, locPi0Mass_Kin_mismatch};
		       		histList[id] = {("pi0Mass_KinMismatch_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {15, applyAccSub, locPi0Mass_Kin_mismatch};
		       		histList[id] = {("pi0MassFCAL_KinMismatch_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InFCAL_mismatch; 
		       		++id; histVals[id] = {15, applyAccSub, locPi0Mass_Kin_mismatch};
		       		histList[id] = {("pi0MassBCAL_KinMismatch_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InBCAL_mismatch; 
		       		++id; histVals[id] = {15, applyAccSub, locPi0Mass_Kin_mismatch};
		       		histList[id] = {("pi0MassSplit_KinMismatch_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pPi0InSplit_mismatch; 
		       		++id; histVals[id] = {16, applyAccSub, locEtaMass_Kin_mismatch};
		       		histList[id] = {("pi0Mass_KinMismatch_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {16, applyAccSub, locEtaMass_Kin_mismatch};
		       		histList[id] = {("pi0MassFCAL_KinMismatch_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(),etaBinRange[0], etaBinRange[1], etaBinRange[2] , "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInFCAL_mismatch; 
		       		++id; histVals[id] = {16, applyAccSub, locEtaMass_Kin_mismatch};
		       		histList[id] = {("pi0MassBCAL_KinMismatch_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInBCAL_mismatch; 
		       		++id; histVals[id] = {16, applyAccSub, locEtaMass_Kin_mismatch};
		       		histList[id] = {("pi0MassSplit_KinMismatch_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
		       		histCuts[id] = cutsToApply*pEtaInSplit_mismatch; 
			}
		}
	}


	
	if (outputBkgSubOnBkgAndSignalRegions){
		// ********** SECTION TO CHECK HOW THE BKG SUB REMOVES THE BKGS LIKE BARYON RESONANCES
		for (int i=0; i<7; ++i){
			// The effect of using applyAccSub=weight is basially like only selecting the red and yellow ellipses since the weight is = 0 when not in those regions.
			// Sig_weigthA should be the same as Sig_weightB in terms of hist shape but the number of events might differ since in Sig_weightA we keep the 0-weight events. This can be argued since the || opertor is a union of sets
			// 	and we can include the 0-weight region in for free since the weight=0 there. Recombine the 3 bkgSub regions to get the full set which is always true.
			// BUT!!!! THIS IS NOT THE CASE! THIS IS due to uniqueness tracking. Since the 0-weight region is so large we have a good chance of selecting a combination here and not being able to use that set of particles anymore.
		       if (i == 0){ pi0etaCutString = "Bkg_weightAS"; cutsToApply = mEllipse_pre*pYellowBKG; cutsApplied="mEllipse_pre*pYellowBKG"; applyAccSub=weightAS;}
		       if (i == 1){ pi0etaCutString = "Sig_weightAS"; cutsToApply = mEllipse_pre*pinsideEllipse; cutsApplied="mEllipse_pre*pinsideEllipse"; applyAccSub=weightAS;}
		       if (i == 2){ pi0etaCutString = "Sig_weight_badTracking"; cutsToApply = mEllipse_pre; cutsApplied="mEllipse_pre"; applyAccSub=weight;}
		       if (i == 3){ pi0etaCutString = "Sig_weight"; cutsToApply = mEllipse; cutsApplied="mEllipse"; applyAccSub=weight;}
		       if (i == 4){ pi0etaCutString = "Sig_weight_tLT1"; cutsToApply = mEllipse*ptLT1; cutsApplied="mEllipse*ptLT1"; applyAccSub=weight;}
		       if (i == 5){ pi0etaCutString = "Bkg_weightAS_tLT1"; cutsToApply = mEllipse_pre*pYellowBKG*ptLT1; cutsApplied="mEllipse_pre*pYellowBKG*ptLT1"; applyAccSub=weightAS;}
		       if (i == 6){ pi0etaCutString = "Sig_weightAS_tLT1"; cutsToApply = mEllipse_pre*pinsideEllipse*ptLT1; cutsApplied="mEllipse_pre*pinsideEllipse*ptLT1"; applyAccSub=weightAS;}
			if (is_pi0eta){
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0eta1D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
		       		histList[id] = {("pi0proton1D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin};
		       		histList[id] = {("etaproton1D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#etaproton) (GeV)").c_str(), "450", "0", "4.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, vanHove_x, vanHove_y};
		       		histList[id] = {("vanHove"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";Van Hove Plot").c_str(), "60","-3","3", "60", "-3", "3", "Events / 0.1 degrees", "Events / 0.1 degrees"};
		       		histCuts[id] = cutsToApply;
        	        	++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_eta_GJ};
        	        	histList[id] = {("eta_cosTheta_GJvsM"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta)Events / 0.01 GeV", "Cos(#theta) of #eta Events / 0.2"};
        	        	histCuts[id] = cutsToApply;
				++id; histVals[id] = {11, applyAccSub, mandelstam_tp};
				histList[id] = {("pi0eta_t"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using recoil+target").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {11, applyAccSub, mandelstam_tp_pe};
				histList[id] = {("pi0eta_t_pe"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using beam+pi0eta").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
				histCuts[id] = cutsToApply;
			}
			else{
		       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
		       		histList[id] = {("pi0pi01D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
		       		histList[id] = {("pi0proton1D_1"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin};
		       		histList[id] = {("pi0proton1D_2"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#etaproton) (GeV)").c_str(), "450", "0", "4.5", "Events / 0.01 GeV"};
		       		histCuts[id] = cutsToApply;
		       		++id; histVals[id] = {8, applyAccSub, vanHove_x, vanHove_y};
		       		histList[id] = {("vanHove"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";Van Hove Plot").c_str(), "60","-3","3", "60", "-3", "3", "Events / 0.1 degrees", "Events / 0.1 degrees"};
		       		histCuts[id] = cutsToApply;
        	        	++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_pi0_GJ};
        	        	histList[id] = {("pi0_cosTheta_GJvsM_1"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta)Events / 0.01 GeV", "Cos(#theta) of #eta Events / 0.2"};
        	        	histCuts[id] = cutsToApply;
        	        	++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_eta_GJ};
        	        	histList[id] = {("pi0_cosTheta_GJvsM_2"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta)Events / 0.01 GeV", "Cos(#theta) of #eta Events / 0.2"};
        	        	histCuts[id] = cutsToApply;
				++id; histVals[id] = {11, applyAccSub, mandelstam_tp};
				histList[id] = {("pi0pi0_t"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using recoil+target").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {11, applyAccSub, mandelstam_tp_pe};
				histList[id] = {("pi0pi0_t_pp"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using beam+pi0eta").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
				histCuts[id] = cutsToApply;
			}
		}
	}


	if(outputCorrelationBetweenMasses){
		// This section will look at the correlation between the M(pi0) vs M(pi0eta) and M(eta) vs M(pi0eta) to see if there exists strong correlation between the variables. This is to see if we
		// can properly use the SPlot method instead of bkg subtraction. SPlot method requires the control variable to not be correlated with the varaibles in the discriminating variable set.
		// we use mEllipse_pre since we want to include all the regions on the M(pi0)vsM(eta) 2D histogram; if we use mEllipse we would already be removing the 0-weight region in between the yellow
		// and red regions.
		if (is_pi0eta){
			cutsToApply = mEllipse_pre; cutString ="_correlation"; cutsApplied="mEllipse_pre"; applyAccSub=weightAS;
			++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin, locPi0Mass_Kin};
			histList[id] = {("pi0etaPi0"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin, locEtaMass_Kin};
			histList[id] = {("pi0etaEta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", etaBinRange[0], etaBinRange[1], etaBinRange[2], "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#eta_{0}) Events / 0.0025 GeV"};
			histCuts[id] = cutsToApply;
		}
	}

	if(outputMassShift){
		if (is_pi0eta){
			cutsToApply = allGeneralCutsPassed; cutString ="_chargedVertex"; cutsApplied="allGeneralCutsPassed"; applyAccSub=weightAS;
			++id; histVals[id] = {5, applyAccSub, locPi0Mass_charged};
			histList[id] = {("pi0Mass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {6, applyAccSub, locEtaMass_charged};
			histList[id] = {("etaMass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
		       	++id; histVals[id] = {4, applyAccSub, locPi0Mass_charged, locEtaMass_charged};
		       	histList[id] = {("pi0eta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], etaBinRange[0], etaBinRange[1], etaBinRange[2], "#pi_{0} Mass (GeV) with Events / 0.001 GeV","#eta Mass (GeV) with Events / 0.0025 GeV"};
		       	histCuts[id] = cutsToApply;

			cutsToApply = allGeneralCutsPassed; cutString ="_targetVertex"; cutsApplied="allGeneralCutsPassed"; applyAccSub=weightAS;
			++id; histVals[id] = {5, applyAccSub, locPi0Mass_target};
			histList[id] = {("pi0Mass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {6, applyAccSub, locEtaMass_target};
			histList[id] = {("etaMass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
		       	++id; histVals[id] = {4, applyAccSub, locPi0Mass_target, locEtaMass_target};
		       	histList[id] = {("pi0eta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], etaBinRange[0], etaBinRange[1], etaBinRange[2], "#pi_{0} Mass (GeV) with Events / 0.001 GeV","#eta Mass (GeV) with Events / 0.0025 GeV"};
		       	histCuts[id] = cutsToApply;
		}
		else {
			cutsToApply = allGeneralCutsPassed; cutString ="_chargedVertex"; cutsApplied="allGeneralCutsPassed"; applyAccSub=weightAS;
			++id; histVals[id] = {5, applyAccSub, locPi0Mass_charged};
			histList[id] = {("pi0Mass"+cutString+"_1").c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
			cutsToApply = allGeneralCutsPassed; cutString ="_chargedVertex"; cutsApplied="allGeneralCutsPassed"; applyAccSub=weightAS;
			++id; histVals[id] = {6, applyAccSub, locEtaMass_charged};
			histList[id] = {("pi0Mass"+cutString+"_2").c_str(), ("Cuts="+cutsApplied).c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
		       	++id; histVals[id] = {4, applyAccSub, locPi0Mass_charged, locEtaMass_charged};
		       	histList[id] = {("pi0pi0"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], etaBinRange[0], etaBinRange[1], etaBinRange[2], "#pi_{0} Mass (GeV) with Events / 0.001 GeV","#eta Mass (GeV) with Events / 0.0025 GeV"};
		       	histCuts[id] = cutsToApply;

			cutsToApply = allGeneralCutsPassed; cutString ="_targetVertex"; cutsApplied="allGeneralCutsPassed"; applyAccSub=weightAS;
			++id; histVals[id] = {5, applyAccSub, locPi0Mass_target};
			histList[id] = {("pi0Mass"+cutString+"_1").c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
			cutsToApply = allGeneralCutsPassed; cutString ="_targetVertex"; cutsApplied="allGeneralCutsPassed"; applyAccSub=weightAS;
			++id; histVals[id] = {6, applyAccSub, locEtaMass_target};
			histList[id] = {("pi0Mass"+cutString+"_2").c_str(), ("Cuts="+cutsApplied).c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
		       	++id; histVals[id] = {4, applyAccSub, locPi0Mass_target, locEtaMass_target};
		       	histList[id] = {("pi0pi0"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], etaBinRange[0], etaBinRange[1], etaBinRange[2], "#pi_{0} Mass (GeV) with Events / 0.001 GeV","#eta Mass (GeV) with Events / 0.0025 GeV"};
		       	histCuts[id] = cutsToApply;
		}
	}




	//if (is_pi0eta){ 
	//	 mEllipse is used here since we are trying to show how the bkg and signal regions from the bkgSub technique changes things. 
	//	cutsToApply = mEllipse*!pMPi0P14; cutsApplied="mEllipse+!pMPi0P14"; applyAccSub=weight; cutString="_Sig_weight_MPi0PGT14";
        //        ++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_eta_GJ};
        //        histList[id] = {("eta_cosTheta_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta)Events / 0.01 GeV", "Cos(#theta) of #eta Events / 0.2"};
        //        histCuts[id] = cutsToApply;
	//	cutsToApply = mEllipse*pMPi0P14; cutsApplied="mEllipse+pMPi0P14"; applyAccSub=weight; cutString="_Sig_weight_MPi0PLT14";
	//	++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
	//	histList[id] = {("pi0eta1D_Mppi014"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
	//	histCuts[id] = cutsToApply;
	//	++id; histVals[id] = {11, applyAccSub, mandelstam_tp};
	//	histList[id] = {("pi0eta_t_Mppi014"+cutString).c_str(), ("Cuts="+cutsApplied+";t momentum transfer of #pi_{0}+#eta").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
	//	histCuts[id] = cutsToApply;
	//}
	// ******************** END SECTION ON BKG SUBTRACTION


	cutVariations = {"", "Cut", "CutOnlydzR", "CutOnlydzRdEdx", "CutOnlyMMChiSq"};
	for (int j=0; j<5; ++j){
	       if (j==0) { cutsToApply = 1; cutsApplied="None"; applyAccSub=noAccSub;}
	       if (j==1) { cutsToApply = mMagP3; cutsApplied="mMagP3"; applyAccSub=weight;}
	       if (j==2) { cutsToApply = pzCutmin*pRProton; cutsApplied="pzCutmin*pRProton"; applyAccSub=weight;}
	       if (j==3) { cutsToApply = pzCutmin*pRProton*pdEdxCDCProton; cutsApplied="pzCutmin*pRProton*pdEdxCDCProton"; applyAccSub=weight;}
	       if (j==4) { cutsToApply = pMissingMassSquared*pChiSq; cutsApplied="pMissingMassSquared*pChiSq"; applyAccSub=weight;}
	       cutString = cutVariations[j];
	       ++id; histVals[id] = {2, applyAccSub, locMagP3Proton}; 
	       histList[id] = {("MagP3Proton"+cutString).c_str(), ("Cuts="+cutsApplied+";Proton Momenetum GeV").c_str(), "200", "0" , "5", "Events / 0.025 GeV"};
	       histCuts[id] = cutsToApply;
	}

	if (outputPhotonXandShowerLoc){
		cutsApplied="None";
		cutsToApply = 1;
		applyAccSub = noAccSub;
		++id; histVals[id] = {13, applyAccSub, photonXs_Kin[0]}; 
		histList[id] = {"X_Kin",("Cuts="+cutsApplied+";position (cm)").c_str(), "160", "0", "160", "Events / 1 cm"};
		histCuts[id] = cutsToApply;
		histCuts[id] = cutsToApply;
		++id; histVals[id] = {13, applyAccSub, photonYs_Kin[0]}; 
		histList[id] = {"Y_Kin",("Cuts="+cutsApplied+";position (cm)").c_str(), "160", "0", "160", "Events / 1 cm"};
		histCuts[id] = cutsToApply;
		++id; histVals[id] = {13, applyAccSub, photonZs_Kin[0]}; 
		histList[id] = {"Z_Kin",("Cuts="+cutsApplied+";position (cm)").c_str(), "250", "0", "1000", "Events / 4 cm"};
		histCuts[id] = cutsToApply;
		++id; histVals[id] = {13, applyAccSub, photonTs_Kin[0]}; 
		histList[id] = {"T_Kin",("Cuts="+cutsApplied+";time (s)").c_str(), "160", "0", "160", "Events / 1 cm"};
		histCuts[id] = cutsToApply;
		++id; histVals[id] = {13, applyAccSub, photonXs_Shower[0]}; 
		histList[id] = {"X_Shower",("Cuts="+cutsApplied+";position (cm)").c_str(), "160", "0", "160", "Events / 1 cm"};
		histCuts[id] = cutsToApply;
		++id; histVals[id] = {13, applyAccSub, photonYs_Shower[0]}; 
		histList[id] = {"Y_Shower",("Cuts="+cutsApplied+";position (cm)").c_str(), "160", "0", "160", "Events / 1 cm"};
		histCuts[id] = cutsToApply;
		++id; histVals[id] = {13, applyAccSub, photonZs_Shower[0]}; 
		histList[id] = {"Z_Shower",("Cuts="+cutsApplied+";position (cm)").c_str(), "160", "0", "160", "Events / 1 cm"};
		histCuts[id] = cutsToApply;
		++id; histVals[id] = {13, applyAccSub, photonTs_Shower[0]}; 
		histList[id] = {"T_Shower",("Cuts="+cutsApplied+";time (s)").c_str(), "160", "0", "160", "Events / 1 cm"};
		histCuts[id] = cutsToApply;
	}



	// ***** MUST CHANGE groupNames INSIDE THE HEADER FILE *****************
	//std::string groupNames[16] = {"beam", "p1", "ph12b1", "ph1234", "ph12", "ph34", "ph34b1", "entire combo", "ph12p1", "ph34p1", "ph1234p1", "FCAL Pairs", "phN", "BCAL Pairs", "ph13", "ph24"};
	if(showOutput) {cout << endl << "Uniqueness tracking: Let photon 1-4, proton 1, 1 beam photon = {ph1,ph2,ph3,ph4,p1,b1}" << endl << "----------------------- Tracking prototypes ---------------------------" << endl;}
	if(showOutput) {cout <<  "Tracking group 1: beam will track only beam photon" << endl; }
	if(showOutput) {cout <<  "Tracking group 2: proton will track {p1}" << endl; }
	if(showOutput) {cout <<  "Tracking group 3: pi0_cm will track {ph1,ph2,b1} - including b1 since Center of Mass uses photonbeam and target" << endl; }
	if(showOutput) {cout <<  "Tracking group 4: pi0eta will track {ph1,ph2,ph3,ph4}" << endl; }
	if(showOutput) {cout <<  "Tracking group 5: pi0 will track {ph1,ph2}" << endl; }
	if(showOutput) {cout <<  "Tracking group 6: eta will track {ph3,ph4}" << endl; }
	if(showOutput) {cout <<  "Tracking group 7: eta_cm will track {ph3,ph4,b1} - include b1 since Center of Mass uses photonbeam and target" << endl; }
	if(showOutput) {cout <<  "Tracking group 8: MissingMassSquared will track entire combo - will also track timing which relates to the RFtime and thus the entire combo" << endl; }
	if(showOutput) {cout <<  "Tracking group 9: pi0proton will track {ph1,ph2,p1}" << endl; }
	if(showOutput) {cout <<  "Tracking group 10: etaproton will track {ph3,ph4,p1}" << endl; }
	if(showOutput) {cout <<  "Tracking group 11: pi0etaproton will track {ph1,ph2,ph3,ph4,p1}" << endl; }
	if(showOutput) {cout <<  "Tracking group 12: (FCAL) pairwise photon will track {phN,phM} agnostic to which particle number only PID" << endl; }
	if(showOutput) {cout <<  "Tracking group 13: photon will track {phN} agnostic to which particle number only PID" << endl; }
	if(showOutput) {cout <<  "Tracking group 14: (BCAL) pairwise photon will track {phN,phM} agnostic to which particle number only PID" << endl; }
	if(showOutput) {cout <<  "Tracking group 15: pi0Mismatch will track {ph1,ph3}" << endl; }
	if(showOutput) {cout <<  "Tracking group 16: etaMismatch will track {ph2,ph4}" << endl; }
	if(showOutput) {cout <<  "Tracking group 17: No uniqueness tracking check" << endl; }
	int totalGroups=17;

	// Split up the ids into groups so we do not have to loop through all the histograms when trying to do our filling.
	for(int groupNum=1; groupNum<totalGroups+1; ++groupNum){
		group_ids.clear();
		// loop through all the histograms
		for(int id_idx=0; id_idx<id+1; ++id_idx){
			if(histVals[id_idx][0]==groupNum){
				group_ids.push_back(id_idx);
				//if(showOutput){ cout << histList[id_idx][0] << " is under the group: " << std::to_string(groupNum) << endl; }
			}
		}
		vec_group_ids.push_back(group_ids);
	}


	if(onlyNamesPi0_1){ cout << "ONLY HISTOGRAMS WITH AFFIX _1 and _2 FOR MERGING PI0_1 AND PI0_2 HISTOGRAMS !!!" << endl; }
	int hist_dim;
	for(std::size_t i = 0; i<vec_group_ids.size(); ++i){
		if((!onlyNamesPi0_1)*showOutput) {cout << " hist Ids for a group: " << std::to_string(i+1) << " - " << groupNames[i] << endl << "=============================" << endl;}
		for(std::size_t j = 0; j<vec_group_ids[i].size(); ++j){
			if(onlyNamesPi0_1){
				if(histList[vec_group_ids[i][j]][0].find("_1") != std::string::npos) {
					if ( histList[vec_group_ids[i][j]].size() > 8 ){ hist_dim = 2; }
					else { hist_dim = 1;}
					if(showOutput){cout << histList[vec_group_ids[i][j]][0] << "," << hist_dim << endl;}
				}
			}
			else {
				if(showOutput){cout << histList[vec_group_ids[i][j]][0] << " with ID:" << histVals[vec_group_ids[i][j]][0] << endl;}
			}
		}
	}

	if(showOutput) { cout << endl << "-------------------------------" << endl << "Total Number of histograms: " << std::to_string(id+1) << endl << "----------------------------------" << endl;}

// **************************************** MAKE THE HISTOGRAMS ********************************************//
	for(int i=0; i<id+1; ++i){// need + 1 since we have to include the last histogram.... 
		if(showOutput) {cout << "Building histogram: " << histList[i][0] << " with id = " << std::to_string(i) << " and group = " << std::to_string(histVals[i][0]) << endl; }
		int histListSize = (int)histList[i].size();
		if(histListSize == 6 || histListSize == 7 ){
			dHist_all1DHists[i] = new TH1F(histList[i][0].c_str(), histList[i][1].c_str(), atof(histList[i][2].c_str()), atof(histList[i][3].c_str()), atof(histList[i][4].c_str()));
			dHist_all1DHists[i]->SetYTitle(histList[i][5].c_str());
			// some have a title making the size 1 larger
			if(histListSize==7){
				dHist_all1DHists[i]->SetTitle(histList[i][6].c_str());
			}
			
		}
		else if (histListSize==10 || histListSize==11){
			dHist_all2DHists[i] = new TH2F(histList[i][0].c_str(), histList[i][1].c_str(), atof(histList[i][2].c_str()), atof(histList[i][3].c_str()), atof(histList[i][4].c_str()), 
				atof(histList[i][5].c_str()), atof(histList[i][6].c_str()), atof(histList[i][7].c_str()));
			dHist_all2DHists[i]->SetXTitle(histList[i][8].c_str());
			dHist_all2DHists[i]->SetYTitle(histList[i][9].c_str());
			// it could exist that a 2D hist has a title making total num of params = 11
			if((int)histList[i].size()==11){
				dHist_all2DHists[i]->SetTitle(histList[i][10].c_str());
			}
		}
		else { cout << "ERROR!" << endl  << "ERROR!" << endl << "ERROR!" << endl; }
	}


        if(showOutput){cout << endl << endl << "--------------------------------" << endl << "Initialized all histograms"  <<  endl;}


	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
	*/

        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("AccWeight"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("uniqueComboID"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mandelstam_tp"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("finalStateComboID"); //fundamental = char, int, float, double, etc.
	if (is_pi0eta){
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Meta"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0eta"); //fundamental = char, int, float, double, etc.
	}
	else{
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0pi0"); //fundamental = char, int, float, double, etc.
	}
	// Introduce some angles to use in the phase space distance calculation
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("cosTheta_X_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phi_X_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("cosTheta_eta_gj"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phi_eta_gj"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("cosThetaHighestEphotonIneta_gj"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("cosThetaHighestEphotonInpi0_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("vanHove_x");
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("vanHove_y");
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("vanHove_omega");
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("pi0_energy");	
        //dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiHighestEphotonIneta_etaFrame"); //fundamental = char, int, float, double, etc.
        //dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("cosThetaHighestEphotonInpi0_pi0Frame"); //fundamental = char, int, float, double, etc.
        //dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phiHighestEphotonInpi0_pi0Frame"); //fundamental = char, int, float, double, etc.
        //dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("HighestEphotonIneta_etaFrame"); //fundamental = char, int, float, double, etc.
        //dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("HighestEphotonIneta_pi0Frame"); //fundamental = char, int, float, double, etc.

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
	



} // end of initialization

Bool_t DSelector_ver20::Process(Long64_t locEntry)
{

	//if(itersToRun<5){ ++itersToRun; //so we can just try to show the outut of one event 
        if(showOutput){cout << "Starting next process looping" << endl;}
	// The Process() function is called for each entry in the tree. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	// Use fStatus to set the return value of TTree::Process().
	// The return value is currently not used.

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();




	// NOT SURE IF THIS WOULD EVEN WORK. SO IT SEEMS LIKE dThrownWrapper MATCHES EACH OF THE THROWN PARTICLES TO A HYPOTHESIS AND ASSIGNS A FIGURE OF MERIT TO IT.
	// THERE IS NO WAY TO REALLY TELL FROM THIS HOW TO SEPARATE SIGNAL AND BKG UNLESS WE JUST TAKE THE HIGHEST FOM COMBINATION AND USE IT, DISCARDING ALL OTHER COMBOS.
	// PROBABLY NOT SO GOOD OF A WAY OF DISCERNING BKG AND SIGNAL. 
	// First we loop through the thrown combo wrapper and save all the parent ids and p4s
	//cout << "Event " << eventIdx << "!\n-----------------------" << endl;
	//TLorentzVector thrownP4s[Get_NumThrown()];
	//Int_t thrownPIDs[Get_NumThrown()];
	//for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	//{
	//	//Set branch array indices corresponding to this particle
	//	dThrownWrapper->Set_ArrayIndex(loc_i);
	//	thrownP4s[loc_i] = dThrownWrapper->Get_P4();
	//	thrownPIDs[loc_i] = dThrownWrapper->Get_PID();
	//	
	//}
	//for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	//{
	//	dThrownWrapper->Set_ArrayIndex(loc_i);
	//	cout << "particle_pid: " <<  thrownPIDs[loc_i] << " with parent_pid: " << thrownPIDs[dThrownWrapper->Get_ParentIndex()] << endl;
	//	cout << "	matching NeutralID/TrackID: " << dThrownWrapper->Get_MatchID() << endl;
	//	cout << "	FOM: " << dThrownWrapper->Get_MatchFOM() << endl;
	//}
	//cout << endl;



	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
        // we must have the following condition or else we run into errors, or maybe rcdb hangs due to the amount of querries or something. Anyways, it wouldn't work. 
        if(locRunNumber != dPreviousRunNumber)
        {
                dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
                hasPolarizationAngle = dAnalysisUtilities.Get_PolarizationAngle(locRunNumber, locPolarizationAngle);
        	if(showOutput){cout << "Getting beam polarization and filling used runs" << endl;}
                if (usedRuns.find(locRunNumber)==usedRuns.end()){
                        if (hasPolarizationAngle) {
                                dHist_BeamAngle->Fill(locPolarizationAngle);
                        }
                        usedRuns.insert(locRunNumber);
                }
        dPreviousRunNumber = locRunNumber;
        }


	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/
	// For uniqueness tracking we have to make a separate uniqueness tracking for each combination. This is partly because our general cuts span numerous quantities that can come from individual particles
	// from particle pairs, and from the entire combo. This makes it difficult to not apply the pi0 cut when filling the M(pi0) hist since two combos can contain the same 2 photons yet one can pass and the
	// other not (I think this could happen). If it can happen then the order in which the algorithm encouters the combo can fill the plot one or twice due to the other particles in the combo affecting the cut.
	// So we have to make filling the histogram and inserting into the uniquess set based on the same condition. 
	// There are 4Gamma + 1Proton + 1Beam in a combo. 

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();
	//dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	bool fSetNotEmpty = false;
	set<Int_t> usedBeam[id+1]; //Int_t: Unique ID for beam particles. set: easy to use, fast to search
	set<Int_t> usedP1[id+1];
	set<Int_t> usedPhN[id+1];

	set<map<Particle_t, set<Int_t> > > usedPh1234[id+1]; 
	set<map<Particle_t, set<Int_t> > > usedPh12[id+1];
	set<map<Particle_t, set<Int_t> > > usedPh34[id+1];
	set<map<Particle_t, set<Int_t> > > usedPh13[id+1];
	set<map<Particle_t, set<Int_t> > > usedPh24[id+1];
	set<map<Particle_t, set<Int_t> > > usedPairIdsFCAL[id+1]; // we dont make an array of them since there are only a few
	set<map<Particle_t, set<Int_t> > > usedPairIdsBCAL[id+1];
	//EXmap<Particle_t, AMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
        //use set<map .... > > even for reactions with the only one type of final state particles as in the case of 4 gammas
	set<map<Particle_t, set<Int_t> > > usedCombo[id+1];
	set<map<Particle_t, set<Int_t> > > usedPh12P1[id+1];
	set<map<Particle_t, set<Int_t> > > usedPh34P1[id+1];
	set<map<Particle_t, set<Int_t> > > usedPh1234P1[id+1];
	set<map<Particle_t, set<Int_t> > > usedPh12B1[id+1];
	set<map<Particle_t, set<Int_t> > > usedPh34B1[id+1];
	// MAKING SURE ALL THE SETS ARE EMPTY SO WE DONT START A LOOP WITH ELEMENTS IN IT ALREADY...
	for(int it3=0; it3<id+1; ++it3){
		if(usedBeam[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedP1[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPhN[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPh1234[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPh12[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPh34[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPairIdsFCAL[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPairIdsBCAL[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedCombo[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPh12P1[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPh34P1[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPh1234P1[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPh12B1[it3].size()!=0){fSetNotEmpty = 1; break;}
		if(usedPh34B1[it3].size()!=0){fSetNotEmpty = 1; break;}
	}
	if(showOutput) { 
		if (fSetNotEmpty) { cout << endl << endl << "******************************************" << endl << "ERROR - fSetNotEmpty!" << endl << "*********************************************" << endl;}
		else { cout << endl << endl << "******************************************" << endl << "Good all uniqueness sets are empty and initalized!\nStarting to loop through combos!" << endl << "*********************************************" << endl;}
	}

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

	/************************************************* LOOP OVER COMBOS *************************************************/


	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
                if(showOutput) { cout << "\n\n\n*********************************************************\n**************************************************\n########    EventIdx, ComboIdx: " << (eventIdx) << ", " << loc_i << "    #############" << endl; }
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
                if(showOutput){cout << "** Getting tracks p4, x4, ids" << endl;}
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		if(showOutput){cout <<"Got beamID"<<endl;}
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();
		if(showOutput){cout <<"Got ProtonID"<<endl;}

		//Step 1
		Int_t locPhoton1NeutralID = dPhoton1Wrapper->Get_NeutralID();
		Int_t locPhoton2NeutralID = dPhoton2Wrapper->Get_NeutralID();
		if(showOutput){cout <<"Got photon1/2 ID"<<endl;}

		//Step 2
		Int_t locPhoton3NeutralID = dPhoton3Wrapper->Get_NeutralID();
		Int_t locPhoton4NeutralID = dPhoton4Wrapper->Get_NeutralID();
		if(showOutput){cout <<"Got photon3/4 ID"<<endl;}
        	if(showOutput){cout << "Got Beam, proton, 4 photon IDs" << endl;}

// ********************************************************************************************************************************
		// Get 4 Momenta!!! 
// ********************************************************************************************************************************

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		if(showOutput){cout << "** Getting kin fitted X4, P4" << endl;}
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4_Measured();
                TLorentzVector locBeamX4 = dComboBeamWrapper->Get_X4_Measured();
        	//if(showOutput){cout << "Got kin fitted beam x4 P4" << endl;}
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4_Measured();
		//Step 1
		// using dDecayingPi0/EtaWrapper->Get_P4 no longer works, gives me an error now, luckily we don't use it anyways
                //TLorentzVector locDecayingPi0P4 = dDecayingPi0Wrapper->Get_P4();
                TLorentzVector locPhoton1X4_Kin = dPhoton1Wrapper->Get_X4();
                TLorentzVector locPhoton2X4_Kin = dPhoton2Wrapper->Get_X4();
        	if(showOutput){cout << "Got kin fitted photon 1,2 x4 and P4" << endl;}
		//Step 2
                //TLorentzVector locDecayingEtaP4 = dDecayingEtaWrapper->Get_P4();
                TLorentzVector locPhoton3X4_Kin = dPhoton3Wrapper->Get_X4();
                TLorentzVector locPhoton4X4_Kin = dPhoton4Wrapper->Get_X4();
        	if(showOutput){cout << "Got kin fitted photon 3,4 x4 and P4" << endl;}


		// Get Measured P4's:
		//Step 0
        	if(showOutput){cout << "** Getting measured x4 p4" << endl;}
		TLorentzVector locBeamP4_Kin = dComboBeamWrapper->Get_P4();
                TLorentzVector locBeamX4_Kin = dComboBeamWrapper->Get_X4();
		TLorentzVector locProtonP4_Kin = dProtonWrapper->Get_P4();
                TLorentzVector locProtonX4_Kin = dProtonWrapper->Get_X4();
		//Step 1
		TLorentzVector locPhoton1P4_Kin = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4_Kin = dPhoton2Wrapper->Get_P4();
		TLorentzVector locPhoton1P4 = dPhoton1Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton2P4 = dPhoton2Wrapper->Get_P4_Measured();
		//Step 2
		TLorentzVector locPhoton3P4_Kin = dPhoton3Wrapper->Get_P4();
		TLorentzVector locPhoton4P4_Kin = dPhoton4Wrapper->Get_P4();
		TLorentzVector locPhoton3P4 = dPhoton3Wrapper->Get_P4_Measured();
		TLorentzVector locPhoton4P4 = dPhoton4Wrapper->Get_P4_Measured();

		// We need to use shower position since the photons do not really have a track. 
                TLorentzVector locPhoton1X4_Shower = dPhoton1Wrapper->Get_X4_Shower();
                TLorentzVector locPhoton2X4_Shower = dPhoton2Wrapper->Get_X4_Shower();
                TLorentzVector locPhoton3X4_Shower = dPhoton3Wrapper->Get_X4_Shower();
                TLorentzVector locPhoton4X4_Shower = dPhoton4Wrapper->Get_X4_Shower();

        	if(showOutput){cout << "Got measured x4 p4 with using shower X4 for photons" << endl;}



		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

        	if(showOutput){cout << "** Combining 4-vectors to get pi0, eta, pi0eta, etaProton, pi0Proton" << endl;}
		// Combine 4-vectors
		TLorentzVector locMissingP4 = locBeamP4 + dTargetP4;
		locMissingP4 -= locProtonP4 + locPhoton1P4 + locPhoton2P4 + locPhoton3P4 + locPhoton4P4;

                //TLorentzVector locPi0P4 = locPhoton1P4 + locPhoton2P4;
                //TLorentzVector locEtaP4 = locPhoton3P4 + locPhoton4P4;

                TLorentzVector locPi0P4_Kin = locPhoton1P4_Kin + locPhoton2P4_Kin;
                TLorentzVector locEtaP4_Kin = locPhoton3P4_Kin + locPhoton4P4_Kin;
                TLorentzVector locPi0P4 = locPhoton1P4 + locPhoton2P4;
                TLorentzVector locEtaP4 = locPhoton3P4 + locPhoton4P4;
		
                TLorentzVector locPi0P4_Kin_mismatch = locPhoton1P4_Kin + locPhoton3P4_Kin;
                TLorentzVector locEtaP4_Kin_mismatch = locPhoton2P4_Kin + locPhoton4P4_Kin;

                TLorentzVector mixingPi0Eta_Kin = locPi0P4_Kin + locEtaP4_Kin;
                TLorentzVector mixingPi0Eta = locPi0P4 + locEtaP4;
                TLorentzVector mixingEtaProton_Kin = locEtaP4_Kin+locProtonP4_Kin;
                TLorentzVector mixingPi0Proton_Kin = locPi0P4_Kin+locProtonP4_Kin;

        	if(showOutput){cout << "Combined 4-vectors" << endl;}

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		
		// Perform_Action gives an error  with our current code probably because it does not have any actions to perform
		//dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
		*/


// ********************************************************************************************************************************
		// Calculating some variables!!! 
// ********************************************************************************************************************************
        	if(showOutput){cout << "\nCalculating variables to fill the histograms with and to put into histVals..." << endl;}
		        //Missing Mass Squared
		locMissingMassSquared = locMissingP4.M2();
                // 4 Momenta related variables  
                locBeamE = locBeamP4_Kin.E();

                locEtaMass_Kin = locEtaP4_Kin.M();
                //double locEtaMass_KinFit = locEtaP4.M();
                locPi0Mass_Kin = locPi0P4_Kin.M();
                //double locPi0Mass_KinFit = locPi0P4.M();

                locEtaMass = locEtaP4.M();
                locPi0Mass = locPi0P4.M();

                locEtaMass_Kin_mismatch = locEtaP4_Kin_mismatch.M();
                locPi0Mass_Kin_mismatch = locPi0P4_Kin_mismatch.M();

                locEtaE_Kin = locEtaP4_Kin.E();
                locPi0E_Kin = locPi0P4_Kin.E();

                locEtaProton_Kin = mixingEtaProton_Kin.M();
                locPi0Proton_Kin = mixingPi0Proton_Kin.M();
                locPi0Eta_Kin = mixingPi0Eta_Kin.M();

		// IN THE FOLLOWING SECTION WE WILL CALCUALTE BY OURSELVES THE MASS OF THE PI0 AND THE ETA WITH DIFFERENT STARTING POINTS ( USING THE PROTON X3 VS USING THE TARGET CENTER)
		//
		TLorentzVector anyVertexMomentum1;
		TLorentzVector anyVertexMomentum2;
		TVector3 vecFromVertexToShower;

		// Calculating the pi0 mass using the charged proton x3 as the initial starting point of the photon
		vecFromVertexToShower = locPhoton1X4_Shower.Vect()-locProtonX4_Kin.Vect();
		vecFromVertexToShower *= locPhoton1P4.E()/vecFromVertexToShower.Mag();
		anyVertexMomentum1.SetVect(vecFromVertexToShower);
		anyVertexMomentum1.SetE(locPhoton1P4.E());
		vecFromVertexToShower = locPhoton2X4_Shower.Vect()-locProtonX4_Kin.Vect();
		vecFromVertexToShower *= locPhoton2P4.E()/vecFromVertexToShower.Mag();
		anyVertexMomentum2.SetVect(vecFromVertexToShower);
		anyVertexMomentum2.SetE(locPhoton2P4.E());
		locPi0Mass_charged = (anyVertexMomentum1+anyVertexMomentum2).M();
		
		// Calculating the eta mass using the charged proton x3 as the initial starting point of the photon
		vecFromVertexToShower = locPhoton3X4_Shower.Vect()-locProtonX4_Kin.Vect();
		vecFromVertexToShower *= locPhoton3P4.E()/vecFromVertexToShower.Mag();
		anyVertexMomentum1.SetVect(vecFromVertexToShower);
		anyVertexMomentum1.SetE(locPhoton3P4.E());
		vecFromVertexToShower = locPhoton4X4_Shower.Vect()-locProtonX4_Kin.Vect();
		vecFromVertexToShower *= locPhoton4P4.E()/vecFromVertexToShower.Mag();
		anyVertexMomentum2.SetVect(vecFromVertexToShower);
		anyVertexMomentum2.SetE(locPhoton4P4.E());
		locEtaMass_charged = (anyVertexMomentum1+anyVertexMomentum2).M();

		// Calculating the pi0 mass using the target center x3 as the initial starting point of the photon
		vecFromVertexToShower = locPhoton1X4_Shower.Vect()-targetCenter;
		vecFromVertexToShower *= locPhoton1P4.E()/vecFromVertexToShower.Mag();
		anyVertexMomentum1.SetVect(vecFromVertexToShower);
		anyVertexMomentum1.SetE(locPhoton1P4.E());
		vecFromVertexToShower = locPhoton2X4_Shower.Vect()-targetCenter;
		vecFromVertexToShower *= locPhoton2P4.E()/vecFromVertexToShower.Mag();
		anyVertexMomentum2.SetVect(vecFromVertexToShower);
		anyVertexMomentum2.SetE(locPhoton2P4.E());
		locPi0Mass_target = (anyVertexMomentum1+anyVertexMomentum2).M();
		 
		// Calculating the eta mass using the target center x3 as the initial starting point of the photon
		vecFromVertexToShower = locPhoton3X4_Shower.Vect()-targetCenter;
		vecFromVertexToShower *= locPhoton3P4.E()/vecFromVertexToShower.Mag();
		anyVertexMomentum1.SetVect(vecFromVertexToShower);
		anyVertexMomentum1.SetE(locPhoton3P4.E());
		vecFromVertexToShower = locPhoton4X4_Shower.Vect()-targetCenter;
		vecFromVertexToShower *= locPhoton4P4.E()/vecFromVertexToShower.Mag();
		anyVertexMomentum2.SetVect(vecFromVertexToShower);
		anyVertexMomentum2.SetE(locPhoton4P4.E());
		locEtaMass_target = (anyVertexMomentum1+anyVertexMomentum2).M();

		if(showMassCalc) {
			cout << "Mass of pi0 for labeled=" << locPi0Mass << "  charged=" << locPi0Mass_charged << "  target=" << locPi0Mass_target << endl;
			cout << "Mass of eta for labeled=" << locEtaMass << "  charged=" << locEtaMass_charged << "  target=" << locEtaMass_target << endl;
			cout << endl;
		}	



                // Event Selction Variables
                //Charged Track
                locPtProton = TMath::Sqrt(TMath::Power(locProtonP4_Kin.Px(),2)+TMath::Power(locProtonP4_Kin.Py(),2));
                locPzProton = locProtonP4_Kin.Pz();
                locPolarAngleProton = locProtonP4_Kin.Theta()*radToDeg;
                locXProton = locProtonX4_Kin.X();
                locYProton = locProtonX4_Kin.Y();
                locRProton = TMath::Sqrt(TMath::Power(locXProton,2)+TMath::Power(locYProton,2));
                locdzProton = locProtonX4_Kin.Z();
                locdEdxCDCProton = dProtonWrapper->Get_dEdx_CDC();
                locdEdxFDCProton = dProtonWrapper->Get_dEdx_FDC();
                locMagP3Proton = TMath::Sqrt(TMath::Power(locProtonP4_Kin.Px(),2)+TMath::Power(locProtonP4_Kin.Py(),2)+TMath::Power(locProtonP4_Kin.Pz(),2));

		if(showOutput) { cout << "Filling photon quantity vectors!" << endl; }
		// Neutral Track. We define it this way so that we can easily loop through the 4 photons when filling.
                photonThetas[0] = locPhoton1P4_Kin.Theta()*radToDeg;
                photonThetas[1] = locPhoton2P4_Kin.Theta()*radToDeg;
                photonThetas[2] = locPhoton3P4_Kin.Theta()*radToDeg;
                photonThetas[3] = locPhoton4P4_Kin.Theta()*radToDeg;
                photonEnergies[0] = locPhoton1P4_Kin.E();
                photonEnergies[1] = locPhoton2P4_Kin.E();
                photonEnergies[2] = locPhoton3P4_Kin.E();
                photonEnergies[3] = locPhoton4P4_Kin.E();
                photonPhis[0] = locPhoton1P4_Kin.Phi()*radToDeg;
                photonPhis[1] = locPhoton2P4_Kin.Phi()*radToDeg;
                photonPhis[2] = locPhoton3P4_Kin.Phi()*radToDeg;
                photonPhis[3] = locPhoton4P4_Kin.Phi()*radToDeg;
                photonXs_Kin[0] = locPhoton1X4_Kin.X();
                photonXs_Kin[1] = locPhoton2X4_Kin.X();
                photonXs_Kin[2] = locPhoton3X4_Kin.X();
                photonXs_Kin[3] = locPhoton4X4_Kin.X();
                photonYs_Kin[0] = locPhoton1X4_Kin.Y();
                photonYs_Kin[1] = locPhoton2X4_Kin.Y();
                photonYs_Kin[2] = locPhoton3X4_Kin.Y();
                photonYs_Kin[3] = locPhoton4X4_Kin.Y();
                photonZs_Kin[0] = locPhoton1X4_Kin.Z();
                photonZs_Kin[1] = locPhoton2X4_Kin.Z();
                photonZs_Kin[2] = locPhoton3X4_Kin.Z();
                photonZs_Kin[3] = locPhoton4X4_Kin.Z();
                photonTs_Kin[0] = locPhoton1X4_Kin.T();
                photonTs_Kin[1] = locPhoton2X4_Kin.T();
                photonTs_Kin[2] = locPhoton3X4_Kin.T();
                photonTs_Kin[3] = locPhoton4X4_Kin.T();
                photonXs_Shower[0] = locPhoton1X4_Shower.X();
                photonXs_Shower[1] = locPhoton2X4_Shower.X();
                photonXs_Shower[2] = locPhoton3X4_Shower.X();
                photonXs_Shower[3] = locPhoton4X4_Shower.X();
                photonYs_Shower[0] = locPhoton1X4_Shower.Y();
                photonYs_Shower[1] = locPhoton2X4_Shower.Y();
                photonYs_Shower[2] = locPhoton3X4_Shower.Y();
                photonYs_Shower[3] = locPhoton4X4_Shower.Y();
                photonZs_Shower[0] = locPhoton1X4_Shower.Z();
                photonZs_Shower[1] = locPhoton2X4_Shower.Z();
                photonZs_Shower[2] = locPhoton3X4_Shower.Z();
                photonZs_Shower[3] = locPhoton4X4_Shower.Z();
                photonTs_Shower[0] = locPhoton1X4_Shower.T();
                photonTs_Shower[1] = locPhoton2X4_Shower.T();
                photonTs_Shower[2] = locPhoton3X4_Shower.T();
                photonTs_Shower[3] = locPhoton4X4_Shower.T();
                photonDeltaTs[0] = locPhoton1X4_Kin.T()-(dComboWrapper->Get_RFTime() + (locPhoton1X4_Kin.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
                photonDeltaTs[1] = locPhoton2X4_Kin.T()-(dComboWrapper->Get_RFTime() + (locPhoton2X4_Kin.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
                photonDeltaTs[2] = locPhoton3X4_Kin.T()-(dComboWrapper->Get_RFTime() + (locPhoton3X4_Kin.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
                photonDeltaTs[3] = locPhoton4X4_Kin.T()-(dComboWrapper->Get_RFTime() + (locPhoton4X4_Kin.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
                photonDetectedSyss[0] = dPhoton1Wrapper->Get_Detector_System_Timing();
                photonDetectedSyss[1] = dPhoton2Wrapper->Get_Detector_System_Timing();
                photonDetectedSyss[2] = dPhoton3Wrapper->Get_Detector_System_Timing();
                photonDetectedSyss[3] = dPhoton4Wrapper->Get_Detector_System_Timing();

		if(showOutput) { cout << "Getting some combo specific variables like CL" << endl; }
                //Kinematic fit variables
                locCLKinFit = dComboWrapper->Get_ConfidenceLevel_KinFit("");
                locUnusedEnergy = dComboWrapper->Get_Energy_UnusedShowers();
		locNumExtraNeutralShowers = Get_NumNeutralHypos()-4;
		locDOFKinFit = dComboWrapper->Get_NDF_KinFit("");
                locChiSqKinFit = dComboWrapper->Get_ChiSq_KinFit("");
                //TString DOF; DOF.Form("%d", dComboWrapper->Get_NDF_KinFit(""));
                //TString title = "Num DOF: "; title+=DOF.Data();

		// Shower shape variables
		if(showOutput) { cout << "Getting photon shower variables" << endl; }
                E1E9_FCAL[0] = dPhoton1Wrapper->Get_E1E9_FCAL();
                E1E9_FCAL[1] = dPhoton2Wrapper->Get_E1E9_FCAL();
                E1E9_FCAL[2] = dPhoton3Wrapper->Get_E1E9_FCAL();
                E1E9_FCAL[3] = dPhoton4Wrapper->Get_E1E9_FCAL();
                E9E25_FCAL[0] = dPhoton1Wrapper->Get_E9E25_FCAL();
                E9E25_FCAL[1] = dPhoton2Wrapper->Get_E9E25_FCAL();
                E9E25_FCAL[2] = dPhoton2Wrapper->Get_E9E25_FCAL();
                E9E25_FCAL[3] = dPhoton3Wrapper->Get_E9E25_FCAL();
                SumU_FCAL[0] = dPhoton1Wrapper->Get_SumU_FCAL();
                SumU_FCAL[1] = dPhoton2Wrapper->Get_SumU_FCAL();
                SumU_FCAL[2] = dPhoton3Wrapper->Get_SumU_FCAL();
                SumU_FCAL[3] = dPhoton4Wrapper->Get_SumU_FCAL();
                SumV_FCAL[0] = dPhoton1Wrapper->Get_SumV_FCAL();
                SumV_FCAL[1] = dPhoton2Wrapper->Get_SumV_FCAL();
                SumV_FCAL[2] = dPhoton3Wrapper->Get_SumV_FCAL();
                SumV_FCAL[3] = dPhoton4Wrapper->Get_SumV_FCAL();
                Energy_BCALPreshower[0] = dPhoton1Wrapper->Get_Energy_BCALPreshower();
                Energy_BCALPreshower[1] = dPhoton2Wrapper->Get_Energy_BCALPreshower();
                Energy_BCALPreshower[2] = dPhoton3Wrapper->Get_Energy_BCALPreshower();
                Energy_BCALPreshower[3] = dPhoton4Wrapper->Get_Energy_BCALPreshower();
                Energy_BCAL[0] = dPhoton1Wrapper->Get_Energy_BCAL();
                Energy_BCAL[1] = dPhoton2Wrapper->Get_Energy_BCAL();
                Energy_BCAL[2] = dPhoton3Wrapper->Get_Energy_BCAL();
                Energy_BCAL[3] = dPhoton4Wrapper->Get_Energy_BCAL();
                SigLong_BCAL[0] = dPhoton1Wrapper->Get_SigLong_BCAL();
                SigLong_BCAL[1] = dPhoton2Wrapper->Get_SigLong_BCAL();
                SigLong_BCAL[2] = dPhoton3Wrapper->Get_SigLong_BCAL();
                SigLong_BCAL[3] = dPhoton4Wrapper->Get_SigLong_BCAL();
                SigTheta_BCAL[0] = dPhoton1Wrapper->Get_SigTheta_BCAL();
                SigTheta_BCAL[1] = dPhoton2Wrapper->Get_SigTheta_BCAL();
                SigTheta_BCAL[2] = dPhoton3Wrapper->Get_SigTheta_BCAL();
                SigTheta_BCAL[3] = dPhoton4Wrapper->Get_SigTheta_BCAL();
		// We will grab the quality values first and then apply a cut to allow just the FCAL photons in.
		if(showOutput) { cout << "Getting photon shower quality" << endl; }
		showerQuality_FCAL[0] = dPhoton1Wrapper->Get_Shower_Quality();
		showerQuality_FCAL[1] = dPhoton2Wrapper->Get_Shower_Quality();
		showerQuality_FCAL[2] = dPhoton3Wrapper->Get_Shower_Quality();
		showerQuality_FCAL[3] = dPhoton4Wrapper->Get_Shower_Quality();
		pShowerQuality = pShowerQuality0*pShowerQuality1*pShowerQuality2*pShowerQuality3;
                SigTrans_BCAL[0] = dPhoton1Wrapper->Get_SigTrans_BCAL();
                SigTrans_BCAL[1] = dPhoton2Wrapper->Get_SigTrans_BCAL();
                SigTrans_BCAL[2] = dPhoton3Wrapper->Get_SigTrans_BCAL();
                SigTrans_BCAL[3] = dPhoton4Wrapper->Get_SigTrans_BCAL();
		DeltaPhi_BCAL[0] = dPhoton1Wrapper->Get_TrackBCAL_DeltaPhi();
		DeltaPhi_BCAL[1] = dPhoton1Wrapper->Get_TrackBCAL_DeltaPhi();
		DeltaPhi_BCAL[2] = dPhoton2Wrapper->Get_TrackBCAL_DeltaPhi();
		DeltaPhi_BCAL[3] = dPhoton3Wrapper->Get_TrackBCAL_DeltaPhi();
		DeltaZ_BCAL[0] = dPhoton1Wrapper->Get_TrackBCAL_DeltaZ();
		DeltaZ_BCAL[1] = dPhoton1Wrapper->Get_TrackBCAL_DeltaZ();
		DeltaZ_BCAL[2] = dPhoton2Wrapper->Get_TrackBCAL_DeltaZ();
		DeltaZ_BCAL[3] = dPhoton3Wrapper->Get_TrackBCAL_DeltaZ();
		DOCA_FCAL[0] = dPhoton1Wrapper->Get_TrackFCAL_DOCA();
		DOCA_FCAL[1] = dPhoton1Wrapper->Get_TrackFCAL_DOCA();
		DOCA_FCAL[2] = dPhoton2Wrapper->Get_TrackFCAL_DOCA();
		DOCA_FCAL[3] = dPhoton3Wrapper->Get_TrackFCAL_DOCA();


		if(showOutput) { cout << "Getting some shower variables" << endl; }
                locE1E9_FCAL_proton = dProtonWrapper->Get_E1E9_FCAL();
                locE9E25_FCAL_proton = dProtonWrapper->Get_E9E25_FCAL();
                locSumU_FCAL_proton = dProtonWrapper->Get_SumU_FCAL();
                locSumV_FCAL_proton = dProtonWrapper->Get_SumV_FCAL();
                locEnergy_BCALPreshower_proton = dProtonWrapper->Get_Energy_BCALPreshower();
                locEnergy_BCAL_proton = dProtonWrapper->Get_Energy_BCAL();
                locSigTheta_BCAL_proton = dProtonWrapper->Get_SigTheta_BCAL();
                locSigTrans_BCAL_proton = dProtonWrapper->Get_SigTrans_BCAL();;
                locSigLong_BCAL_proton = dProtonWrapper->Get_SigLong_BCAL();;





                // Timing variables, first section is about the proton
                // Get_RFTime is quoted at the center of the target. We have to shift it outwards.
                // We shift it  by the time it takes for light to travel between the initial target and the recoiled target.  
                RFtime = dComboWrapper->Get_RFTime() + (locProtonX4_Kin.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458;
                RFtimeProton = locProtonX4_Kin.T() - RFtime;

                // Might not make sense to shift relative to the BeamX4.Z since we cant really track the photon in the detector. Use Protons maybe
                //double locDeltaTRF = locBeamX4.T() - (dComboWrapper->Get_RFTime() + (locBeamX4.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
                locDeltaTRF = locBeamX4_Kin.T() - RFtime;

                // Accidental Subtracting we do a weighted histogram. Where we scale the increment to our histogram by the weight. The weight is 1 for in central beam bunch RF time and -0.5 for off central.
                weightAS = 1;
                if (14>abs(locDeltaTRF) && abs(locDeltaTRF) >2){ weightAS = -0.16667;} //about 1/6 since we capture a sample of 3 outside peaks
                // Shifted realtive to the proton
                inBCAL = 0; inTOF = 0; inFCAL = 0; inSTART = 0; inSYS_NULL = 0;
                if (dProtonWrapper->Get_Detector_System_Timing() == SYS_BCAL){ inBCAL = 1;}
                if (dProtonWrapper->Get_Detector_System_Timing() == SYS_TOF){ inTOF = 1;}
                if (dProtonWrapper->Get_Detector_System_Timing() == SYS_FCAL){ inFCAL = 1;}
                if (dProtonWrapper->Get_Detector_System_Timing() == SYS_START){ inSTART = 1;}
                if (dProtonWrapper->Get_Detector_System_Timing() == SYS_NULL){ inSYS_NULL = 1;}

		if(showOutput) { cout << "Got timing quantites" << endl; }


		// CALCULATING PAIRWISE QUANTITES BETWEEN PHOTONS.
		// the weird thing about this is that I have to do a cut on the dij3 of all the pairwise combinations. So I have to cacluate it before I can start filling most of the other histograms (that uses that dij3Cut)
		// I also have to worry about uniquess tracking when trying to do the filling here otherwise the dij3 noncut and cut graphs would share the same tracking which is wrong. The one with cuts should only fill the 
		// unique set when the cut has been passed. SO! We have to do our uniqueness tracking out here rather than before we fill. Since we will have done the uniquenss tracking here we do not have to do it in the future
		// when we are trying to fill the histograms.
		int countBothInFCAL = 0;
		int countBothInBCAL = 0;
		int countNotInEither = 0;
		int hist_id=0;
		double dij3;// for the FCAL
		double angle_ij;
		double deltaZ_ij;
		double deltaPhi_ij;
		double dij;// for the BCAL
		double R_BCALi;
		double z_BCALi;
		double theta_BCALi;
		double R_BCALj;
		double z_BCALj;
		double theta_BCALj;
		double delta_Rij_BCAL;
		double delta_Zij_BCAL;
		double delta_Thetaij_BCAL;
                std::vector<double> dij3VecFCAL[id+1];
                std::vector<double> dijVecBCAL[id+1];
                std::vector<double> angle_ijVec[id+1];
                std::vector<double> deltaZ_ijVec[id+1];
                std::vector<double> deltaPhi_ijVec[id+1];
		std::vector<std::vector<Int_t>> photonBCALPIDs_ijVec[id+1]; // used to match out dij3, angle_ij,... to the corresponding pid pairs. Does not persists through the event, combo specific.
		std::vector<std::vector<Int_t>> photonFCALPIDs_ijVec[id+1];
		std::vector<Int_t> photonBCALPIDs_ij;
		std::vector<Int_t> photonFCALPIDs_ij;
		map<Particle_t, set<Int_t>> usingPairIdsFCAL;// will be used to do our uniqueness trackingo over the entire event
		map<Particle_t, set<Int_t>> usingPairIdsBCAL;
		// This will just be used to help do the pN uniqueness tracking.	
                std::vector<Int_t> beingUsedNeutralIds; // should be reset after each section uniqueness tracking 
		if(showOutput){ cout << "\nCalculating pairwise quantities between photons, doing uniqueness tracking now instead of later since there will be a vector of data to fill\n-----------------------------------------" << endl;}
		std::vector<TLorentzVector> photonX4 = {locPhoton1X4_Shower, locPhoton2X4_Shower, locPhoton3X4_Shower, locPhoton4X4_Shower};
		std::vector<Int_t> photonIds = {locPhoton1NeutralID, locPhoton2NeutralID, locPhoton3NeutralID, locPhoton4NeutralID} ;
		if(showOutput){ cout << "Set up photonX4 and photonIds to use for pairwise quantity calculating!" << endl;}
		for (int i=0;i<4;++i){
			for(int j=i+1; j<4;j++){
				if(showOutput){ cout << "-- ith,jth photon: " << std::to_string(i) << ", " << std::to_string(j) << endl;}
				// if both the photons are not in the resppective detector there is no point on continuing to check the histograms + cuts + uniquenss tracking
				if(photonDetectedSyss[i]==SYS_FCAL && photonDetectedSyss[j]==SYS_FCAL){
					++countBothInFCAL;
					// usingPairIdsFCAL persists throughout the event since it it created outside it
					// photonFCALPIDs_ijVec persists only through the combo which is fine we will will just use it to fill.
					usingPairIdsFCAL.clear();
					usingPairIdsFCAL[Gamma].insert(photonIds[i]);
					usingPairIdsFCAL[Gamma].insert(photonIds[j]);
					photonFCALPIDs_ij.clear();
					photonFCALPIDs_ij.push_back(photonIds[i]);
					photonFCALPIDs_ij.push_back(photonIds[j]);
		                	dij3 = (photonX4[i]-photonX4[j]).Vect().Mag();
					if(showOutput){ cout << "BOTH IN FCAL - Got dij3! but is it unique?" << endl;}
					groupVec_it = 11;
					for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
						hist_id = vec_group_ids[groupVec_it][groupHist_it];
		                		if (usedPairIdsFCAL[hist_id].find(usingPairIdsFCAL)==usedPairIdsFCAL[hist_id].end()){
							if (histCuts[hist_id]){
		                				usedPairIdsFCAL[hist_id].insert(usingPairIdsFCAL);
								dij3VecFCAL[hist_id].push_back(dij3);
								photonFCALPIDs_ijVec[hist_id].push_back(photonFCALPIDs_ij);
								if(showOutput){cout << "Is Unique and Passed Cuts, Filling - Name of FCAL tracking hist: " << histList[hist_id][0] << endl;}
							}
							else if(showOutput){cout << "Not Passed - Name of FCAL tracking hist: " << histList[hist_id][0] << endl;}
						}
						if(showOutput){cout<<"Is not unique, checking next photon pair" << endl; }
					}
				}
				else if(photonDetectedSyss[i]==SYS_BCAL && photonDetectedSyss[j]==SYS_BCAL){
					++countBothInBCAL;
					usingPairIdsBCAL.clear();
					usingPairIdsBCAL[Gamma].insert(photonIds[i]);
					usingPairIdsBCAL[Gamma].insert(photonIds[j]);
					photonBCALPIDs_ij.clear();
					photonBCALPIDs_ij.push_back(photonIds[i]);
					photonBCALPIDs_ij.push_back(photonIds[j]);
		                	angle_ij = (photonX4[i].Vect().Angle(photonX4[j].Vect()))*radToDeg;
		                	if (angle_ij > 180) { angle_ij = 360-angle_ij; }
		                	deltaZ_ij = abs(photonX4[i].Z()-photonX4[j].Z());
		                	deltaPhi_ij = abs(photonX4[i].Phi()-photonX4[j].Phi())*radToDeg;
		                	if (deltaPhi_ij > 180) { deltaPhi_ij = 360-deltaPhi_ij; }
					if(showOutput){ cout << "BOTH IN BCAL - Got BCAL angle, deltaZ, deltaPhi, distanceOnCylinder!" << endl;}
		                	//dij3 = (photonX4[i]-photonX4[j]).Vect().Mag();

					// a better metric to use would be to use dij which is the distance conneting two points on the cylindrical BCAL BUT constraining to be on the surface to tell us about separation
					// 65 < R < 90 cm is the inner and outer dimensions of the BCAL
					R_BCALi = TMath::Sqrt(TMath::Sq(photonX4[i].X())+TMath::Sq(photonX4[i].Y()));
					z_BCALi = photonX4[i].Z();
					theta_BCALi = TMath::ATan(photonX4[i].Y()/photonX4[i].X());
					R_BCALj = TMath::Sqrt(TMath::Sq(photonX4[j].X())+TMath::Sq(photonX4[j].Y()));
					z_BCALj = photonX4[j].Z();
					theta_BCALj = TMath::ATan(photonX4[j].Y()/photonX4[j].X());
					delta_Rij_BCAL = R_BCALj - R_BCALi;
					delta_Zij_BCAL = z_BCALj - z_BCALi;
					delta_Thetaij_BCAL = theta_BCALj - theta_BCALi;
					dij = TMath::Sqrt(TMath::Sq(delta_Rij_BCAL)+TMath::Sq((65+delta_Rij_BCAL)*delta_Thetaij_BCAL)+TMath::Sq(delta_Zij_BCAL));
					groupVec_it = 13;
					for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
						hist_id = vec_group_ids[groupVec_it][groupHist_it];
		                		if (usedPairIdsBCAL[hist_id].find(usingPairIdsBCAL)==usedPairIdsBCAL[hist_id].end()){
							if (histCuts[hist_id]){
								dijVecBCAL[hist_id].push_back(dij);
								angle_ijVec[hist_id].push_back(angle_ij);
								deltaZ_ijVec[hist_id].push_back(deltaZ_ij);
								deltaPhi_ijVec[hist_id].push_back(deltaPhi_ij);
		                        			usedPairIdsBCAL[hist_id].insert(usingPairIdsBCAL);
								photonBCALPIDs_ijVec[hist_id].push_back(photonBCALPIDs_ij);
								if(showOutput){ cout << "Is Unique and Passed Cuts, Filling - Name of BCAL tracking hist: " << histList[hist_id][0] << endl;}
							}
							else if(showOutput){cout << "Not Passed - Name of BCAL tracking hist: " << histList[hist_id][0] << endl;}
						}
						if(showOutput){cout<<"Is not unique, checking next photon pair" << endl; }
					}
				}
				else if((photonDetectedSyss[i]==SYS_BCAL && photonDetectedSyss[j]==SYS_FCAL) || (photonDetectedSyss[i]==SYS_FCAL && photonDetectedSyss[j]==SYS_BCAL)){
					++countNotInEither;
					if(showOutput){ cout << "ONE IN FCAL ONE IN BCAL" << endl;}
				}
				else { cout << "\n\n*************************\nERROR - Pairwise both in FCAL, both in FCAL, in either, is not complete!\n***********************************" << endl;}
			}
		}
		if((countBothInBCAL+countNotInEither+countBothInFCAL)==6){ if(showOutput) {cout << "\n\n*********************\nPairwise photons counted correctly!\n*******************\nNum both in FCAL: " + std::to_string(countBothInFCAL) +"\nNum both in BCAL: " + std::to_string(countBothInBCAL) + "\nNum one in either: " + std::to_string(countNotInEither)<< endl; }}
		else { cout << "\n\n*********************\nERROR - Pairwise photons not counted correctly!\n*******************\n" << endl;}



		if(showOutput){ cout << endl << endl; }
                if(showOutput){ cout << "Calculating prodplanephi and cosTheta in different system" << endl;}
                // calculating phi
                //double decayPlanePhi = dAnalysisUtilities.Calc_DecayPlanePsi_Vector_3BodyDecay(locBeamE, Proton, locProtonP4_Kin, mixingPi0Eta_Kin, locPi0P4_Kin, locEtaP4_Kin, locDecayPlaneTheta);
                double prodPlanePhi = dAnalysisUtilities.Calc_ProdPlanePhi_Pseudoscalar(locBeamE, Proton, mixingPi0Eta_Kin);
                double polarizationAngle;
                // we kind of have to do this since we use sed to search and replace locDegAngle/45/90/135...
                string deg = degAngle.substr(degAngle.length()-3);
                if (deg == "000") { polarizationAngle = 0; }
                else if (deg == "045") { polarizationAngle = 45; }
                else if (deg == "090") { polarizationAngle = 90; }
                // there is uniqueness in the second digit. By the construction of using the last two characters in the string we have to ignmore 1 in 135 but it is unambiguous.
                else if (deg == "135") { polarizationAngle = 135; }
                // if it isnt one of the above then degAngle must be deg000
                else { polarizationAngle = -999; } // just wanted a number negative enough to make locPhi, no matter the value, always negative so there is some indication to look for.
                locPhi = polarizationAngle - prodPlanePhi;

                // Calculating kinematic variables like t and cosTheta
                mandelstam_t = (locProtonP4_Kin-dTargetP4).M2();
                mandelstam_t_pe = (locBeamP4_Kin-mixingPi0Eta).M2();

		// We will also calculate the cos theta in the lab frame. This is used to see if we are actually getting the right amount of events in the detector
		theta_pi0_lab = locPi0P4.Theta()*radToDeg;
		theta_eta_lab = locEtaP4.Theta()*radToDeg;

		// Determine the CM vectors to be used to calculate the angles in the Hel frame
		// also calculating the cosTheta of pi0eta system, pi0, eta  in the CM frame
		TLorentzVector cm_vec = locBeamP4_Kin+dTargetP4;
		TLorentzVector mixingPi0Eta_cm = mixingPi0Eta_Kin;
		TLorentzVector pi0_cm = locPi0P4_Kin;
		TLorentzVector eta_cm = locEtaP4_Kin;
                TLorentzVector beam_cm = locBeamP4_Kin;
                TLorentzVector recoil_cm = locProtonP4_Kin;
		TLorentzVector largestEPhotonInPi0_CM;
		if (locPhoton1P4_Kin.E() > locPhoton2P4_Kin.E()){
			largestEPhotonInPi0_CM = locPhoton1P4_Kin;
		}
		else {
			largestEPhotonInPi0_CM = locPhoton2P4_Kin;
		}
		largestEPhotonInPi0_CM.Boost(-cm_vec.BoostVector());
		cosTheta_largestEinPi0_CM = largestEPhotonInPi0_CM.CosTheta();
	
		mixingPi0Eta_cm.Boost(-cm_vec.BoostVector());
		pi0_cm.Boost(-cm_vec.BoostVector());
		eta_cm.Boost(-cm_vec.BoostVector());
		beam_cm.Boost(-cm_vec.BoostVector());
		recoil_cm.Boost(-cm_vec.BoostVector());

		cosTheta_pi0eta_CM = mixingPi0Eta_cm.CosTheta();
		cosTheta_pi0_CM = pi0_cm.CosTheta();
		cosTheta_eta_CM = eta_cm.CosTheta();
		theta_pi0_CM = pi0_cm.Theta();
		theta_eta_CM = eta_cm.Theta();
		phi_pi0eta_CM = mixingPi0Eta_cm.Phi()*radToDeg;
		phi_pi0_CM = pi0_cm.Phi()*radToDeg;
		phi_eta_CM = eta_cm.Phi()*radToDeg;

		//mandelstam_t0 = TMath::Power((locProtonP4_Kin.M2()-mixingPi0Eta_Kin.M2()-dTargetP4.M2())/(2*(locBeamP4_Kin+dTargetP4).M()),2)-(beam_cm-mixingPi0Eta_cm).M2();
		mandelstam_t0 = TMath::Power((locBeamP4_Kin.M2()-mixingPi0Eta_Kin.M2()-dTargetP4.M2()+locProtonP4_Kin.M2())/(2*(locBeamP4_Kin+dTargetP4).M()),2)-(beam_cm-mixingPi0Eta_cm).M2();
		mandelstam_tp = abs(mandelstam_t-mandelstam_t0);
		mandelstam_tp_pe = abs(mandelstam_t_pe-mandelstam_t0);

		// Calculate cosTheta, phi in maybe the helicity axes.
                TVector3 z = mixingPi0Eta_cm.Vect().Unit();
		// this y should be the normal of the production plane. If we do a boost in a direction in the production plane the perp direction doesn't change. We could use the beam and the recoiled proton to define the
		// production plane in this new frame. 
                TVector3 y = beam_cm.Vect().Cross(z).Unit();
                TVector3 x = y.Cross(z).Unit();

		TVector3 pi0_cm_unit = pi0_cm.Vect().Unit();	
		TVector3 eta_cm_unit = eta_cm.Vect().Unit();	
		TVector3 pi0eta_cm_unit = mixingPi0Eta_cm.Vect().Unit();	

                TVector3 norm = eta_cm.Vect().Cross(pi0_cm.Vect()).Unit();
                TVector3 angles(   norm.Dot(x),
                                   norm.Dot(y),
                                   norm.Dot(z) );
		// this tvector3 angles is not really the angles. Each component would be the cos(angle) where angle is the angle needed to rotate the vector into one of the axes. This new unit vector is not defined in the new axes
		TVector3 angles_pi0 ( pi0_cm_unit.Dot(x), pi0_cm_unit.Dot(y), pi0_cm_unit.Dot(z) );
		TVector3 angles_eta ( eta_cm_unit.Dot(x), eta_cm_unit.Dot(y), eta_cm_unit.Dot(z) );
		TVector3 angles_pi0eta ( pi0eta_cm_unit.Dot(x), pi0eta_cm_unit.Dot(y), pi0eta_cm_unit.Dot(z) );

                cosTheta_decayPlane_hel = angles.CosTheta();
                phi_decayPlane_hel = angles.Phi()*radToDeg;
		cosTheta_pi0_hel = angles_pi0.CosTheta();
		cosTheta_eta_hel = angles_eta.CosTheta();
		cosTheta_pi0eta_hel = angles_pi0eta.CosTheta();
		theta_pi0_hel = angles_pi0.Theta();
		theta_eta_hel = angles_eta.Theta();
		phi_pi0_hel = angles_pi0.Phi()*radToDeg;
		phi_eta_hel = angles_eta.Phi()*radToDeg;
		phi_pi0eta_hel = angles_pi0eta.Phi()*radToDeg;

		
		// defining the vectors to be used for the GJ frame
                //TLorentzVector beam_res = locBeamP4_Kin;
                //TLorentzVector recoil_res = locProtonP4_Kin;
                //TLorentzVector eta_res = locEtaP4_Kin;
                //TLorentzVector pi0_res = locPi0P4_Kin;
		//TLorentzVector mixingPi0Eta_res = mixingPi0Eta_Kin;

                //beam_res.Boost(-mixingPi0Eta_Kin.BoostVector());
                //recoil_res.Boost(-mixingPi0Eta_Kin.BoostVector());
                //eta_res.Boost(-mixingPi0Eta_Kin.BoostVector());
                //pi0_res.Boost(-mixingPi0Eta_Kin.BoostVector());
		//mixingPi0Eta_res.Boost(-mixingPi0Eta_Kin.BoostVector());

                TLorentzVector beam_res = beam_cm;
                TLorentzVector recoil_res = recoil_cm;
                TLorentzVector eta_res = eta_cm;
                TLorentzVector pi0_res = pi0_cm;
		TLorentzVector mixingPi0Eta_res = mixingPi0Eta_cm;

		// Boost using mixingPi0Eta_cm instead of res... since that is what we are changing.
                beam_res.Boost(-mixingPi0Eta_cm.BoostVector());
                recoil_res.Boost(-mixingPi0Eta_cm.BoostVector());
                eta_res.Boost(-mixingPi0Eta_cm.BoostVector());
                pi0_res.Boost(-mixingPi0Eta_cm.BoostVector());
		mixingPi0Eta_res.Boost(-mixingPi0Eta_cm.BoostVector());
		TLorentzVector largestEPhotonInEta_GJ;
		if (locPhoton3P4_Kin.E() > locPhoton4P4_Kin.E()){
			largestEPhotonInEta_GJ = locPhoton3P4_Kin;
		}
		else {
			largestEPhotonInEta_GJ = locPhoton4P4_Kin;
		}
		largestEPhotonInEta_GJ.Boost(-cm_vec.BoostVector());
		largestEPhotonInEta_GJ.Boost(-mixingPi0Eta_cm.BoostVector());
		TVector3 largestEPhotonInEta_res_unit = largestEPhotonInEta_GJ.Vect().Unit();

		TVector3 pi0_res_unit = pi0_res.Vect().Unit();	
		TVector3 eta_res_unit = eta_res.Vect().Unit();	
		TVector3 pi0eta_res_unit = mixingPi0Eta_res.Vect().Unit();	
		// Calculate cosTheta, phi in maybe the GJ axes.
		// since we already defined the x,y,z as TVector3 we don't have to do it again.
                z = beam_res.Vect().Unit();
		// this y should be the normal of the production plane. If we do a boost in a direction in the production plane the perp direction doesn't change. We could use the beam and the recoiled proton to define the
		// production plane in this new frame. Let us define it in the CM frame. 
                y = mixingPi0Eta_cm.Vect().Cross(beam_cm.Vect()).Unit();
		locYDotZ_GJ=y.Dot(z);
                x = y.Cross(z).Unit();

		// recylcing our use of norm and the defined angles
                norm = eta_res.Vect().Cross(pi0_res.Vect()).Unit();
                angles.SetXYZ(   norm.Dot(x),
                          norm.Dot(y),
                          norm.Dot(z) );
		angles_pi0.SetXYZ ( pi0_res_unit.Dot(x), pi0_res_unit.Dot(y), pi0_res_unit.Dot(z) );
		angles_eta.SetXYZ ( eta_res_unit.Dot(x), eta_res_unit.Dot(y), eta_res_unit.Dot(z) );
		angles_pi0eta.SetXYZ ( pi0eta_res_unit.Dot(x), pi0eta_res_unit.Dot(y), pi0eta_res_unit.Dot(z) );
		TVector3 angles_largestEinEta_GJ = largestEPhotonInEta_res_unit; 
		angles_largestEinEta_GJ.SetXYZ( angles_largestEinEta_GJ.Dot(x), angles_largestEinEta_GJ.Dot(y), angles_largestEinEta_GJ.Dot(z) );

                cosTheta_decayPlane_GJ = angles.CosTheta();
                phi_decayPlane_GJ = angles.Phi()*radToDeg;
		theta_pi0_GJ = angles_pi0.Theta()*radToDeg;
		theta_eta_GJ = angles_eta.Theta()*radToDeg;
		cosTheta_pi0_GJ = angles_pi0.CosTheta();
		cosTheta_eta_GJ = angles_eta.CosTheta();
		cosTheta_largestEinEta_GJ = angles_largestEinEta_GJ.CosTheta();
		//cout << "*cosTheta*" << endl;
		//cout << "Boosting first into the CM frame then boosting into the resonance frame with the resonances momentum vector" << endl;
		//cout << "The resonance's momentum vector in the boosted frame without axis redefinition:" << endl;
		//cout << mixingPi0Eta_res.Vect().X() << ","<< mixingPi0Eta_res.Vect().Y() << ","<< mixingPi0Eta_res.Vect().Z() << endl;
		//cout << "The resonance's momentum UNIT vector in the boosted frame without axis redefinition:" << endl;
		//cout << mixingPi0Eta_res.Vect().Unit().X() << ","<< mixingPi0Eta_res.Vect().Unit().Y() << ","<< mixingPi0Eta_res.Vect().Unit().Z() << endl;
		//cout << "Projecting into the the GJ angles we get the following Vector:" << endl;
		//cout << angles_pi0eta.X() << ","<< angles_pi0eta.Y() << ","<< angles_pi0eta.Z() << endl;
		//cout << angles_pi0eta.CosTheta() << endl;
		//cout << mixingPi0Eta_res.Vect().Dot(z)/mixingPi0Eta_res.Vect().Mag() << endl;
		//cout << "*cosTheta*" << endl;

		phi_pi0_GJ = angles_pi0.Phi()*radToDeg;
		phi_eta_GJ = angles_eta.Phi()*radToDeg;	

		// Angle() returns the angle between two vectors where the angle is in radians.
		angleBetweenPi0Eta = pi0_res_unit.Angle(eta_res_unit)*radToDeg;
		if(showOutput) { cout << "Filling inCone bools" << endl; }
		for (int iDelta=0; iDelta<2; ++iDelta){
			// THere should be some redundancy here. since mixingPi0Eta is the sum fo the individual 4 vectors it must be that the pi0 and eta decaay back to back in the resonance rest frame.
			// Therefore largeAngle=180 always and only 1 of the mesons have to be in the cone. As a check it is good to see if the number of times it passes these 3 degenerate conditions are the same!
			pi0_inCone[iDelta] = abs(theta_pi0_GJ)>90-delta[iDelta] && abs(theta_pi0_GJ)<90+delta[iDelta] && abs(phi_pi0_GJ)>90-delta[iDelta] && abs(phi_pi0_GJ)<90+delta[iDelta]; 
			eta_inCone[iDelta] = abs(theta_eta_GJ)>90-delta[iDelta] && abs(theta_eta_GJ)<90+delta[iDelta] && abs(phi_eta_GJ)>90-delta[iDelta] && abs(phi_eta_GJ)<90+delta[iDelta]; 
			largeAngle[iDelta] = angleBetweenPi0Eta>180-2*delta[iDelta]; // that is the closest they can be
			withinCone[iDelta] = largeAngle[iDelta] && pi0_inCone[iDelta] && eta_inCone[iDelta]; 
		}
		q = sqrt(pi0_cm.Pz()*pi0_cm.Pz() + eta_cm.Pz()*eta_cm.Pz() + recoil_cm.Pz()*recoil_cm.Pz());
		pi0_cmZ = pi0_cm.Pz();
		eta_cmZ = eta_cm.Pz();
		recoil_cmZ = recoil_cm.Pz();
		omega = 0;
   		if(recoil_cmZ > 0 && pi0_cmZ < 0 && eta_cmZ > 0) {
   		        // in quadrant I, so asin returns correct answer
   		        omega = asin(sqrt(3/2.)*recoil_cmZ/q);
   		}
   		if(recoil_cmZ > 0 && pi0_cmZ < 0 && eta_cmZ < 0) {
   		        // use pi0_cmZ as X-axis since gives simpler quadrant definition
   		        omega = asin(sqrt(3/2.)*pi0_cmZ/q) + 2*TMath::Pi()/3;
   		}
   		if(recoil_cmZ > 0 && pi0_cmZ > 0 && eta_cmZ < 0) {
   		        // in quadrant II, so need to subtract asin from PI to get correct angle
   		        omega = TMath::Pi() - asin(sqrt(3/2.)*recoil_cmZ/q);
   		}
   		if(recoil_cmZ < 0 && pi0_cmZ > 0 && eta_cmZ < 0) {
   		        // in quadrant III, so need to subtract asin from PI to get correct angle
   		        omega = TMath::Pi() - asin(sqrt(3/2.)*recoil_cmZ/q);
   		}
   		if(recoil_cmZ < 0 && pi0_cmZ > 0 && eta_cmZ > 0) {
   		        // use pi0_cmZ as X-axis since gives simpler quadrant definition
   		        omega = (TMath::Pi() - asin(sqrt(3/2.)*pi0_cmZ/q)) + 2*TMath::Pi()/3;
   		}
   		if(recoil_cmZ < 0 && pi0_cmZ < 0 && eta_cmZ > 0) {
   		        // in quadrant IV, so asin returns the correct answer
   		        omega = asin(sqrt(3/2.)*recoil_cmZ/q);
   		}

		// compute vanHove_x and y before changing omega from rad to degrees
		vanHove_x = q*cos(omega);
		vanHove_y = q*sin(omega);
		
		omega = omega*radToDeg;
		if(showOutput) { cout << "Got Vanhove quantities" << endl; }

		// eow that we gotten the PIDs we can also do our uniqueness tracking setup now.
		if(showOutput){cout << "Setting up the sets and maps to be used in uniqueness tracking" << endl;}

// ********************************************************************************************************************************
// 			// UNIQUENESS TRACKING 
// ********************************************************************************************************************************
		//Uniqueness tracking:
		map<Particle_t, set<Int_t> > usingCombo;
		usingCombo[Unknown].insert(locBeamID); //beam
		usingCombo[Proton].insert(locProtonTrackID);
		usingCombo[Gamma].insert(locPhoton1NeutralID);
		usingCombo[Gamma].insert(locPhoton2NeutralID);
		usingCombo[Gamma].insert(locPhoton3NeutralID);
		usingCombo[Gamma].insert(locPhoton4NeutralID);

		map<Particle_t, set<Int_t> > usingPh12;
		usingPh12[Gamma].insert(locPhoton1NeutralID);
		usingPh12[Gamma].insert(locPhoton2NeutralID);

		map<Particle_t, set<Int_t> > usingPh34;
		usingPh34[Gamma].insert(locPhoton3NeutralID);
		usingPh34[Gamma].insert(locPhoton4NeutralID);

		map<Particle_t, set<Int_t> > usingPh13;
		usingPh13[Gamma].insert(locPhoton1NeutralID);
		usingPh13[Gamma].insert(locPhoton3NeutralID);

		map<Particle_t, set<Int_t> > usingPh24;
		usingPh24[Gamma].insert(locPhoton2NeutralID);
		usingPh24[Gamma].insert(locPhoton4NeutralID);

		map<Particle_t, set<Int_t> > usingPh1234;
		usingPh1234[Gamma].insert(locPhoton1NeutralID);
		usingPh1234[Gamma].insert(locPhoton2NeutralID);
		usingPh1234[Gamma].insert(locPhoton3NeutralID);
		usingPh1234[Gamma].insert(locPhoton4NeutralID);

		map<Particle_t, set<Int_t> > usingPh12P1;
		usingPh12P1[Gamma].insert(locPhoton1NeutralID);
		usingPh12P1[Gamma].insert(locPhoton2NeutralID);
		usingPh12P1[Gamma].insert(locProtonTrackID);

		map<Particle_t, set<Int_t> > usingPh34P1;
		usingPh34P1[Proton].insert(locProtonTrackID);
		usingPh34P1[Gamma].insert(locPhoton3NeutralID);
		usingPh34P1[Gamma].insert(locPhoton4NeutralID);

		map<Particle_t, set<Int_t> > usingPh12B1;
		usingPh34P1[Unknown].insert(locProtonTrackID);
		usingPh34P1[Gamma].insert(locPhoton1NeutralID);
		usingPh34P1[Gamma].insert(locPhoton2NeutralID);

		map<Particle_t, set<Int_t> > usingPh34B1;
		usingPh34P1[Unknown].insert(locProtonTrackID);
		usingPh34P1[Gamma].insert(locPhoton3NeutralID);
		usingPh34P1[Gamma].insert(locPhoton4NeutralID);

		map<Particle_t, set<Int_t> > usingPh1234P1;
		usingPh1234P1[Proton].insert(locProtonTrackID);
		usingPh1234P1[Gamma].insert(locPhoton1NeutralID);
		usingPh1234P1[Gamma].insert(locPhoton2NeutralID);
		usingPh1234P1[Gamma].insert(locPhoton3NeutralID);
		usingPh1234P1[Gamma].insert(locPhoton4NeutralID);

		
// ********************************************************************************************************************************
		// Recomputing some booleans since they could not have been sucessfully initialized without running though the loop first.
// ********************************************************************************************************************************
                if(showOutput){cout << "Re-evaluating bool variables to fill into histCuts" << endl;}
		
		// Beam cuts
		pBeamAsymE = locBeamE > 8.0 && locBeamE < 8.7;
		pBeamE30to46 = locBeamE > 3.0 && locBeamE < 4.6;
		pBeamE46to62 = locBeamE > 4.6 && locBeamE < 6.2;
		pBeamE62to78 = locBeamE > 6.2 && locBeamE < 7.8;
		pBeamE78to94 = locBeamE > 7.8 && locBeamE < 9.4;
		pBeamE94to11 = locBeamE > 9.4 && locBeamE < 11.0;
		pBeamE8GeVPlus = locBeamE > 8.0;
		
		dEdxCut = TMath::Power(10,-6)*(0.9+TMath::Exp(3.0-3.5*(locMagP3Proton+0.05)/.93827)); // The bigger then number multiplying MagP3 the sharper the cut.
		// pi0Eta specifc cuts
		pEtaProtonBaryonCut = locEtaProton_Kin >= etaProtonBaryonCut;
		ppi0ProtonBaryonCut = locPi0Proton_Kin >= pi0ProtonBaryonCut;
		pBeamE = locBeamE  >= beamECut;

		iLowMass = lowMass; // reinitialize
		iUpMass = lowMass+binScale; // reinitialize
		for (int bin=0; bin<numBinsMass; ++bin){
		        p_phiMassBinned[bin] = locPi0Eta_Kin>iLowMass && locPi0Eta_Kin<iUpMass; 
		        iLowMass+=binScale;
		        iUpMass+=binScale;
		}
		iLowMass_t = lowMass_t; // reinitialize
		iUpMass_t = lowMass_t+binScale_t; // reinitialize
		for (int bin=0; bin<numBinsMass_t; ++bin){
			// we use this set of bools to bin the t distribution in bins of M(pi0eta)
		        p_tMassBinned[bin] = locPi0Eta_Kin>iLowMass_t && locPi0Eta_Kin<iUpMass_t; 
		        iLowMass_t+=binScale_t;
		        iUpMass_t+=binScale_t;
		}
		iUpUE = 0.1; // such that after running through 20 bins we will reach a UE cut of 1GeV.
		for (int bin=0; bin<numRegions_UE; ++bin){
			p_pi0MassEtaMassUEregion[bin] = locUnusedEnergy<iUpUE;
                	iUpUE+=0.1;
		}
		iUpChiSq = 5; 
		for (int bin=0; bin<numRegions_ChiSq; ++bin){
			p_pi0MassEtaMassChiSqregion[bin] = locChiSqKinFit<iUpChiSq;
                	iUpChiSq+=5;
		}
		// THIS IS FOR THE pi0Mass binned in E(pi0) selecting out the f2 peak
		for (int bin=0; bin<numRegions_E; ++bin){
            		p_pi0MassPi0Eregion_1[bin] = locPi0E_Kin<iUpE[bin] && locPi0E_Kin>iLowE[bin];
            		p_pi0MassPi0Eregion_2[bin] = locEtaE_Kin<iUpE[bin] && locEtaE_Kin>iLowE[bin];
		}
		pSelectf2 = locPi0Eta_Kin<1.5 && locPi0Eta_Kin>1;

		outsideEllipse_loose=0;
                if (locPi0Mass_Kin<=(ellipseX-ellipseXr_loose) || locPi0Mass_Kin>=(ellipseX+ellipseXr_loose)){ outsideEllipse = 1;}
                double ellipseDeltaY = TMath::Sqrt((1 - TMath::Power((locPi0Mass_Kin-ellipseX)/ellipseXr_loose,2))*TMath::Power(ellipseYr_loose,2));
                if (locEtaMass_Kin>=(ellipseDeltaY+ellipseY) || locEtaMass_Kin<=(-ellipseDeltaY+ellipseY)){ outsideEllipse = 1; }
		pinsideEllipse_loose = !outsideEllipse_loose;

		outsideEllipse = 0;
                if (locPi0Mass_Kin<=(ellipseX-ellipseXr) || locPi0Mass_Kin>=(ellipseX+ellipseXr)){ outsideEllipse = 1;}
                ellipseDeltaY = TMath::Sqrt((1 - TMath::Power((locPi0Mass_Kin-ellipseX)/ellipseXr,2))*TMath::Power(ellipseYr,2));
                if (locEtaMass_Kin>=(ellipseDeltaY+ellipseY) || locEtaMass_Kin<=(-ellipseDeltaY+ellipseY)){ outsideEllipse = 1; }
		pinsideEllipse = !outsideEllipse;

		outsideEllipseBS1 = 0;
                if (locPi0Mass_Kin<=(ellipseXBS1-ellipseXrBS1) || locPi0Mass_Kin>=(ellipseXBS1+ellipseXrBS1)){ outsideEllipseBS1 = 1;}
                double ellipseDeltaYBS1 = TMath::Sqrt((1 - TMath::Power((locPi0Mass_Kin-ellipseXBS1)/ellipseXrBS1,2))*TMath::Power(ellipseYrBS1,2));
                if (locEtaMass_Kin>=(ellipseDeltaYBS1+ellipseYBS1) || locEtaMass_Kin<=(-ellipseDeltaYBS1+ellipseYBS1)){ outsideEllipseBS1 = 1; }
		pinsideEllipseBS1 = !outsideEllipseBS1;

		outsideEllipseBS2 = 0;
                if (locPi0Mass_Kin<=(ellipseXBS2-ellipseXrBS2) || locPi0Mass_Kin>=(ellipseXBS2+ellipseXrBS2)){ outsideEllipseBS2 = 1;}
                double ellipseDeltaYBS2 = TMath::Sqrt((1 - TMath::Power((locPi0Mass_Kin-ellipseXBS2)/ellipseXrBS2,2))*TMath::Power(ellipseYrBS2,2));
                if (locEtaMass_Kin>=(ellipseDeltaYBS2+ellipseYBS2) || locEtaMass_Kin<=(-ellipseDeltaYBS2+ellipseYBS2)){ outsideEllipseBS2 = 1; }
		pinsideEllipseBS2 = !outsideEllipseBS2;

		pYellowBKG=pinsideEllipseBS2*outsideEllipseBS1;

		// The sharp transition in the weights would only be reflected in the hists using kin values!
		if(pYellowBKG){ 
			// this histogram is just to check if the region I am selecting is good
			dHist_checkEllipseBS[0]->Fill(locPi0Mass_Kin, locEtaMass_Kin);
			weightBS=-areaRatio; 
		} 
		else if (pinsideEllipse){
			dHist_checkEllipseBS[1]->Fill(locPi0Mass_Kin, locEtaMass_Kin);
			weightBS=1;
		}
		else { 
			weightBS=0;
			dHist_checkEllipseBS[2]->Fill(locPi0Mass_Kin, locEtaMass_Kin);
		}
		
		// now that we have defined both the weights we can multiply them together
		//weight = weightAS*weightBS;
		weight = weightAS;

		// General Cuts
		pUnusedEnergy = locUnusedEnergy <= unusedEnergyCut;
		//pCLKinFit1 = locCLKinFit >= CLCut1;
		//pCLKinFit = locChiSqKinFit <= ChiSqCut;	
		pChiSq = locChiSqKinFit <= ChiSqCut;	
		chiSq100 = locChiSqKinFit <= 100;
		//pCLKinFit3 = locCLKinFit >= CLCut3;
		//pCLKinFit4 = locCLKinFit >= CLCut4;
		//pCLKinFit5 = locCLKinFit >= CLCut5;
		//pCLKinFit6 = locCLKinFit >= CLCut6;
		pDeltaTRF = abs(locDeltaTRF) <= RFCut;
		pMissingMassSquared = locMissingMassSquared <= MMsqCut;
		
		// Neutral Cuts
                pdij3pass = 1;
		// we must index the dij3Vec which is a vector with size id. id_noCutDij3 would be the index of the no cut dij3 histogram
		// which we would use to calculate whether the dij3 cut passes. Probably the most fair way to calculate it rather than
		// using the one where we have already uniquenss trakced it.
                for (std::size_t i=0; i<dij3VecFCAL[id_noCutDij3].size();i++){
			// since we want all the dij3FCAL to be <= 12.5 cm we can check if any of them is > 12.5 since its simpler. 
                        if ( dij3VecFCAL[id_noCutDij3][i] <= dijCut){
                                pdij3pass = 0;
                                break;
                        }
                }
		pPhoton1E = photonEnergies[0] >= ECut;
		pPhoton2E = photonEnergies[1] >= ECut;
		pPhoton3E = photonEnergies[2] >= ECut;
		pPhoton4E = photonEnergies[3] >= ECut;
		pPhotonE = pPhoton1E*pPhoton2E*pPhoton3E*pPhoton4E;
		pPhoton1Theta = (photonThetas[0]*radToDeg>=thetaCutMin && photonThetas[0]*radToDeg<=thetaCutMax1) || photonThetas[0]*radToDeg>=thetaCutMax2;
		pPhoton2Theta = (photonThetas[1]*radToDeg>=thetaCutMin && photonThetas[1]*radToDeg<=thetaCutMax1) || photonThetas[1]*radToDeg>=thetaCutMax2;
		pPhoton3Theta = (photonThetas[2]*radToDeg>=thetaCutMin && photonThetas[2]*radToDeg<=thetaCutMax1) || photonThetas[2]*radToDeg>=thetaCutMax2;
		pPhoton4Theta = (photonThetas[3]*radToDeg>=thetaCutMin && photonThetas[3]*radToDeg<=thetaCutMax1) || photonThetas[3]*radToDeg>=thetaCutMax2;
		pPhotonTheta = pPhoton1Theta*pPhoton2Theta*pPhoton3Theta*pPhoton4Theta;
		
		pShowerQuality0 = showerQuality_FCAL[0] > 0.5;
		pShowerQuality1 = showerQuality_FCAL[1] > 0.5;
		pShowerQuality2 = showerQuality_FCAL[2] > 0.5;
		pShowerQuality3 = showerQuality_FCAL[3] > 0.5;
		// Charged Cuts
		pMagP3Proton = locMagP3Proton >= P3Cut;
		pzCutmin = zCutmin <= locdzProton && locdzProton <= zCutmax;
		pRProton = locRProton <= Rcut;
		pdEdxCDCProton = locdEdxCDCProton >= dEdxCut;
		pReg1 = locPzProton>Reg1Xmin && locPtProton>Reg1Ymin && locPzProton<Reg1Xmax && locPtProton<Reg1Ymax;
		pReg2 = locPzProton>Reg2Xmin && locPtProton>Reg2Ymin && locPzProton<Reg2Xmax && locPtProton<Reg2Ymax;
		pReg3 = locPzProton>Reg3Xmin && locPtProton>Reg3Ymin && locPzProton<Reg3Xmax && locPtProton<Reg3Ymax;
		pReg4 = locPzProton>Reg4Xmin && locPtProton>Reg4Ymin && locPzProton<Reg4Xmax && locPtProton<Reg4Ymax;
		
		// locWherePhoton will be set equal to photonDetectedSyss[N] 
		for (int i = 0; i<4; ++i){
			pPhotonInBCALorFCAL[i] = photonDetectedSyss[i] == SYS_BCAL || photonDetectedSyss[i] == SYS_FCAL;
			pPhotonInFCAL[i] = photonDetectedSyss[i] == SYS_FCAL;
			pPhotonInBCAL[i] = photonDetectedSyss[i] == SYS_BCAL;
		}

		pPi0InFCAL = photonDetectedSyss[0]==SYS_FCAL && photonDetectedSyss[1]==SYS_FCAL;
		pPi0InBCAL = photonDetectedSyss[0]==SYS_BCAL && photonDetectedSyss[1]==SYS_BCAL;
		pPi0InSplit = (photonDetectedSyss[0]==SYS_FCAL && photonDetectedSyss[1]==SYS_BCAL) || (photonDetectedSyss[1]==SYS_FCAL && photonDetectedSyss[0]==SYS_BCAL);
		pEtaInFCAL = photonDetectedSyss[2]==SYS_FCAL && photonDetectedSyss[3]==SYS_FCAL;
		pEtaInBCAL = photonDetectedSyss[2]==SYS_BCAL && photonDetectedSyss[3]==SYS_BCAL;
		pEtaInSplit = (photonDetectedSyss[2]==SYS_FCAL && photonDetectedSyss[3]==SYS_BCAL) || (photonDetectedSyss[3]==SYS_FCAL && photonDetectedSyss[2]==SYS_BCAL);
		
		pPi0InFCAL_mismatch = photonDetectedSyss[0]==SYS_FCAL && photonDetectedSyss[2]==SYS_FCAL;
		pPi0InBCAL_mismatch = photonDetectedSyss[0]==SYS_BCAL && photonDetectedSyss[2]==SYS_BCAL;
		pPi0InSplit_mismatch = (photonDetectedSyss[0]==SYS_FCAL && photonDetectedSyss[2]==SYS_BCAL) || (photonDetectedSyss[2]==SYS_FCAL && photonDetectedSyss[0]==SYS_BCAL);
		pEtaInFCAL_mismatch = photonDetectedSyss[1]==SYS_FCAL && photonDetectedSyss[3]==SYS_FCAL;
		pEtaInBCAL_mismatch = photonDetectedSyss[1]==SYS_BCAL && photonDetectedSyss[3]==SYS_BCAL;
		pEtaInSplit_mismatch = (photonDetectedSyss[1]==SYS_FCAL && photonDetectedSyss[3]==SYS_BCAL) || (photonDetectedSyss[3]==SYS_FCAL && photonDetectedSyss[1]==SYS_BCAL);

		if (omega < 0) {
			pVanHove = omega > -120 && omega < -60; 
		}
		else{
			pVanHove = omega > 240 && omega < 300; 
		}

		ptLT1 = mandelstam_tp<1; 
		// cut to accept the delta peak in M(pi0proton)
		// we will actually use it to reject the delta peak in the cuts. pi0pi0 should also have this resonance right?
		pMPi0P14 = locPi0Proton_Kin<1.4;

		// Various combinations of cuts, the majority of them will be used just for a few histograms like when showing unused energy graph we will use mUE which
		// removes the UE cut from allGeneralCutsPassed. m prefix basically stands for minus

		pShowerQuality=1;
		allGeneralCutsPassed = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mMPi0P14 = pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mBeamE = !pMPi0P14*pShowerQuality*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mMMSq = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pdEdxCDCProton*pinsideEllipse;
		//pDiffCL = pBeamE8GeVPlus*pUnusedEnergy*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pdEdxCDCProton*pinsideEllipse*pMissingMassSquared;
		mRProton = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mRProtonZMin = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mdEdxCDC = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pinsideEllipse;
		mZMin = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mMagP3 = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mPhotonE = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mPhotonTheta = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mdij3 = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mUE = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mUEChiSq = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
		mChiSq = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pdEdxCDCProton*pinsideEllipse*pMissingMassSquared;
		// ------ 
		mEllipse_pre = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
		mEllipseUE_pre = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pChiSq*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
		mEllipseUEChiSq_pre = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
		mEllipseChiSq_pre = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pDeltaTRF*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
		// ------ Rejects the 0-weight region since that would mess with the uniqueness tracking. We use the OR logic to sum the subsets of the yellow and the red region.
		mEllipse = mEllipse_pre*pYellowBKG || allGeneralCutsPassed;
		mEllipseUE = mEllipseUE_pre*pYellowBKG || mEllipseUE_pre*pinsideEllipse;
		mEllipseUEChiSq = mEllipseUEChiSq_pre*pYellowBKG || mEllipseUEChiSq_pre*pinsideEllipse;
		mEllipseChiSq = mEllipseChiSq_pre*pYellowBKG || mEllipseChiSq_pre*pinsideEllipse;

		if (allGeneralCutsPassed*withinCone[1]) {
			if(showOutput){ cout << "$$$Checking Angles of the combo in GJ!!!\n" << angles_pi0.X() << ","<< angles_pi0.Y() << ","<< angles_pi0.Z() << "," << angles_eta.X() << ","<< angles_eta.Y() << ","<< angles_eta.Z() << endl;}
		}
		if(pShowerQuality){ ++count_ShowerQuality; dHist_Cuts->Fill(cutNames[0],1);}
		if(pBeamE8GeVPlus){ ++count_BeamE8GeVPlus; dHist_Cuts->Fill(cutNames[1],1);} 
                if(pUnusedEnergy){ ++count_UnusedEnergy;dHist_Cuts->Fill(cutNames[2],1);}
                if(pChiSq){ ++count_ChiSq;dHist_Cuts->Fill(cutNames[3],1);}
                if(pDeltaTRF){ ++count_DeltaTRF;dHist_Cuts->Fill(cutNames[4],1);}
                if(pdij3pass){ ++count_dij3pass;dHist_Cuts->Fill(cutNames[5],1);}
                if(pPhotonE){ ++count_PhotonE;dHist_Cuts->Fill(cutNames[6],1);}
                if(pPhotonTheta){ ++count_PhotonTheta;dHist_Cuts->Fill(cutNames[7],1);}
                if(pMagP3Proton){ ++count_MagP3Proton;dHist_Cuts->Fill(cutNames[8],1);}
                if(pzCutmin){ ++count_zCutmin;dHist_Cuts->Fill(cutNames[9],1);}
                if(pRProton){ ++count_RProton;dHist_Cuts->Fill(cutNames[10],1);}
                if(pMissingMassSquared){ ++count_MissingMassSquared;dHist_Cuts->Fill(cutNames[11],1);}
                if(pdEdxCDCProton){ ++count_dEdxCDCProton;dHist_Cuts->Fill(cutNames[12],1);}
                if(pinsideEllipse){ ++count_insideEllipse;dHist_Cuts->Fill(cutNames[13],1);}
		if(allGeneralCutsPassed){ ++count_allGeneralCutsPassed;dHist_Cuts->Fill(cutNames[14],1);}


		if(showOutput) {cout << "Start Filling histVals and histCuts" << endl;}

		
		// Have to re-evaluate the booleans
		id = -1; // have to reinitialize id
		for (int i = 0; i < 2; ++i){
		        if (i<1) { cutString = ""; cutsToApply = 1; applyAccSub=noAccSub;}
		        else { cutString = "Cut"; cutsToApply = allGeneralCutsPassed; applyAccSub=weight;}
                	++id; histVals[id] = {8, applyAccSub, locMissingMassSquared};
			//histList[id] = {("MissingMassSquared"+cutString).c_str(), ";Missing Mass Squared (GeV/c^{2})^{2}", "200", "-0.2", "0.2", "Events / 0.002 GeV/c^{2}"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mMMSq; }
			++id; histVals[id] = {1, applyAccSub, locBeamE};
                	//histList[id] = {("BeamEnergy"+cutString).c_str(), ";Beam Energy (GeV)", "120", "0", "12", "Events / 0.1 GeV" };
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mBeamE; }
	       		++id; histVals[id] = {8, applyAccSub, locYDotZ_GJ};
	       		//histList[id] = {("yDotZ_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";Cos(theta) angle between Z and Y in GJ").c_str(), "100", "0", "1", "Events / 0.01"};
	       		histCuts[id] = cutsToApply;
			// Invariant Mass Hists
			++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
			//histList[id] = {("pi0proton1D"+cutString).c_str(), ";M(#pi_{0}proton) (GeV)", "400", "0", "4", "Events / 0.01 GeV"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin};
			//histList[id] = {("etaproton1D"+cutString).c_str(), ";M(#etaproton) (GeV)", "450", "0", "4.5", "Events / 0.01 GeV"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locPi0Proton_Kin};
			//histList[id] = {("pi0etaPi0Proton"+cutString).c_str(), ";", "160", "0.", "4", "90", "0.", "4.5", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#pi_{0}Proton) (GeV) with Events / 0.05 GeV"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locEtaProton_Kin};
			//histList[id] = {("pi0etaEtaProton"+cutString).c_str(), ";", "160", "0.", "4", "100", "0.", "5", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#etaProton) (GeV) with Events / 0.05 GeV"};
			histCuts[id] = cutsToApply;
            		if (!is_pi0eta){
	    		   	++id; histVals[id] = {5, applyAccSub, locPi0E_Kin};
	    		   	//histList[id] = {("pi0E_1"+cutString).c_str(), ("Cuts="+cutsApplied+";E (GeV)").c_str(), "100", "0", "4", "Events / 0.04 GeV"};
	    		   	histCuts[id] = cutsToApply*pSelectf2;
	    		   	++id; histVals[id] = {6, applyAccSub, locEtaE_Kin};
	    		   	//histList[id] = {("pi0E_2"+cutString).c_str(), ("Cuts="+cutsApplied+";E (GeV)").c_str(), "100", "0", "4", "Events / 0.04 GeV"};
	    		   	histCuts[id] = cutsToApply*pSelectf2;
            		}
		
			// Kinematic Hists
			++id; histVals[id] = {8, applyAccSub, locPhi}; 
			//histList[id] = {("prodPlanePSphi"+cutString).c_str(), ";#phi", "180", "-540", "540", "Events / 6 degrees"};
			histCuts[id] = cutsToApply;
			//++id; histVals[id] = {8, applyAccSub, cosTheta_decayPlane_hel};
			////histList[id] = {("decayPlane_cosTheta_hel"+cutString).c_str(), "cos(#theta) of #pi_{0}+#eta", "100","-1","1", "Events / 0.02"};
			//histCuts[id] = cutsToApply;
			//++id; histVals[id] = {8, applyAccSub, phi_decayPlane_hel};
			////histList[id] = {("decayPlane_phi_hel"+cutString).c_str(), "#phi of #pi_{0}+#eta", "160","-1.6","1.6", "Events / 0.02"};
			//histCuts[id] = cutsToApply;
	       		++id; histVals[id] = {8, applyAccSub, omega};
	       		//histList[id] = {("vanHove_omega"+cutString).c_str(), ("Cuts="+cutsApplied+";Van Hove Omega Plot").c_str(), "120","-360","360", "60", "Events / 6 degrees"};
	       		histCuts[id] = cutsToApply;
			if(outputGJ_HEL){
	       			++id; histVals[id] = {8, applyAccSub, cosTheta_pi0_hel};
	       			//histList[id] = {("pi0_cosTheta_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#eta").c_str(), "100","-1","1", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, cosTheta_eta_hel};
	       			//histList[id] = {("eta_cosTheta_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#eta").c_str(), "100","-1","1", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, cosTheta_pi0eta_hel};
	       			//histList[id] = {("pi0eta_cosTheta_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#eta").c_str(), "100","-1","1", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, phi_eta_hel};
	       			//histList[id] = {("eta_phi_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}+#eta").c_str(), "160","-1.6","1.6", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, phi_pi0_hel};
	       			//histList[id] = {("pi0_phi_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}+#eta").c_str(), "160","-1.6","1.6", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, phi_pi0eta_hel};
	       			//histList[id] = {("pi0eta_phi_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}+#eta").c_str(), "160","-1.6","1.6", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
				
	       			++id; histVals[id] = {8, applyAccSub, theta_pi0_GJ};
	       			//histList[id] = {("pi0_theta_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#eta").c_str(), "100","-1","1", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, theta_eta_GJ};
	       			//histList[id] = {("eta_theta_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#eta").c_str(), "100","-1","1", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, cosTheta_eta_GJ};
	       			//histList[id] = {("eta_cosTheta_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#eta").c_str(), "100","-1","1", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, cosTheta_pi0_GJ};
	       			//histList[id] = {("pi0_cosTheta_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of #pi_{0}+#eta").c_str(), "100","-1","1", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_eta_GJ};
	       			//histList[id] = {("eta_cosTheta_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+"cosTheta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "Events / 0.01 GeV", "Events / 0.2"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_pi0_GJ};
	       			//histList[id] = {("pi0_cosTheta_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+"cosTheta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "Events / 0.01 GeV", "Events / 0.2"};
	       			histCuts[id] = cutsToApply;
                        	++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, phi_pi0_GJ};
                        	//histList[id] = {("pi0_phi_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+" phi of #pi_{0} vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "90","-180","180", "M(#pi_{0}#eta) Events / 0.01 GeV", "Events / 4 degrees"};
                        	histCuts[id] = cutsToApply;
                        	++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, phi_eta_GJ};
                        	//histList[id] = {("eta_phi_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+" phi of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "90","-180","180", "M(#pi_{0}#eta)Events / 0.01 GeV", "Events / 4 degrees"};
                        	histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, phi_eta_GJ};
	       			//histList[id] = {("eta_phi_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}+#eta").c_str(), "160","-1.6","1.6", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;
	       			++id; histVals[id] = {8, applyAccSub, phi_pi0_GJ};
	       			//histList[id] = {("pi0_phi_GJ"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of #pi_{0}+#eta").c_str(), "160","-1.6","1.6", "Events / 0.02"};
	       			histCuts[id] = cutsToApply;

				++id; histVals[id] = {8, applyAccSub, cosTheta_pi0eta_CM};
				//histList[id] = {("pi0eta_cosTheta_CM"+cutString).c_str(), "cos(#theta) of #pi_{0}+#eta", "100","-1","1", "Events / 0.02"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {3, applyAccSub, cosTheta_pi0_CM};
				//histList[id] = {("pi0_cosTheta_CM"+cutString).c_str(), "cos(#theta) of #pi_{0}", "100","-1","1", "Events / 0.02"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {7, applyAccSub, cosTheta_eta_CM};
				//histList[id] = {("eta_cosTheta_CM"+cutString).c_str(), "cos(#theta) of #eta", "100","-1","1", "Events / 0.02"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {8, applyAccSub, phi_pi0eta_CM};
				//histList[id] = {("pi0eta_phi_CM"+cutString).c_str(), "#phi of #pi_{0}+#eta", "160","-1.6","1.6", "Events / 0.02"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {3, applyAccSub, phi_pi0_CM};
				//histList[id] = {("pi0_phi_CM"+cutString).c_str(), "#phi of #pi_{0}", "160", "-1.6", "1.6", "Events / 0.02"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {7, applyAccSub, phi_eta_CM};
				//histList[id] = {("eta_phi_CM"+cutString).c_str(), "#phi of #eta", "160", "-1.6", "1.6", "Events / 0.02"};
				histCuts[id] = cutsToApply;
			}

			// Invariant Mass Hists EBeam > 8 GeV
			//++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
			////histList[id] = {("pi0proton1D8GeVPlus"+cutString).c_str(), ";M(#pi_{0}proton) (GeV)", "400", "0", "4", "Events / 0.01 GeV"};
			//histCuts[id] = cutsToApply*pBeamE8GeVPlus;
			//++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin}; 
			////histList[id] = {("etaproton1D8GeVPlus"+cutString).c_str(), ";M(#etaproton) (GeV)", "450", "0", "4.5", "Events / 0.01 GeV"};
			//histCuts[id] = cutsToApply*pBeamE8GeVPlus;
			//++id; histVals[id] = {11, applyAccSub, locEtaProton_Kin,locPi0Eta_Kin}; 
			////histList[id] = {("pi0etaProton8GeVPlus"+cutString).c_str(), ";", "90", "0.", "4.5", "160", "0.", "4", "M(#etaProton) (GeV) with Events / 0.05 GeV", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV"};
			//histCuts[id] = cutsToApply*pBeamE8GeVPlus;
			//++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locPi0Proton_Kin};
			////histList[id] = {("pi0etaPi0Proton8GeVPlus"+cutString).c_str(), ";", "160", "0.", "4", "90", "0.", "4.5", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#pi_{0}Proton) (GeV) with Events / 0.05 GeV"};
			//histCuts[id] = cutsToApply*pBeamE8GeVPlus;
			//++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locEtaProton_Kin};
			////histList[id] = {("pi0etaEtaProton8GeVPlus"+cutString).c_str(), ";", "160", "0.", "4", "100", "0.", "5", "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#etaProton) (GeV) with Events / 0.05 GeV"};
			//histCuts[id] = cutsToApply*pBeamE8GeVPlus;
		
			// Charged Track Hists
			++id; histVals[id] = {2, applyAccSub, locRProton}; 
			//histList[id] = {("RadiusProton"+cutString).c_str(), ";Radius(proton) (cm)", "200", "0" , "10", "Events / 0.05 cm"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mRProton; }
			++id; histVals[id] = {2, applyAccSub, locdzProton};
			//histList[id] = {("dzProton"+cutString).c_str(), ";z(Proton) (cm)", "160", "0" , "160", "Events / 1 cm"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mZMin; }
			++id; histVals[id] = {2, applyAccSub, locdEdxCDCProton};
			//histList[id] = {("dEdxProtonCDC"+cutString).c_str(), ";dEdx(proton) GeV/cm", "200", "0." , "0.00003", "Events / 1.5E-7 GeV/cm"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mdEdxCDC; }
			++id; histVals[id] = {2, applyAccSub, locPzProton}; 
			//histList[id] = {("PzProton"+cutString).c_str(), ";Pz GeV", "200", "0" , "5", "Events / 0.025 GeV"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {2, applyAccSub, locXProton,locYProton}; 
			//histList[id] = {("XYplaneProton"+cutString).c_str(), ";", "100", "-4", "4", "100", "-4", "4", "x(proton) (cm) with Events / 0.08 cm", "y(proton) (cm) with Events / 0.08 cm"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mRProton; }
			++id; histVals[id] = {2, applyAccSub, locRProton,locdzProton};
			//histList[id] = {("RZplaneProton"+cutString).c_str(), ";", "100", "0", "4", "100", "0", "100", "R(proton) (cm) wtih Events / 0.04 cm","z(proton) (cm) with Events / 1 cm"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mRProtonZMin; }
			++id; histVals[id] = {2, applyAccSub, locMagP3Proton,locdEdxCDCProton};
			//histList[id] = {("P3dEdxCDCProton"+cutString).c_str(), ";", "100", "0", "4", "100", "0", "0.00003", "Momentum(proton) CDC (GeV/c) with Events / 0.04 GeV/c", "dEdx(proton) CDC (GeV/cm) with Events / 3E-7"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = pMagP3Proton*pzCutmin*pRProton; } // mdEdxCDC was used before but now we do things sequentially....
			++id; histVals[id] = {2, applyAccSub, locMagP3Proton,locdEdxFDCProton};
			//histList[id] = {("P3dEdxFDCProton"+cutString).c_str(), ";", "100", "0", "4", "100", "0", "0.00003", "Momentum(proton) FDC (GeV/c) with Events / 0.04 GeV/c", "dEdx(proton) FDC (GeV/cm) with Events / 3E-7"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {2, applyAccSub, locPzProton,locPtProton}; 
			//histList[id] = {("PzPtProton"+cutString).c_str(), ";", "125", "0", "5", "125", "0", "2.5", "Pz(proton) (GeV/c) with Events / 0.04 GeV/c", "Pt(proton) (GeV/c) with Events / 0.02 GeV/c"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {2, applyAccSub, locPzProton,locPolarAngleProton};
			//histList[id] = {("PzThetaProton"+cutString).c_str(), ";", "125", "0", "5", "250", "0", "100", "Pz(proton) (GeV/c) with Events / 0.04 GeV/c", "Theta(proton) (degrees) with Events / 0.4 degrees"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {2, applyAccSub, locPtProton,locPolarAngleProton};
			//histList[id] = {("PtThetaProton"+cutString).c_str(), ";", "125", "0", "5", "250", "0", "100", "Pz(proton) (GeV/c) with Events / 0.04 GeV/c", "Theta(proton) (degrees) with Events / 0.4 degrees"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
			//histList[id] = {("PolarAngleProton"+cutString).c_str(), "#theta(Proton) degrees"," 200", "0", "100", "Events / 0.5 degrees"};
			histCuts[id] = cutsToApply;
			if (outputThetaRegion){
				++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
				//histList[id] = {("ThetaRegion1Proton"+cutString).c_str(), ";#theta(degrees)", "125", "0", "100", "Events / 0.08 degrees"};
				histCuts[id] = cutsToApply*pReg1; 
				++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
				//histList[id] = {("ThetaRegion2Proton"+cutString).c_str(), ";#theta(degrees)", "125", "0", "100", "Events / 0.08 degrees"};
				histCuts[id] = cutsToApply*pReg2; 
				++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
				//histList[id] = {("ThetaRegion3Proton"+cutString).c_str(), ";#theta(degrees)", "125", "0", "100", "Events / 0.08 degrees"};
				histCuts[id] = cutsToApply*pReg3; 
				++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
				//histList[id] = {("ThetaRegion4Proton"+cutString).c_str(), ";#theta(degrees)", "125", "0", "100", "Events / 0.08 degrees"};
				histCuts[id] = cutsToApply*pReg4; 
			}
		
		
			// Neutral Track Hists
			++id; histVals[id] = {13, applyAccSub, photonEnergies[0]};
			histVecVals[id] = {photonEnergies};
			//histList[id] = {("PhotonShowerE1"+cutString).c_str(), ";Energy(#gamma) GeV", "200", "0" , "10", "Events / 0.05 GeV"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mPhotonE; }
			++id; histVals[id] = {13, applyAccSub, photonThetas[0]}; 
			histVecVals[id] = {photonThetas};
			//histList[id] = {("PhotonShowerTheta1"+cutString).c_str(), ";PolarAngle(#gamma) degrees", "300", "0" , "150", "Events / 0.5 degrees"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mPhotonTheta; }
			++id; histVals[id] = {13, applyAccSub, photonEnergies[0], photonThetas[0]};
			histVecVals[id] = {photonEnergies, photonThetas};
			//histList[id] = {("thetaEPhotonShower"+cutString).c_str(), ";", "200", "0", "10", "400", "0", "150", "E(#gamma) (GeV) wtih Events / 0.05", "#theta(#gamma) (radians) with Events / 0.25"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {13, applyAccSub, photonXs_Shower[0], photonYs_Shower[0]}; 
			histVecVals[id] = {photonXs_Shower, photonYs_Shower};
			//histList[id] = {("XYPhotonShower"+cutString).c_str(), ";", "260", "-130", "130", "260", "-130", "130", "x(cm) with Events / 1 cm", "y(cm) with Events / 1 cm"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {13, applyAccSub, photonXs_Shower[0], photonYs_Shower[0]}; 
			histVecVals[id] = {photonXs_Shower, photonYs_Shower};
			//histList[id] = {("XYPhotonShowerBCAL"+cutString).c_str(), ";", "260", "-130", "130", "260", "-130", "130", "x(cm) with Events / 1 cm", "y(cm) with Events / 1 cm"};
			histCuts[id] = cutsToApply;
			medBool[0] = cutsToApply*pPhotonInBCAL[0];
			medBool[1] = cutsToApply*pPhotonInBCAL[1];
			medBool[2] = cutsToApply*pPhotonInBCAL[2];
			medBool[3] = cutsToApply*pPhotonInBCAL[3];
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
			++id; histVals[id] = {13, applyAccSub, photonPhis[0], photonThetas[0]};
			histVecVals[id] = {photonPhis, photonThetas};
			//histList[id] = {("ThetaPhiPhotonShower"+cutString).c_str(), ";", "360", "-180", "180", "130", "0", "130", "#phi(degrees) with Events / 1 cm", "#theta(degrees) with Events / 1 cm"};
		        if (cutString=="") { histCuts[id] = cutsToApply; }
		        else { histCuts[id] = mPhotonTheta; }
			// These below are tracked differently than the charged tracks and neutral tracks so they have their own unique identifiers.
			++id; histVals[id] = {12, applyAccSub, locPhotonDijFCAL}; 
			histVecVals[id] = {dij3VecFCAL[id]};
			//histList[id] = {("3DistanceBetweenPhotonsFCAL"+cutString).c_str(), ";FCAL - Distance between Pairs of Photons (cm)", "500", "0", "250", "Events / 0.5 cm"};
			if (cutString=="") { histCuts[id] = cutsToApply; }
			else { histCuts[id] = mdij3; } 
			++id; histVals[id] = {14, applyAccSub, locPhotonDijBCAL}; 
			histVecVals[id] = {dijVecBCAL[id]};
			//histList[id] = {("3DistanceBetweenPhotonsBCAL"+cutString).c_str(), ";BCAL - Distance between Pairs of Photons (cm)", "500", "0", "250", "Events / 0.5 cm"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {14, applyAccSub, locPhotonAij}; 
			histVecVals[id] = {angle_ijVec[id]};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {14, applyAccSub, locPhotonZij, locPhotonPij}; 
			histVecVals[id] = {deltaZ_ijVec[id], deltaPhi_ijVec[id]};
			//histList[id] = {("deltaZvsPhiPhotonBCAL"+cutString).c_str(), ";", "200", "0", "400", "100", "0", "200", "#DeltaZ  (cm) with Events / 2 cm",  "#Delta#phi between showers in BCAL (degrees) with Events / 2 degrees"};
			histCuts[id] = cutsToApply;
	


// *******************************************************************************************************************************************
// *******************************************************************************************************************************************
// ********************************************** ===== TIMING HISTOGRAMS === **************************************************************
// *******************************************************************************************************************************************

			if(outputTimingHists){
				++id; histVals[id] = {8, applyAccSub, locPolarAngleProton,RFtimeProton}; 
				//histList[id] = {("timeProtonRFvsTheta"+cutString).c_str(),"", "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees","T_{Proton} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
				//histList[id] = {("timeProtonRFvsP3"+cutString).c_str(),"", "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c", "T_{Proton} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply;
				// ********************** protonInBCAL ********************
				++id; histVals[id] = {8, applyAccSub, locPolarAngleProton,RFtimeProton}; 
				//histList[id] = {("timeBCALRFvsTheta"+cutString).c_str(), "", "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees","T_{BCAL} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				medBool[0]=cutsToApply*inBCAL;
				histVecCuts[id] = {medBool[0]};
				++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
				//histList[id] = {("timeBCALRFvsP3"+cutString).c_str(), "", "200", "0", "5", "120","-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{BCAL} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				histVecCuts[id] = {medBool[0]};
				// ********************** protonInTOF ********************
				++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
				//histList[id] = {("timeTOFRFvsTheta"+cutString).c_str(), "", "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees", "T_{TOF} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				medBool[0]=cutsToApply*inTOF;
				histVecCuts[id] = {medBool[0]};
				++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
				//histList[id] = {("timeTOFRFvsP3"+cutString).c_str(), "", "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{TOF} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				histVecCuts[id] = {medBool[0]};
				// ********************** protonInFCAL ********************
				++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
				//histList[id] = {("timeFCALRFvsTheta"+cutString).c_str(), "", "200"," 0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 cm", "T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				medBool[0]=cutsToApply*inFCAL;
				histVecCuts[id] = {medBool[0]};
				++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
				//histList[id] = {("timeFCALRFvsP3"+cutString).c_str(), "", "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c", "T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				histVecCuts[id] = {medBool[0]};
				// ********************** protonInSTART ********************
				++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
				//histList[id] = {("timeSTARTRFvsTheta"+cutString).c_str(), "", "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees", "T_{START} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				medBool[0]=cutsToApply*inSTART;
				histVecCuts[id] = {medBool[0]};
				++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
				//histList[id] = {("timeSTARTRFvsP3"+cutString).c_str(), "", "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{START} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				histVecCuts[id] = {medBool[0]};
				// ********************** protonInSYS_NULL ********************
				++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
				//histList[id] = {("timeSYSNULLRFvsTheta"+cutString).c_str(), "", "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees","T_{SYSNULL} - T_{RF} (ns) with Events / 0.1 ns "};
				histCuts[id] = cutsToApply; 
				medBool[0]=cutsToApply*inSYS_NULL;
				histVecCuts[id] = {medBool[0]};
				++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
				//histList[id] = {("timeSYSNULLRFvsP3"+cutString).c_str(), "", "200", "0", "5", "120", "-6" , "6","Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{SYSNULL} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				histVecCuts[id] = {medBool[0]};


				// ********************** photonInBCALorFCAL ********************
				++id; histVals[id] = {8, applyAccSub, photonThetas[0], photonDeltaTs[0]}; 
				histVecVals[id] = {photonThetas, photonDeltaTs};
				//histList[id] = {("timePhotonBCALFCALRFvsTheta"+cutString).c_str(),"","300","0","120","120","-6","6", "#theta(Proton) degrees with Events / 0.4 degrees","T_{BCAL/FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				medBool[0] = cutsToApply*pPhotonInBCALorFCAL[0];
				medBool[1] = cutsToApply*pPhotonInBCALorFCAL[1];
				medBool[2] = cutsToApply*pPhotonInBCALorFCAL[2];
				medBool[3] = cutsToApply*pPhotonInBCALorFCAL[3];
				histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
				++id; histVals[id] = {8, applyAccSub, photonDeltaTs[0]}; 
				histVecVals[id] = {photonDeltaTs};
				//histList[id] = {("timePhotonBCALFCALRF"+cutString).c_str(),"T_{BCAL/FCAL} - T_{RF} (ns)","120","-6","6", "Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};

				// ********************** photonInFCAL ********************
				++id; histVals[id] = {8, applyAccSub, photonThetas[0], photonDeltaTs[0]}; 
				histVecVals[id] = {photonThetas, photonDeltaTs};
				//histList[id] = {("timePhotonFCALRFvsTheta"+cutString).c_str(),"","60","0","12","120","-6","6", "#theta(Proton) degrees with Events / 0.2 degrees","T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				medBool[0] = cutsToApply*pPhotonInFCAL[0];
				medBool[1] = cutsToApply*pPhotonInFCAL[1];
				medBool[2] = cutsToApply*pPhotonInFCAL[2];
				medBool[3] = cutsToApply*pPhotonInFCAL[3];
				histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
				++id; histVals[id] = {8, applyAccSub, photonEnergies[0], photonDeltaTs[0]};
				histVecVals[id] = {photonEnergies, photonDeltaTs};
				//histList[id] = {("timePhotonFCALRFvsEnergy"+cutString).c_str(),"","200","0","12","120","-6","6", "Photon Energy (GeV) with Events / 0.06 GeV/c", "T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
				++id; histVals[id] = {8, applyAccSub, photonDeltaTs[0]};                                                                                    
				histVecVals[id] = {photonDeltaTs};
				//histList[id] = {("timePhotonFCALRF"+cutString).c_str(),"T_{FCAL} - T_{RF} (ns)","120","-6","6", "Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};

				// ********************** photonInBCAL ********************
				++id; histVals[id] = {8, applyAccSub, photonEnergies[0], photonDeltaTs[0]};
				histVecVals[id] = {photonEnergies, photonDeltaTs};
				//histList[id] = {("timePhotonBCALRFvsEnergy"+cutString).c_str(),"","200","0","6","120","-6","6", "Photon Energy (GeV) with Events / 0.03 GeV/c","T_{BCAL} - T_{RF} (ns) with Events / 0.1 ns"};
				histCuts[id] = cutsToApply; 
				medBool[0] = cutsToApply*pPhotonInBCAL[0];
				medBool[1] = cutsToApply*pPhotonInBCAL[1];
				medBool[2] = cutsToApply*pPhotonInBCAL[2];
				medBool[3] = cutsToApply*pPhotonInBCAL[3];
				histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
				++id; histVals[id] = {8, applyAccSub, photonDeltaTs[0]}; 
				histVecVals[id] = {photonDeltaTs};
				//histList[id] = {("timePhotonBCALRF"+cutString).c_str(),"T_{BCAL} - T_{RF} (ns)","120","-6","6", "Events / 0.1 ns"};
				histCuts[id] = cutsToApply;                                                                                                  
				histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
			}


// =========================================================== END OF TIMING HISTOGRAMS ===================================================
// =========================================================== END OF TIMING HISTOGRAMS ===================================================

		
			// Non track related hists
			++id; histVals[id] = {8, applyAccSub, locCLKinFit}; 
			//histList[id] = {("P4CLKinFit"+cutString).c_str(), ";Confidence Level" , "200", "0", ".30", "Entries / 0.0015"};
			if (cutString=="") { histCuts[id] = cutsToApply; }
			else { histCuts[id] = mChiSq; } 
			++id; histVals[id] = {8, applyAccSub, locUnusedEnergy};
			//histList[id] = {("UnusedEnergy"+cutString).c_str(), ";Unused Energy (GeV)", "450", "-1", "8", "Entries / 0.02 GeV"};
			if (cutString=="") { histCuts[id] = cutsToApply; }
			else { histCuts[id] = mUE; } 
			++id; histVals[id] = {8, applyAccSub, locChiSqKinFit}; 
	       		if (cutString=="") { histCuts[id] = cutsToApply; }
	       		else { histCuts[id] = mChiSq;} 
			//histList[id] = {("P4ChiSqKinFit"+cutString).c_str(), "Chi Squared;", "150", "0", "300", "Entries / 1"};
	       		++id; histVals[id] = {8, applyAccSub, locDOFKinFit}; 
	       		//histList[id] = {("P4DOFKinFit"+cutString).c_str(), ("Cuts="+cutsApplied+";DOF").c_str(), "30", "0", "30", "Entries / 1"};
	       		if (cutString=="") { histCuts[id] = cutsToApply; }
	       		else { histCuts[id] = mChiSq;} 

			//Shower shape variables
			// neutral showers
			medBool[0] = cutsToApply*pPhotonInFCAL[0];
			medBool[1] = cutsToApply*pPhotonInFCAL[1];
			medBool[2] = cutsToApply*pPhotonInFCAL[2];
			medBool[3] = cutsToApply*pPhotonInFCAL[3];
		        ++id; histVals[id] = {13, applyAccSub, E1E9_FCAL[0]}; 
		        //histList[id] = {("E1E9_FCAL"+cutString).c_str(), "E1E9_FCAL;", "100", "-5", "5", "Entries / 0.1"};
		        histCuts[id] = cutsToApply;
			histVecVals[id] = {E1E9_FCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
		        ++id; histVals[id] = {13, applyAccSub, E9E25_FCAL[0]}; 
		        //histList[id] = {("E9E25_FCAL"+cutString).c_str(), "E9E25_FCAL;", "100", "-5", "5", "Entries / 0.1"};
		        histCuts[id] = cutsToApply;
			histVecVals[id] = {E9E25_FCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
		        ++id; histVals[id] = {13, applyAccSub, SumU_FCAL[0]}; 
		        //histList[id] = {("SumU_FCAL"+cutString).c_str(), "SumU_FCAL;", "100", "-5", "5", "Entries / 0.1"};
		        histCuts[id] = cutsToApply;
			histVecVals[id] = {SumU_FCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
		        ++id; histVals[id] = {13, applyAccSub, SumV_FCAL[0]}; 
		        //histList[id] = {("SumV_FCAL"+cutString).c_str(), "SumV_FCAL;", "100", "-5", "5", "Entries / 0.1"};
		        histCuts[id] = cutsToApply;
			histVecVals[id] = {SumV_FCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
	       		++id; histVals[id] = {13, applyAccSub, DOCA_FCAL[0]}; 
	       		//histList[id] = {("DOCA_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";DOCA_FCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	       		histCuts[id] = cutsToApply;
			histVecVals[id] = {DOCA_FCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
	       		++id; histVals[id] = {13, applyAccSub, showerQuality_FCAL[0]}; 
	       		//histList[id] = {("showerQuality_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";showerQuality_FCAL").c_str(), "110", "-10", "100", "Entries / 1"};
	       		histCuts[id] = cutsToApply;
			histVecVals[id] = {showerQuality_FCAL};
			//histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
			medBool[0] = cutsToApply*pPhotonInBCAL[0];
			medBool[1] = cutsToApply*pPhotonInBCAL[1];
			medBool[2] = cutsToApply*pPhotonInBCAL[2];
			medBool[3] = cutsToApply*pPhotonInBCAL[3];
		        ++id; histVals[id] = {13, applyAccSub, Energy_BCALPreshower[0]}; 
		        //histList[id] = {("Energy_BCALPreshower"+cutString).c_str(), "Energy_BCALPreshower;", "100", "-5", "5", "Entries / 0.1 GeV"};
		        histCuts[id] = cutsToApply;
			histVecVals[id] = {Energy_BCALPreshower};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
		        ++id; histVals[id] = {13, applyAccSub, Energy_BCAL[0]}; 
		        //histList[id] = {("Energy_BCAL"+cutString).c_str(), "Energy_BCAL;", "100", "-5", "5", "Entries / 0.1 GeV"};
		        histCuts[id] = cutsToApply;
			histVecVals[id] = {Energy_BCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
	       		++id; histVals[id] = {13, applyAccSub, SigLong_BCAL[0]}; 
	       		//histList[id] = {("sigLong_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SigLong_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	       		histCuts[id] = cutsToApply;
			histVecVals[id] = {SigLong_BCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
	       		++id; histVals[id] = {13, applyAccSub, SigTrans_BCAL[0]}; 
	       		//histList[id] = {("sigTrans_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTrans_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	       		histCuts[id] = cutsToApply;
			histVecVals[id] = {SigTrans_BCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
	       		++id; histVals[id] = {13, applyAccSub, SigTheta_BCAL[0]}; 
	       		//histList[id] = {("sigTheta_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTheta_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	       		histCuts[id] = cutsToApply;
			histVecVals[id] = {SigTheta_BCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
	       		++id; histVals[id] = {13, applyAccSub, DeltaPhi_BCAL[0]}; 
	       		//histList[id] = {("DeltaPhi_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";DeltaPhi_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	       		histCuts[id] = cutsToApply;
			histVecVals[id] = {DeltaPhi_BCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
	       		++id; histVals[id] = {13, applyAccSub, DeltaZ_BCAL[0]}; 
	       		//histList[id] = {("DeltaZ_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";DeltaZ_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	       		histCuts[id] = cutsToApply;
			histVecVals[id] = {DeltaZ_BCAL};
			histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3]};
		         
		         // charged showers
		        ++id; histVals[id] = {2, applyAccSub, locE1E9_FCAL_proton}; 
		        //histList[id] = {("E1E9_FCALproton"+cutString).c_str(), "E1E9_FCAL;", "100", "-5", "5", "Entries / 0.1"};
		        histCuts[id] = cutsToApply*inFCAL;
		        ++id; histVals[id] = {2, applyAccSub, locE9E25_FCAL_proton}; 
		        //histList[id] = {("E9E25_FCALproton"+cutString).c_str(), "E9E25_FCAL;", "100", "-5", "5", "Entries / 0.1"};
		        histCuts[id] = cutsToApply*inFCAL;
		        ++id; histVals[id] = {2, applyAccSub, locSumU_FCAL_proton}; 
		        //histList[id] = {("SumU_FCALproton"+cutString).c_str(), "SumU_FCAL;", "100", "-5", "5", "Entries / 0.1"};
		        histCuts[id] = cutsToApply*inFCAL;
		        ++id; histVals[id] = {2, applyAccSub, locSumV_FCAL_proton}; 
		        //histList[id] = {("SumV_FCALproton"+cutString).c_str(), "SumV_FCAL;", "100", "-5", "5", "Entries / 0.1"};
		        histCuts[id] = cutsToApply*inFCAL;
		        ++id; histVals[id] = {2, applyAccSub, locEnergy_BCALPreshower_proton}; 
		        //histList[id] = {("Energy_BCALPreshowerproton"+cutString).c_str(), "Energy_BCALPreshower;", "100", "-5", "5", "Entries / 0.1 GeV"};
		        histCuts[id] = cutsToApply*inBCAL;
		        ++id; histVals[id] = {2, applyAccSub, locEnergy_BCAL_proton}; 
		        //histList[id] = {("Energy_BCALproton"+cutString).c_str(), "Energy_BCAL;", "100", "-5", "5", "Entries / 0.1 GeV"};
		        histCuts[id] = cutsToApply*inBCAL;
	       		++id; histVals[id] = {2, applyAccSub, locSigLong_BCAL_proton}; 
	       		//histList[id] = {("sigLong_BCAL_proton"+cutString).c_str(), ("Cuts="+cutsApplied+";SigLong_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	       		histCuts[id] = cutsToApply*inBCAL;
	       		++id; histVals[id] = {2, applyAccSub, locSigTrans_BCAL_proton}; 
	       		//histList[id] = {("sigTrans_BCAL_proton"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTrans_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	       		histCuts[id] = cutsToApply*inBCAL;
	       		++id; histVals[id] = {2, applyAccSub, locSigTheta_BCAL_proton}; 
	       		//histList[id] = {("sigTheta_BCAL_proton"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTheta_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	       		histCuts[id] = cutsToApply*inBCAL;
		}

		if (outputVanHoveCutEffect){
			for (int j=0; j<3; ++j){
			       if (j==0) { cutsToApply = 1;  applyAccSub=noAccSub;}
			       if (j==1) { cutsToApply = allGeneralCutsPassed; applyAccSub=weight;}
			       if (j==2) { cutsToApply = allGeneralCutsPassed*pVanHove; applyAccSub=weight;}
				++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, mandelstam_tp};
				//histList[id] = {("pi0eta_tVspi0etaMass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "100", "0.8", "1.8", "100" ,"0", "2", "M(#pi_{0}#eta) with Events / 0.01  GeV", "t momentum transfer of #pi_{0}+#eta with Events / 0.02 GeV"};
				histCuts[id] = cutsToApply;
			       	++id; histVals[id] = {11, applyAccSub, mandelstam_tp};
			       	//histList[id] = {("pi0eta_t"+cutString).c_str(), ("Cuts="+cutsApplied+";t momentum transfer of #pi_{0}+#eta").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
			       	histCuts[id] = cutsToApply;

	       			++id; histVals[id] = {8, applyAccSub, vanHove_x, vanHove_y};
	       			//histList[id] = {("vanHove"+cutString).c_str(), ("Cuts="+cutsApplied+";Van Hove Plot").c_str(), "60","-3","3", "60", "-3", "3", "Events / 0.1 degrees", "Events / 0.1 degrees"};
	       			histCuts[id] = cutsToApply;
			}
		}

		for (int j=0; j<4; ++j){
		       if (j==0) { cutsToApply = 1; applyAccSub=noAccSub;}
		       if (j==1) { cutsToApply = allGeneralCutsPassed; applyAccSub=weight;}
		       if (j==2) { cutsToApply = mUE; applyAccSub=weight;}
		       if (j==3) { cutsToApply = mUEChiSq*chiSq100; applyAccSub=weight;}
		       ++id; histVals[id] = {8, applyAccSub, locNumExtraNeutralShowers}; 
		       //histList[id] = {("numExtraNeutralShowers"+cutString).c_str(),("Cuts="+cutsApplied+";Number of extra showers").c_str(),"10","0","10", "Events / 1"};
		       histCuts[id] = cutsToApply;                                                                                                  
		}
		
		std::string massBin;
		if(outputProdPlanePS_BeamAsym){
			for (int i = 0; i < numBinsMass; ++i) {
				applyAccSub = weight;
				// these will always have the general cuts applied, we already have the noCuts version of this (not binned though) 
				cutsToApply = mBeamE*pBeamAsymE;
			        //massBin = "Bin"+std::to_string(i);
				++id; histVals[id] = {8, applyAccSub, locPhi}; 
				//histList[id] = {("prodPlanePSphi"+massBin).c_str(), ";phi binned M(#pi_{0}#eta), General Cuts + 8 GeV < EBeam < 8.7 GeV", "180", "-540","540", "Events / 6 degrees"};
				histCuts[id] = mBeamE*pBeamAsymE*p_phiMassBinned[i];
				++id; histVals[id] = {8, applyAccSub, locPhi, locPi0Eta_Kin}; 
				//histList[id] = {("prodPlanePSphiVsMass"+massBin).c_str(), ";phi vs mass binned in M(#pi_{0}#eta), General Cuts + 8 GeV < EBeam < 8.7 GeV", "180", "-540","540", "300", "0", "3", "Events / 6 degrees", "Events / 0.1 GeV"};
				histCuts[id] = cutsToApply*p_phiMassBinned[i];

				for (int iDelta=0; iDelta<2; ++iDelta){
					++id; histVals[id] = {8, applyAccSub, 1};
					applyAccSub = weight;
					cutsToApply = allGeneralCutsPassed;
					//histList[id] = {"numPassConeDelta","",3,0,3};
					histVecVals[id] = {countCone};
					histCuts[id] = cutsToApply;
					medBool[0] = cutsToApply*pi0_inCone[iDelta]*p_phiMassBinned[i];
					medBool[1] = cutsToApply*eta_inCone[iDelta]*p_phiMassBinned[i];
					medBool[2] = cutsToApply*largeAngle[iDelta]*p_phiMassBinned[i];
					medBool[3] = cutsToApply*withinCone[iDelta]*p_phiMassBinned[i];
					medBool[4] = cutsToApply*!withinCone[iDelta]*p_phiMassBinned[i];
					histVecCuts[id] = {medBool[0], medBool[1], medBool[2], medBool[3],medBool[4]};	
				}
			}
		}
		
		if(outputMandel_tBinnedInMass){
			for (int i = 0; i < numBinsMass_t; ++i){
				// these will always have the general cuts applied, we already have the noCuts version of this (not binned though) 
				applyAccSub = weight;
				cutsToApply = allGeneralCutsPassed;
			        massBin = "Bin"+std::to_string(i);
				++id; histVals[id] = {11, applyAccSub, mandelstam_tp};  
				//histList[id] = {("pi0eta_t"+massBin).c_str(), "t momentum transfer binned in M(#pi_{0}#eta)", "100" ,"-10","0", "Events / 0.1 GeV"};
				histCuts[id] = cutsToApply*p_tMassBinned[i];
				++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, mandelstam_tp};
				//histList[id] = {("pi0eta_tVspi0etaMass"+massBin).c_str(), "", "100", "0.8", "1.8", "100" ,"-2", "0", "M(#pi_{0}#eta) with Events / 0.01  GeV","t momentum transfer of #pi_{0}+#eta with Events / 0.02 GeV"};
				histCuts[id] = cutsToApply*p_tMassBinned[i];
			}
		}
		
		//std::string rfCutString;
		//// we will use i=0 to be the noCut graph. i>0 would be differnet CL cuts applied where i is the order.
		//for (int i = 0; i < 7; i++){
		//	if (i == 0){ cutsToApply = 1; applyAccSub=noAccSub;}
		//	if (i == 1){ cutsToApply = pCLKinFit1*pDiffCL; applyAccSub=weight;}
		//	if (i == 3){ cutsToApply = pCLKinFit3*pDiffCL; applyAccSub=weight;}
		//	if (i == 4){ cutsToApply = pCLKinFit4*pDiffCL; applyAccSub=weight;}
		//	if (i == 5){ cutsToApply = pCLKinFit5*pDiffCL; applyAccSub=weight;}
		//	if (i == 6){ cutsToApply = pCLKinFit6*pDiffCL; applyAccSub=weight;}
		//       	rfCutString = "CutCLOrderNeg"+std::to_string(i);
		//       	++id; histVals[id] = {8, applyAccSub, locDeltaTRF}; 
		//       	//histList[id] = {("RFTime"+rfCutString).c_str(), "; RF Time (ns)", "400", "-6.0", "6.0", "Entries / 0.03 ns"};
		//	histCuts[id] = cutsToApply;
		//}
		
		if (outputUERegion){
			//std::string ueCutString;
			// we will use i=0 to be the noCut graph. i>0 would be differnet CL cuts applied where i is the order.
			for (int i = 0; i < numRegions_UE; i++){
			       //ueCutString = "UEregion"+std::to_string(i);
			       cutsToApply = mEllipseUE_pre*chiSq100; applyAccSub=weight;
			       ++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
			       //histList[id] = {("pi0Mass"+ueCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
			       histCuts[id] = p_pi0MassEtaMassUEregion[i]*cutsToApply;
			       ++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
			       //histList[id] = {("etaMass"+ueCutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
			       histCuts[id] = p_pi0MassEtaMassUEregion[i]*cutsToApply;
			}
		}
    
		if(outputPi0Resolutions){
	    		for (int i = 0; i < numRegions_E; i++){
	    		        //reuse ueCutString to set the region name
	        	        cutsToApply = mEllipse_pre; applyAccSub=weight;
	    		        if (!is_pi0eta){
                	        	++id; histVals[id] = {5, applyAccSub, theta_pi0_lab};
                	        	//histList[id] = {("thetaLabEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";cos(#theta) lab of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
                	        	// these are not the final cuts we apply! Its in the second section
	        	           	histCuts[id] = cutsToApply*p_pi0MassPi0Eregion_1[i]*pSelectf2;
                	        	++id; histVals[id] = {6, applyAccSub, theta_eta_lab};
                	        	//histList[id] = {("tThetaLabEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";cos(#theta) lab of #pi_{0}").c_str(), "100","-1.00","1.00", "Events / 0.02"};
                	        	// these are not the final cuts we apply! Its in the second section
	        	           	histCuts[id] = cutsToApply*p_pi0MassPi0Eregion_2[i]*pSelectf2;


	        			++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
	        			//histList[id] = {("pi0MassKinEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
					// these are not the final cuts we apply! Its in the second section
	        	           	histCuts[id] = cutsToApply*p_pi0MassPi0Eregion_1[i]*pSelectf2;
	        			++id; histVals[id] = {5, applyAccSub, locPi0Mass};
	        			//histList[id] = {("pi0MassEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
					// these are not the final cuts we apply! Its in the second section
	        	           	histCuts[id] = cutsToApply*p_pi0MassPi0Eregion_1[i]*pSelectf2;
	        			++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
	        			//histList[id] = {("pi0MassKinEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	        	           	histCuts[id] = cutsToApply*p_pi0MassPi0Eregion_2[i]*pSelectf2;
	        			++id; histVals[id] = {6, applyAccSub, locEtaMass};
	        			//histList[id] = {("pi0MassEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	        	           	histCuts[id] = cutsToApply*p_pi0MassPi0Eregion_2[i]*pSelectf2;


	       				++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
	       				//histList[id] = {("pi0MassFCALKinEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pPi0InFCAL*p_pi0MassPi0Eregion_1[i]*pSelectf2; 
	       				++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
	       				//histList[id] = {("pi0MassBCALKinEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pPi0InBCAL*p_pi0MassPi0Eregion_1[i]*pSelectf2; 
	       				++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
	       				//histList[id] = {("pi0MassFCALKinEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(),pi0BinRange[0], pi0BinRange[1], pi0BinRange[2] , "Events / 0.0025 GeV"};
	       				histCuts[id] = cutsToApply*pEtaInFCAL*p_pi0MassPi0Eregion_2[i]*pSelectf2; 
	       				++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
	       				//histList[id] = {("pi0MassBCALKinEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.0025 GeV"};
	       				histCuts[id] = cutsToApply*pEtaInBCAL*p_pi0MassPi0Eregion_2[i]*pSelectf2; 

	       				++id; histVals[id] = {5, applyAccSub, locPi0Mass};
	       				//histList[id] = {("pi0MassFCALEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pPi0InFCAL*p_pi0MassPi0Eregion_1[i]*pSelectf2; 
	       				++id; histVals[id] = {5, applyAccSub, locPi0Mass};
	       				//histList[id] = {("pi0MassBCALEregion"+std::to_string(i)+"_1").c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pPi0InBCAL*p_pi0MassPi0Eregion_1[i]*pSelectf2; 
	       				++id; histVals[id] = {6, applyAccSub, locEtaMass};
	       				//histList[id] = {("pi0MassFCALEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(),pi0BinRange[0], pi0BinRange[1], pi0BinRange[2] , "Events / 0.0025 GeV"};
	       				histCuts[id] = cutsToApply*pEtaInFCAL*p_pi0MassPi0Eregion_2[i]*pSelectf2; 
	       				++id; histVals[id] = {6, applyAccSub, locEtaMass};
	       				//histList[id] = {("pi0MassBCALEregion"+std::to_string(i)+"_2").c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.0025 GeV"};
	       				histCuts[id] = cutsToApply*pEtaInBCAL*p_pi0MassPi0Eregion_2[i]*pSelectf2; 
	    		        }
	    		}
		}

		if(outputChiSqRegion){
			//std::string chiSqCutString;
			// we will use i=0 to be the noCut graph. i>0 would be differnet CL cuts applied where i is the order.
			for (int i = 0; i < numRegions_ChiSq; i++){
			       //chiSqCutString = "ChiSqregion"+std::to_string(i);
			       cutsToApply = mEllipseUEChiSq_pre; applyAccSub=weight;
			       ++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
			       //histList[id] = {("pi0Mass"+chiSqCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
			       histCuts[id] = p_pi0MassEtaMassChiSqregion[i]*cutsToApply;
			       ++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
			       //histList[id] = {("etaMass"+chiSqCutString).c_str(), ("Cuts="+cutsApplied+";M(#eta) (GeV)").c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
			       histCuts[id] = p_pi0MassEtaMassChiSqregion[i]*cutsToApply;
			}
		}

		//std::string pi0etaCutString;
		for (int i = 0; i < 5; i++){
			// 0 = noCut; 1 = Gen Cuts + EBeam > 6; 2 = Gen Cuts + BaryBkg
		        if (i==0) { cutsToApply = 1; applyAccSub=noAccSub;}
		        if (i==1) { cutsToApply = allGeneralCutsPassed; applyAccSub=weight;}
		        if (i==2) { cutsToApply = allGeneralCutsPassed*pEtaProtonBaryonCut*ppi0ProtonBaryonCut; applyAccSub=weight;}
	       		if (i==3) { cutsToApply = allGeneralCutsPassed*pVanHove;  applyAccSub=weight;}
	       		if (i==4) { cutsToApply = mEllipse;  applyAccSub=weight;}
		        //pi0etaCutString = "NoCut_AccSubEBeam_AccSubBaryBkg"+std::to_string(i);
			++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
			//histList[id] = {("pi0eta1D"+pi0etaCutString).c_str(), ";M(#pi_{0}#eta) (GeV)", "350", "0", "3.5", "Events / 0.01 GeV"};
			histCuts[id] = cutsToApply;
			//++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
			////histList[id] = {("pi0eta1D8GeVPlus"+pi0etaCutString).c_str(), ";M(#pi_{0}#eta) (GeV)", "350", "0", "3.5", "Events / 0.01 GeV"};
			//histCuts[id] = cutsToApply*pBeamE8GeVPlus;
		}


		if(outputMassBinnedE){
			//Did not merge with the one above just due to the naming convention is different.
			for (int i = 0; i < 3; i++){
				// 0 = noCut; 1 = Gen Cuts + EBeam > 6; 2 = Gen Cuts + BaryBkg
			        if (i==0) { cutsToApply = 1; applyAccSub=noAccSub;}
				// this section does not use pBeamE since this section bins the M(pi0eta) in beam energy. 
			        if (i==1) { cutsToApply = mBeamE; applyAccSub=weight;}
			        if (i==2) { cutsToApply = mBeamE*pEtaProtonBaryonCut*ppi0ProtonBaryonCut; applyAccSub=weight;}
			        //pi0etaCutString = std::to_string(i);
				++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
				//histList[id] = {("pi0eta1D30to46ID"+pi0etaCutString).c_str(), ";M(#pi_{0}#eta) (GeV)", "350", "0", "3.5", "Events / 0.01 GeV", "Beam Energy = 3.0 to 4.6 Gev"};
				histCuts[id] = cutsToApply*pBeamE30to46;
				++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
				//histList[id] = {("pi0eta1D46to62ID"+pi0etaCutString).c_str(), ";M(#pi_{0}#eta) (GeV)", "350", "0", "3.5", "Events / 0.01 GeV","Beam Energy = 4.6 to 6.2 Gev"};
				histCuts[id] = cutsToApply*pBeamE46to62;
				++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
				//histList[id] = {("pi0eta1D62to78ID"+pi0etaCutString).c_str(), ";M(#pi_{0}#eta) (GeV)", "350", "0", "3.5", "Events / 0.01 GeV","Beam Energy = 6.2 to 7.8 Gev"};
				histCuts[id] = cutsToApply*pBeamE62to78;
				++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
				//histList[id] = {("pi0eta1D78to94ID"+pi0etaCutString).c_str(), ";M(#pi_{0}#eta) (GeV)", "350", "0", "3.5", "Events / 0.01 GeV", "Beam Energy = 7.8 to 9.4 Gev"};
				histCuts[id] = cutsToApply*pBeamE78to94;
				++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
				//histList[id] = {("pi0eta1D94to11ID"+pi0etaCutString).c_str(), ";M(#pi_{0}#eta) (GeV)", "350", "0", "3.5", "Events / 0.01 GeV", "Beam Energy = 9.4 to 11 Gev"};
				histCuts[id] = cutsToApply*pBeamE94to11;
			}
		}
		
		if(output_mEllipse_UE_ChiSq){
			for (int i = 0; i<5; i++){
			        if (i == 0){ cutString = ""; cutsToApply = 1; applyAccSub=noAccSub;}
			        if (i == 1){ cutString = "Cut"; cutsToApply = mEllipse_pre; applyAccSub=weight;}
			        if (i == 2){ cutString = "noUECut"; cutsToApply = mEllipseUE_pre; applyAccSub=weight;}
	       			if (i == 3){ cutString = "noUEChiSqCut"; cutsToApply = mEllipseUEChiSq_pre; applyAccSub=weight;}
				// This one should not have any accsub since it would always basicalyl filly with a negative weight due to being in the BKG regional ellipse
	       			if (i == 4){ cutString = "CutBKG"; cutsToApply = mEllipse*pYellowBKG; applyAccSub=noAccSub;}
				++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
				//histList[id] = {("pi0Mass_Kin"+cutString).c_str(), ";M(#pi_{0}) (GeV)", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
				//histList[id] = {("pi0MassFCAL_Kin"+cutString).c_str(), ";M(#pi_{0}) (GeV)", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				histCuts[id] = cutsToApply*pPi0InFCAL; 
				++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
				//histList[id] = {("pi0MassBCAL_Kin"+cutString).c_str(), ";M(#pi_{0}) (GeV)", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				histCuts[id] = cutsToApply*pPi0InBCAL; 
				++id; histVals[id] = {5, applyAccSub, locPi0Mass_Kin};
				//histList[id] = {("pi0MassSplit_Kin"+cutString).c_str(), ";M(#pi_{0}) (GeV)", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				histCuts[id] = cutsToApply*pPi0InSplit; 
				++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
				//histList[id] = {("etaMass_Kin"+cutString).c_str(), ";M(#eta) (GeV)", etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
				//histList[id] = {("etaMassFCAL_Kin"+cutString).c_str(), ";M(#eta) (GeV)",etaBinRange[0], etaBinRange[1], etaBinRange[2] , "Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply*pEtaInFCAL; 
				++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
				//histList[id] = {("etaMassBCAL_Kin"+cutString).c_str(), ";M(#eta) (GeV)", etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply*pEtaInBCAL; 
				++id; histVals[id] = {6, applyAccSub, locEtaMass_Kin};
				//histList[id] = {("etaMassSplit_Kin"+cutString).c_str(), ";M(#eta) (GeV)", etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply*pEtaInSplit; 
				// both pi0eta and pi0pi0 have a kin fit version of the mass plots so we can just put it here
				++id; histVals[id] = {5, applyAccSub, locPi0Mass};
				//histList[id] = {("pi0Mass"+cutString).c_str(), ";M(#pi_{0}) (GeV)", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {5, applyAccSub, locPi0Mass};
				//histList[id] = {("pi0MassFCAL"+cutString).c_str(), ";M(#pi_{0}) (GeV)", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				histCuts[id] = cutsToApply*pPi0InFCAL; 
				++id; histVals[id] = {5, applyAccSub, locPi0Mass};
				//histList[id] = {("pi0MassBCAL"+cutString).c_str(), ";M(#pi_{0}) (GeV)", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				histCuts[id] = cutsToApply*pPi0InBCAL; 
				++id; histVals[id] = {5, applyAccSub, locPi0Mass};
				//histList[id] = {("pi0MassSplit"+cutString).c_str(), ";M(#pi_{0}) (GeV)", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
				histCuts[id] = cutsToApply*pPi0InSplit; 
				++id; histVals[id] = {6, applyAccSub, locEtaMass};
				//histList[id] = {("etaMass"+cutString).c_str(), ";M(#eta) (GeV)", etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {6, applyAccSub, locEtaMass};
				//histList[id] = {("etaMassFCAL"+cutString).c_str(), ";M(#eta) (GeV)",etaBinRange[0], etaBinRange[1], etaBinRange[2] , "Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply*pEtaInFCAL; 
				++id; histVals[id] = {6, applyAccSub, locEtaMass};
				//histList[id] = {("etaMassBCAL"+cutString).c_str(), ";M(#eta) (GeV)", etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply*pEtaInBCAL; 
				++id; histVals[id] = {6, applyAccSub, locEtaMass};
				//histList[id] = {("etaMassSplit"+cutString).c_str(), ";M(#eta) (GeV)", etaBinRange[0], etaBinRange[1], etaBinRange[2], "Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply*pEtaInSplit; 


				// 2D plot
				++id; histVals[id] = {4, applyAccSub, locPi0Mass_Kin, locEtaMass_Kin};
				//histList[id] = {("pi0eta"+cutString).c_str(), ";", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], etaBinRange[0], etaBinRange[1], etaBinRange[2], "#pi_{0} Mass (GeV) with Events / 0.001 GeV","#eta Mass (GeV) with Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply;
				// LOOKING AT THE MISMATCHING OF THE PHOTONS FOR ONLY THE PI0PI0
				if (!is_pi0eta){
	       				++id; histVals[id] = {15, applyAccSub, locPi0Mass_Kin_mismatch};
	       				//histList[id] = {("pi0MassMismatch_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply;
	       				++id; histVals[id] = {15, applyAccSub, locPi0Mass_Kin_mismatch};
	       				//histList[id] = {("pi0MassFCALMismatch_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pPi0InFCAL_mismatch; 
	       				++id; histVals[id] = {15, applyAccSub, locPi0Mass_Kin_mismatch};
	       				//histList[id] = {("pi0MassBCALMismatch_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pPi0InBCAL_mismatch; 
	       				++id; histVals[id] = {15, applyAccSub, locPi0Mass_Kin_mismatch};
	       				//histList[id] = {("pi0MassSplitMismatch_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pPi0InSplit_mismatch; 
	       				++id; histVals[id] = {16, applyAccSub, locEtaMass_Kin_mismatch};
	       				//histList[id] = {("pi0MassMismatch_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply;
	       				++id; histVals[id] = {16, applyAccSub, locEtaMass_Kin_mismatch};
	       				//histList[id] = {("pi0MassFCALMismatch_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(),pi0BinRange[0], pi0BinRange[1], pi0BinRange[2] , "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pEtaInFCAL_mismatch; 
	       				++id; histVals[id] = {16, applyAccSub, locEtaMass_Kin_mismatch};
	       				//histList[id] = {("pi0MassBCALMismatch_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pEtaInBCAL_mismatch; 
	       				++id; histVals[id] = {16, applyAccSub, locEtaMass_Kin_mismatch};
	       				//histList[id] = {("pi0MassSplitMismatch_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}) (GeV)").c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "Events / 0.001 GeV"};
	       				histCuts[id] = cutsToApply*pEtaInSplit_mismatch; 
				}
			}
		}

		if(outputBkgSubOnBkgAndSignalRegions){
			// ********** SECTION TO CHECK HOW THE BKG SUB REMOVES THE BKGS LIKE BARYON RESONANCES
			for (int i=0; i<7; ++i){
			       if (i == 0){ cutsToApply = mEllipse_pre*pYellowBKG; applyAccSub=weightAS;}
			       if (i == 1){ cutsToApply = mEllipse_pre*pinsideEllipse; applyAccSub=weightAS;}
				// The effect of using applyAccSub=weight is basially like only selecting the red and yellow ellipses since the weight is = 0 when not in those regions.
				// Sig_weigthA should be the same as Sig_weightB in terms of hist shape but the number of events might differ since in Sig_weightA we keep the 0-weight events. This can be argued since the || opertor is a union of sets
				// 	and we can include the 0-weight region in for free since the weight=0 there. Recombine the 3 bkgSub regions to get the full set which is always true.
				// BUT!!!! THIS IS NOT THE CASE! THIS IS due to uniqueness tracking. Since the 0-weight region is so large we have a good chance of selecting a combination here and not being able to use that set of particles anymore.
			       if (i == 2){ cutsToApply = mEllipse_pre; applyAccSub=weight;}
			       if (i == 3){ cutsToApply = mEllipse; applyAccSub=weight;}
			       if (i == 4){ cutsToApply = mEllipse*ptLT1; applyAccSub=weight;}
			       if (i == 5){ cutsToApply = mEllipse_pre*pYellowBKG*ptLT1; applyAccSub=weightAS;}
			       if (i == 6){ cutsToApply = mEllipse_pre*pinsideEllipse*ptLT1; applyAccSub=weightAS;}
				if (is_pi0eta){
			       		++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
			       		//histList[id] = {("pi0eta1D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
			       		histCuts[id] = cutsToApply;
			       		++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
			       		//histList[id] = {("pi0proton1D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
			       		histCuts[id] = cutsToApply;
			       		++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin};
			       		//histList[id] = {("etaproton1D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#etaproton) (GeV)").c_str(), "450", "0", "4.5", "Events / 0.01 GeV"};
			       		histCuts[id] = cutsToApply;
			       		++id; histVals[id] = {8, applyAccSub, vanHove_x, vanHove_y};
			       		//histList[id] = {("vanHove"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";Van Hove Plot").c_str(), "60","-3","3", "60", "-3", "3", "Events / 0.1 degrees", "Events / 0.1 degrees"};
			       		histCuts[id] = cutsToApply;
        			        ++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_eta_GJ};
        			        //histList[id] = {("eta_cosTheta_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta)Events / 0.01 GeV", "Cos(#theta) of #eta Events / 0.2"};
        			        histCuts[id] = cutsToApply;
					++id; histVals[id] = {11, applyAccSub, mandelstam_tp};
					//histList[id] = {("pi0eta_t"+cutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using recoil+target").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
					histCuts[id] = cutsToApply;
					++id; histVals[id] = {11, applyAccSub, mandelstam_tp_pe};
					//histList[id] = {("pi0eta_t_pe"+cutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using beam+pi0eta").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
					histCuts[id] = cutsToApply;
				}
				else{
	       				++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin};
	       				//histList[id] = {("pi0pi01D"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}#eta) (GeV)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
	       				histCuts[id] = cutsToApply;
	       				++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
	       				//histList[id] = {("pi0proton1D_1"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#pi_{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
	       				histCuts[id] = cutsToApply;
	       				++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin};
	       				//histList[id] = {("pi0proton1D_2"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";M(#etaproton) (GeV)").c_str(), "450", "0", "4.5", "Events / 0.01 GeV"};
	       				histCuts[id] = cutsToApply;
	       				++id; histVals[id] = {8, applyAccSub, vanHove_x, vanHove_y};
	       				//histList[id] = {("vanHove"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";Van Hove Plot").c_str(), "60","-3","3", "60", "-3", "3", "Events / 0.1 degrees", "Events / 0.1 degrees"};
	       				histCuts[id] = cutsToApply;
                			++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_pi0_GJ};
                			//histList[id] = {("pi0_cosTheta_GJvsM_1"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta)Events / 0.01 GeV", "Cos(#theta) of #eta Events / 0.2"};
                			histCuts[id] = cutsToApply;
                			++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_eta_GJ};
                			//histList[id] = {("pi0_cosTheta_GJvsM_2"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta)Events / 0.01 GeV", "Cos(#theta) of #eta Events / 0.2"};
                			histCuts[id] = cutsToApply;
					++id; histVals[id] = {11, applyAccSub, mandelstam_tp};
					//histList[id] = {("pi0eta_t"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using recoil+target").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
					histCuts[id] = cutsToApply;
					++id; histVals[id] = {11, applyAccSub, mandelstam_tp_pe};
					//histList[id] = {("pi0eta_t_pp"+pi0etaCutString).c_str(), ("Cuts="+cutsApplied+";-t' momentum transfer using beam+pi0eta").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
					histCuts[id] = cutsToApply;
				}
			}
		}



		if(outputCorrelationBetweenMasses){
			// This section will look at the correlation between the M(pi0) vs M(pi0eta) and M(eta) vs M(pi0eta) to see if there exists strong correlation between the variables. This is to see if we
			// can properly use the SPlot method instead of bkg subtraction. SPlot method requires the control variable to not be correlated with the varaibles in the discriminating variable set.
			// we use mEllipse_pre since we want to include all the regions on the M(pi0)vsM(eta) 2D histogram; if we use mEllipse we would already be removing the 0-weight region in between the yellow
			// and red regions.
			if (is_pi0eta){
				cutsToApply = mEllipse_pre; applyAccSub=weightAS;
				++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin, locPi0Mass_Kin};
				//histList[id] = {("pi0etaPi0"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#pi_{0}) Events / 0.001 GeV"};
				histCuts[id] = cutsToApply;
				++id; histVals[id] = {4, applyAccSub, locPi0Eta_Kin, locEtaMass_Kin};
				//histList[id] = {("pi0etaEta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", etaBinRange[0], etaBinRange[1], etaBinRange[2], "M(#pi_{0}#eta) (GeV) with Events / 0.025 GeV", "M(#eta_{0}) Events / 0.0025 GeV"};
				histCuts[id] = cutsToApply;
			}
		}

		if(outputMassShift){
			// Doesnt matter if it is pi0eta or pi0pi0 since the way we will them in this case is indifferent
			cutsToApply = mEllipse_pre; applyAccSub=weightAS;
			++id; histVals[id] = {5, applyAccSub, locPi0Mass_charged};
			//histList[id] = {("pi0Mass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {6, applyAccSub, locEtaMass_charged};
			//histList[id] = {("etaMass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
		       	++id; histVals[id] = {4, applyAccSub, locPi0Mass_charged, locEtaMass_charged};
		       	//histList[id] = {("pi0eta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], etaBinRange[0], etaBinRange[1], etaBinRange[2], "#pi_{0} Mass (GeV) with Events / 0.001 GeV","#eta Mass (GeV) with Events / 0.0025 GeV"};
		       	histCuts[id] = cutsToApply;

			++id; histVals[id] = {5, applyAccSub, locPi0Mass_target};
			//histList[id] = {("pi0Mass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {6, applyAccSub, locEtaMass_target};
			//histList[id] = {("etaMass"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), etaBinRange[0], etaBinRange[1], etaBinRange[2], "M(#pi_{0}) Events / 0.001 GeV"};
			histCuts[id] = cutsToApply;
		       	++id; histVals[id] = {4, applyAccSub, locPi0Mass_target, locEtaMass_target};
		       	//histList[id] = {("pi0eta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), pi0BinRange[0], pi0BinRange[1], pi0BinRange[2], etaBinRange[0], etaBinRange[1], etaBinRange[2], "#pi_{0} Mass (GeV) with Events / 0.001 GeV","#eta Mass (GeV) with Events / 0.0025 GeV"};
		       	histCuts[id] = cutsToApply;
		}

		//if (is_pi0eta){ 
		//	// mEllipse is used here since we are trying to show how the bkg and signal regions from the bkgSub technique changes things. 
		//	cutsToApply = mEllipse*!pMPi0P14; applyAccSub=weight;
        	//        ++id; histVals[id] = {8, applyAccSub, locPi0Eta_Kin, cosTheta_eta_GJ};
        	//        //histList[id] = {("eta_cosTheta_GJvsM"+cutString).c_str(), ("Cuts="+cutsApplied+" cosTheta of #eta vs M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "100","-1.00","1.00", "M(#pi_{0}#eta)Events / 0.01 GeV", "Cos(#theta) of #eta Events / 0.2"};
        	//        histCuts[id] = cutsToApply;
		//	cutsToApply = mEllipse*pMPi0P14; applyAccSub=weight;
		//	++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin};
		//	//histList[id] = {"pi0eta1D_Mppi014", ("Cuts="+cutsApplied+";M(#pi_{0}#eta)").c_str(), "350", "0", "3.5", "Events / 0.01 GeV"};
		//	histCuts[id] = cutsToApply;
		//	++id; histVals[id] = {11, applyAccSub, mandelstam_tp};
		//	//histList[id] = {"pi0eta_t_Mppi014", ("Cuts="+cutsApplied+";t momentum transfer of #pi_{0}+#eta").c_str(), "100" ,"0","10", "Events / 0.1 GeV"};
		//	histCuts[id] = cutsToApply;
		//}
		/// ******************** END SECTION ON BKG SUBTRACTION



		//std::vector<std::string> cutVariations = {"", "Cut", "CutOnlydzR", "CutOnlydzRdEdx", "CutOnlyMMcl"};
		for (int j=0; j<5; ++j){
		        if (j==0) { cutsToApply = 1; applyAccSub=noAccSub;}
		        if (j==1) { cutsToApply = mMagP3; applyAccSub=weight;}
		        if (j==2) { cutsToApply = pzCutmin*pRProton; applyAccSub=noAccSub;}
		        if (j==3) { cutsToApply = pzCutmin*pRProton*pdEdxCDCProton; applyAccSub=noAccSub;}
		        if (j==4) { cutsToApply = pMissingMassSquared*pChiSq; applyAccSub=noAccSub;}
		        //cutString = cutVariations[i];
			++id; histVals[id] = {2, applyAccSub, locMagP3Proton}; 
			//histList[id] = {("MagP3Proton"+cutString).c_str(), ";Proton Momenetum GeV", "200", "0" , "5", "Events / 0.025 GeV"};
			histCuts[id] = cutsToApply;
		}

		if(outputPhotonXandShowerLoc){
			cutsToApply = 1;
			applyAccSub = noAccSub;
			++id; histVals[id] = {13, applyAccSub, photonXs_Kin[0]}; 
			histVecVals[id] = {photonXs_Kin};
			//histList[id] = {"X_Kin","position (cm)", "160", "0", "160", "Events / 1 cm"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {13, applyAccSub, photonYs_Kin[0]}; 
			histVecVals[id] = {photonYs_Kin};
			//histList[id] = {"Y_Kin","position (cm)", "160", "0", "160", "Events / 1 cm"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {13, applyAccSub, photonZs_Kin[0]}; 
			histVecVals[id] = {photonZs_Kin};
			//histList[id] = {"Z_Kin","position (cm)", "250", "0", "1000", "Events / 4 cm"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {13, applyAccSub, photonTs_Kin[0]}; 
			histVecVals[id] = {photonTs_Kin};
			//histList[id] = {"T_Kin","time (s)", "160", "0", "160", "Events / 1 cm"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {13, applyAccSub, photonXs_Shower[0]}; 
			histVecVals[id] = {photonXs_Shower};
			//histList[id] = {"X_Shower","position (cm)", "160", "0", "160", "Events / 1 cm"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {13, applyAccSub, photonYs_Shower[0]}; 
			histVecVals[id] = {photonYs_Shower};
			//histList[id] = {"Y_Shower","position (cm)", "160", "0", "160", "Events / 1 cm"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {13, applyAccSub, photonZs_Shower[0]}; 
			histVecVals[id] = {photonZs_Shower};
			//histList[id] = {"Z_Shower","position (cm)", "160", "0", "160", "Events / 1 cm"};
			histCuts[id] = cutsToApply;
			++id; histVals[id] = {13, applyAccSub, photonTs_Shower[0]}; 
			histVecVals[id] = {photonTs_Shower};
			//histList[id] = {"T_Shower","time (s)", "160", "0", "160", "Events / 1 cm"};
			histCuts[id] = cutsToApply;
		}


		//if(showOutput){
		//	cout << "**** RECHECK GROUPS *****" << endl << "***********************" << endl;
		//	for(std::size_t i = 0; i<vec_group_ids.size(); ++i){
		//		if(showOutput) {cout  << " hist Ids for a group: " << std::to_string(i+1) << " - " << groupNames[i] << endl << "=============================" << endl;}
		//		for(std::size_t j = 0; j<vec_group_ids[i].size(); ++j){
		//			if(showOutput){cout << histList[vec_group_ids[i][j]][0] << " with ID:" << histVals[vec_group_ids[i][j]][0]<<endl;}
		//		}
		//	}
		//}

// ******************************************************************************************************************************
		 //FIlling the histograms!!! 
// ******************************************************************************************************************************

                //locUsedThisCombo_Pi0Mass is a map which can be inserted into locUsedSoFar_Pi0Mass which is a set of maps with the same form. locUsedThisCombo_Pi0Mass is a map with map<key_type,data_type> 
                //where key_type is Particle_t which is a file inside the GlueX software which labels all particles with a number i.e. Gamma=1 and Proton=14. The data_type is a set of integers to hold the 
                //unique particle ids. When we do this uniqueness tracking we insert locUsedThisCombo_Pi0Mass into locUsedSoFar_Pi0Mass, so we are keeping track of entire maps of particle ids. As long as the
                //entire map is not the same, as in the case of containing a new particle, this condition is passed. Weird how this is even needed since if the Process function just loops over all combos 
                //this condition would be passed every time, maybe it is just a check? 
                // ANSWER : This does actually matter since each combo can have the same particles paired with different beam photons. This actually becomes a problem when we are deciding on when to fill the
                // locUsedSoFar_Pi0Mass with locUsedThisCombo_Pi0Mass. 

                //When filling differnet histograms we still need to watch out for the using the same tracked particle ID. i.e. if we have a 5 photon shower 12345 and the first time through we make a combo with
                //1234 and then after we use 2345. photons 234 will have more filled more times for the same event, giving it some more weight. Atleast for invariant mass calculations this is not a problem and
                //is needed or else we might be biasing. But, for like photonAngle it would be a problem. And obviously there could be multiple charged tracks that could be labeled a proton and multiple tagged
                //photons that could be in each RF bunch.  


                if(showOutput){cout << endl << endl << "\n\nFilling histograms with uniqueness tracking!" << endl << "******************************************" << endl;}
		//int idxCut = 0; 
		//	
		//if(vecLocUsedSoFar_Pi0Mass[idxCut].find(locUsedThisCombo_Pi0Mass) == vecLocUsedSoFar_Pi0Mass[idxCut].end())
		//{
		//	// EACH COMBO OF 1 PROTON + 4 PHOTONS CAN BE PAIRED WITH 1 BEAM PHOTON. LETS SAY IF THERE ARE 2 BEAM PHOTONS WHERE ONE OF THEM WILL HAVE COMBO PARAMETERS THAT PASS THE CUT AND THE OTHER ONE DOES (I.E. DIFF
		//	// UNUSED ENERGIES OR CONFIDENCE LEVELS MAYBE?) IF WE WERE TO SELECT THE ONE THAT DOESNT PASS THE CUTS AND WE INSERT THE COMBO FOR UNIQUENESS TRACKING THEN WHEN WE CHECK THE COMBO WITH THE OTHER PHOTON
		//	// THEN WE WONT GET A CHANCE TO ACCEPT IT SINCE THE COMBO WAS ALREADFY USED. 
		//	if(allGeneralCutsPassed){ 
		//		++count_allGeneralCutsPassedPlusTracked; 
		//        	// FILL IN UNIQUENESS TRACKING FOR ALL VARIATION OF CUTS
		//        	vecLocUsedSoFar_Pi0Mass[idxCut].insert(locUsedThisCombo_Pi0Mass);
		//	}
		//}
		
		// *********************** DOES NOT REQUIRE DIFFERNET UNIQUENESS TRACKING ******************//
		// loop over all the histograms
		groupVec_it=0;	
		if(showOutput){cout << "\n\nFilling usedBeam hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedBeam[hist_id].find(locBeamID) == usedBeam[hist_id].end()) {
				if(histCuts[hist_id]){
		        		usedBeam[hist_id].insert(locBeamID);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		groupVec_it=1;
		if(showOutput){cout << "\n\nFilling usedP1 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedP1[hist_id].find(locProtonTrackID) == usedP1[hist_id].end()) {
				if(histCuts[hist_id]){
		        		usedP1[hist_id].insert(locProtonTrackID);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		groupVec_it=2;
		if(showOutput){cout << "\n\nFilling usedPh12B1 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh12B1[hist_id].find(usingPh12B1) == usedPh12B1[hist_id].end()) {
				if(histCuts[hist_id]){
		        		usedPh12B1[hist_id].insert(usingPh12B1);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] <<" with id: " << std::to_string(hist_id) <<endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		groupVec_it=3;
		if(showOutput){cout << "\n\nFilling usedPh1234 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh1234[hist_id].find(usingPh1234) == usedPh1234[hist_id].end()) {
                                // We include this counter here since this will track our uniqueFinalState ids to be used in q-value subtracting
                                ++finalStateComboID;
				if(histCuts[hist_id]){
		        		usedPh1234[hist_id].insert(usingPh1234);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		groupVec_it=4;
		if(showOutput){cout << "\n\nFilling usedPh12 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh12[hist_id].find(usingPh12) == usedPh12[hist_id].end()) {
				if(histCuts[hist_id]){
					usedPh12[hist_id].insert(usingPh12);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		
		groupVec_it=5;
		if(showOutput){cout << "\n\nFilling usedPh34 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh34[hist_id].find(usingPh34) == usedPh34[hist_id].end()) {
				if(histCuts[hist_id]){
					usedPh34[hist_id].insert(usingPh34);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		
		groupVec_it=6;
		if(showOutput){cout << "\n\nFilling usedPh34B1 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh34B1[hist_id].find(usingPh34B1) == usedPh34B1[hist_id].end()) {
				if(histCuts[hist_id]){
					usedPh34B1[hist_id].insert(usingPh34B1);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		
		groupVec_it=7;
		if(showOutput){cout << "\n\nFilling usedCombo hists - contains photon timings so we will track the entire combo and fill every photon's timing no matter what" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedCombo[hist_id].find(usingCombo) == usedCombo[hist_id].end()) {
				// here we check to see if histVecCuts is non zero so that we will histVecVals[i] with the cut histVecCuts[i]
				if(histCuts[hist_id]){
					usedCombo[hist_id].insert(usingCombo);
					// if not filling a vector 
					if (histVecVals[hist_id].size()==0){
						if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
						// since most hists are 1D we will verify against that. Else it is 2D and use the else condition....
						if (histVals[hist_id].size()==3){
							dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
						}
						else {
							dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
						}
					}
					// if filling 1D vector into hists
					else if (histVecVals[hist_id].size()==1){
						if(histVecCuts[hist_id].size()==0){ // all of the vector of values to fill will share the same cut condtion which was defined in previous condition with histCutsn
							if(showOutput){cout << "-- IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM WITH ENTIRE VECTOR: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
							for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
								dHist_all1DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals],histVals[hist_id][1]);	
							}
						}
						else { //it should be that in this case the vector of values to fill should have their won cut condition
							for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
								if(histVecCuts[hist_id][itHistVecVals]){
									if(showOutput){cout << "-- IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM WITH A SINGLE ELEMENT OF THE VECTOR: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
									dHist_all1DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals],histVals[hist_id][1]);	
								}
							}
						
						}
					}
					// if filling 2D vector into hists
					else if (histVecVals[hist_id].size()==2){
						if(histVecCuts[hist_id].size()==0){ // all of the vector of values to fill will share the same cut condtion which was defined in previous condition with histCuts
							if(showOutput){cout << "-- IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM WITH ENTIRE VECTOR: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
							for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
								dHist_all2DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals], histVecVals[hist_id][1][itHistVecVals], histVals[hist_id][1]);	
							}
						}
						else { //it should be that in this case the vector of values to fill should have their won cut condition
							for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
								if(histVecCuts[hist_id][itHistVecVals]){
									if(showOutput){cout << "-- IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM WITH ENTIRE VECTOR: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
									dHist_all2DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals], histVecVals[hist_id][1][itHistVecVals], histVals[hist_id][1]);	
								}
							}
						}
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		
		groupVec_it=8;
		if(showOutput){cout << "\n\nFilling usedPh12P1 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh12P1[hist_id].find(usingPh12P1) == usedPh12P1[hist_id].end()) {
				if(histCuts[hist_id]){
					usedPh12P1[hist_id].insert(usingPh12P1);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		
		groupVec_it=9;
		if(showOutput){cout << "\n\nFilling usedPh34P1 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh34P1[hist_id].find(usingPh34P1) == usedPh34P1[hist_id].end()) {
				if(histCuts[hist_id]){
					usedPh34P1[hist_id].insert(usingPh34P1);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		
		groupVec_it=10;
		if(showOutput){cout << "\n\nFilling usedPh1234P1 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh1234P1[hist_id].find(usingPh1234P1) == usedPh1234P1[hist_id].end()) {
				if(histCuts[hist_id]){
					usedPh1234P1[hist_id].insert(usingPh1234P1);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		
		groupVec_it=12;
		if(showOutput){cout << "\n\nFilling neutral track hists" << endl << "-----------------------------------" << endl;}
		// ****************************** TRACKING NEUTRAL TRACKS **********************************//
		beingUsedNeutralIds.push_back(locPhoton1NeutralID);
		beingUsedNeutralIds.push_back(locPhoton2NeutralID);
		beingUsedNeutralIds.push_back(locPhoton3NeutralID);
		beingUsedNeutralIds.push_back(locPhoton4NeutralID);
		//for (int whichPhoton = 0; whichPhoton<4; ++whichPhoton){  // whichPhoton = 0 actually corresponds to photon 1
		// grab appropriate value to fill
		//locPhotonTheta = photonThetas[whichPhoton];
		//if(showOutput) { cout << "=======\nlocPhotonTheta = " << std::to_string(locPhotonTheta) << endl; }
		//locPhotonEnergy= photonEnergies[whichPhoton];
		//if(showOutput) { cout << "locPhotonEnergy = " << std::to_string(locPhotonEnergy) << endl; }
		//locPhotonPhi = photonPhis[whichPhoton];
		//if(showOutput) { cout << "locPhotonPhi = " << std::to_string(locPhotonPhi) << endl; }
		//locPhotonX = photonXs_Shower[whichPhoton];
		//if(showOutput) { cout << "locPhotonX = " << std::to_string(locPhotonX) << endl; }
		//locPhotonY = photonYs_Shower[whichPhoton];
		//if(showOutput) { cout << "locPhotonY = " << std::to_string(locPhotonY) << endl; }
		//locPhotonZ = photonZs_Shower[whichPhoton];
		//if(showOutput) { cout << "locPhotonZ = " << std::to_string(locPhotonZ) << endl; }
		//locPhotonT = photonTs_Shower[whichPhoton];
		//if(showOutput) { cout << "locPhotonT = " << std::to_string(locPhotonT) << endl; }
		//locPhotonX = photonXs_Shower[whichPhoton];
		//if(showOutput) { cout << "locPhotonX = " << std::to_string(locPhotonX) << endl; }
		//locPhotonY = photonYs_Shower[whichPhoton];
		//if(showOutput) { cout << "locPhotonY = " << std::to_string(locPhotonY) << endl; }
		//locPhotonZ = photonZs_Shower[whichPhoton];
		//if(showOutput) { cout << "locPhotonZ = " << std::to_string(locPhotonZ) << endl; }
		//locPhotonT = photonTs_Shower[whichPhoton];
		//if(showOutput) { cout << "locPhotonT = " << std::to_string(locPhotonT) << endl; }
		//locPhotonDeltaT = photonDeltaTs[whichPhoton];
		//if(showOutput) { cout << "locPhotonDeltaT = " << std::to_string(locPhotonDeltaT) << endl << endl;} 
		//locPhotonDetectedSys = photonDetectedSyss[whichPhoton];

		// loop over all the histograms
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if(histCuts[hist_id]){ // There is a single cut for each of these phN histograms 
				if(showOutput){cout << "---> PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
				// since most hists are 1D we will verify against that first
				if (histVecVals[hist_id].size()==1){
					if(showOutput){cout << "Attempting to filling above hist with IDs beingUsedNeutralId: "; }
					for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
						if(showOutput){cout << std::to_string(beingUsedNeutralIds[itHistVecVals]);}
						if (usedPhN[hist_id].find(beingUsedNeutralIds[itHistVecVals])==usedPhN[hist_id].end()){
							if(showOutput){cout<<" is UNIQUE ";}	
						        usedPhN[hist_id].insert(beingUsedNeutralIds[itHistVecVals]); //we get a iterator which references the element of the set so we need to dereference.              
							if(histVecCuts[hist_id].size()==0){
								if(showOutput){cout<<"AND PASSED WITH NO SUBCUT, ";}
								dHist_all1DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals],histVals[hist_id][1]);	
							}
							else{
								if(histVecCuts[hist_id][itHistVecVals]){	
									if(showOutput){cout<<"AND PASSED WITH A SUBCUT, ";}
									dHist_all1DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals],histVals[hist_id][1]);	
								}
								else { if(showOutput){cout<<"AND FAILED WITH A SUBCUT, ";} } 
							}
						}
						else { if(showOutput){ cout << " is not UNIQUE, ";} }
					}
					if(showOutput){cout << "\nwhere only UNIQUE selections are filled and uniqueness sets\n";}
				}
				else if(histVecVals[hist_id].size()==2){ // 2D histograms
					if(showOutput){cout << "Attempting to filling above hist with IDs beingUsedNeutralId: ";}
					for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
						if(showOutput){cout << std::to_string(beingUsedNeutralIds[itHistVecVals]);}
						if (usedPhN[hist_id].find(beingUsedNeutralIds[itHistVecVals])==usedPhN[hist_id].end()){
							if(showOutput){cout<<" is UNIQUE, ";}	
						        usedPhN[hist_id].insert(beingUsedNeutralIds[itHistVecVals]); //we get a iterator which references the element of the set so we need to dereference.              
							if(histVecCuts[hist_id].size()==0){
								if(showOutput){cout<<"AND PASSED WITH NO SUBCUT, ";}
								dHist_all2DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals], histVecVals[hist_id][1][itHistVecVals], histVals[hist_id][1]);	
							}
							else{
								if(histVecCuts[hist_id][itHistVecVals]){	
									if(showOutput){cout<<"AND PASSED WITH A SUBCUT, ";}
									dHist_all2DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals], histVecVals[hist_id][1][itHistVecVals], histVals[hist_id][1]);	
								}
								else { if(showOutput){cout<<"AND FAILED WITH A SUBCUT, ";} }
							}
						}
						else { if(showOutput){ cout << " is not UNIQUE, ";} }
					}
					if(showOutput){cout << "\nwhere only UNIQUE selections fills histograms and uniqueness sets\n";}
				}
				else { if(showOutput) {cout << "histVecVals not the right size since this group is only for photons, all values should be vectors" << endl; } }
			}
			else { if(showOutput){cout << "---> DID NOT PASS CUT, NOT FILLING HISTOGRAM NOR UNIQUENESS TRACKING: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;} } 
		}
		//} /// closes the for loop for the uniqueness of single photons
		//
		//
		if(showOutput){ cout << "\n\nFilling pairwise quantities between photons, * * * ALREADY PRE-UNIQUENESS TRACKED * * *\n-----------------------------" << endl;}
		groupVec_it = 11;
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			//if(showOutput){cout << "groupVec_it, groupHist_it: " << std::to_string(groupVec_it) << ", " << std::to_string(groupHist_it) << endl;}
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if(showOutput){ cout << "Attempting to fill: " + histList[hist_id][0] + " using photonFCALPIDs_ijVec of size: " + std::to_string(photonFCALPIDs_ijVec[hist_id].size())<< endl;}
			for (std::size_t vec_idx=0; vec_idx<photonFCALPIDs_ijVec[hist_id].size(); ++vec_idx){
				if(histCuts[hist_id]){
					dij3 = dij3VecFCAL[hist_id][vec_idx];
					if(showOutput){cout << "	PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVecVals[hist_id].size()==1){
						for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
							dHist_all1DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals],histVals[hist_id][1]);	
						}
					}
					else {
						for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
							dHist_all2DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals], histVecVals[hist_id][1][itHistVecVals], histVals[hist_id][1]);	
						}
					}
				}
				else {
					if(showOutput){cout << "	DID NOT PASS CUTS, NOT FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
				}
				
			}
		}
		groupVec_it = 13;
		// Fill the histogram
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if(showOutput){ cout << "Attempting to fill: " + histList[hist_id][0] + " using photonBCALPIDs_ijVec of size: " + std::to_string(photonBCALPIDs_ijVec[hist_id].size())<< endl;}
			for (std::size_t vec_idx=0; vec_idx<photonBCALPIDs_ijVec[hist_id].size(); ++vec_idx){
				// since we know the pairs are already unique we dont have to track anymore..
				if(histCuts[hist_id]){
					angle_ij = angle_ijVec[hist_id][vec_idx];
					deltaZ_ij = deltaZ_ijVec[hist_id][vec_idx];
					deltaPhi_ij = deltaPhi_ijVec[hist_id][vec_idx];
					if(showOutput){cout << "	PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVecVals[hist_id].size()==1){
						for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
							dHist_all1DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals],histVals[hist_id][1]);	
						}
					}
					else {
						for (std::size_t itHistVecVals=0; itHistVecVals<histVecVals[hist_id][0].size(); ++itHistVecVals) {
							dHist_all2DHists[hist_id]->Fill(histVecVals[hist_id][0][itHistVecVals], histVecVals[hist_id][1][itHistVecVals], histVals[hist_id][1]);	
						}
					}
				}
				else {
					if(showOutput){cout << "	DID NOT PASS CUTS, NOT FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
				}
			}
		}
		groupVec_it=14;
		if(showOutput){cout << "\n\nFilling usedPh13 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh13[hist_id].find(usingPh13) == usedPh13[hist_id].end()) {
				if(histCuts[hist_id]){
					usedPh13[hist_id].insert(usingPh13);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}
		
		groupVec_it=15;
		if(showOutput){cout << "\n\nFilling usedPh24 hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if (usedPh24[hist_id].find(usingPh24) == usedPh24[hist_id].end()) {
				if(histCuts[hist_id]){
					usedPh24[hist_id].insert(usingPh24);
					if(showOutput){cout << "IS UNIQUE AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
					// since most hists are 1D we will verify against that.
					if (histVals[hist_id].size()==3){
						dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
					}
					else {
						dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
					}
				}
				else{ if(showOutput){cout << "IS UNIQUE BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
			}
			else { if(showOutput){ cout << "Not unique so no filling of histograms nor track set, continuing..." << endl; }} 
		}

		groupVec_it=16;
		if(showOutput){cout << "\n\nFilling no uniqueness tracking hists" << endl << "-----------------------------------" << endl;}
		for(std::size_t groupHist_it=0; groupHist_it<vec_group_ids[groupVec_it].size(); ++groupHist_it){
			hist_id = vec_group_ids[groupVec_it][groupHist_it];
			if(histCuts[hist_id]){
				if(showOutput){cout << "NO UNIQUE TRACKING AND PASSED CUTS, FILLING HISTOGRAM: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) << endl;}
				// since most hists are 1D we will verify against that.
				if (histVals[hist_id].size()==3){
					dHist_all1DHists[hist_id]->Fill(histVals[hist_id][2],histVals[hist_id][1]);	
				}
				else {
					dHist_all2DHists[hist_id]->Fill(histVals[hist_id][2], histVals[hist_id][3], histVals[hist_id][1]);	
				}
			}
			else{ if(showOutput){cout << "NO UNIQUE TRACKING BUT DID NOT PASS, NOT FILLING HISTOGRAM AND TRACKING SET: " + histList[hist_id][0] << " @ group:" << histVals[hist_id][0] << " with id: " << std::to_string(hist_id) <<endl;} } 
		}
                // Cut out this combo from the output tree.
                //if (!allGeneralCutsPassed) { dComboWrapper->Set_IsComboCut(true); continue; }

                // METHOD 1: We want to remove the elliptical cut (which doesn't affect our histograms since we filled then already) so that we have a flat tree that we can pass to another program
                // 	that will do the bkg subtraction. 
                // METHOD 2: Here we will apply a loose elliptical cut so that the data that we pass into SPlot has a better form that we can fit to
                if (!mEllipse_pre*pinsideEllipse_loose) { 
			dComboWrapper->Set_IsComboCut(true); continue; 
		}

		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

		/*
		//FILL ANY CUSTOM BRANCHES FIRST!!
		Int_t locMyInt_Flat = 7;
		dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
		{
			dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
			TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
			dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
		}
		*/

                dFlatTreeInterface->Fill_Fundamental<Double_t>("AccWeight", weight);
                dFlatTreeInterface->Fill_Fundamental<Double_t>("uniqueComboID", uniqueComboID);
                dFlatTreeInterface->Fill_Fundamental<Int_t>("finalStateComboID", finalStateComboID);

		++uniqueComboID;
		if (is_pi0eta){
                	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0", locPi0Mass_Kin);
                	dFlatTreeInterface->Fill_Fundamental<Double_t>("Meta", locEtaMass_Kin);
                	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0eta", locPi0Eta_Kin);
		}
		else{
                	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0", locPi0Mass_Kin);
                	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0", locEtaMass_Kin);
                	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0pi0", locPi0Eta_Kin);
		}
                dFlatTreeInterface->Fill_Fundamental<Double_t>("mandelstam_tp", mandelstam_tp);

		// Introduce some angles to use in the phase space distance calculation
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("cosTheta_X_cm", cosTheta_pi0eta_CM); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("phi_X_cm", phi_pi0eta_CM); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("cosTheta_eta_gj", cosTheta_pi0_GJ); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("phi_eta_gj", phi_eta_GJ); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("cosThetaHighestEphotonIneta_gj", cosTheta_largestEinEta_GJ); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("cosThetaHighestEphotonInpi0_cm", cosTheta_largestEinPi0_CM); //fundamental = char, int, float, double, etc.
                dFlatTreeInterface->Fill_Fundamental<Double_t>("vanHove_x",vanHove_x);
                dFlatTreeInterface->Fill_Fundamental<Double_t>("vanHove_y",vanHove_y);
                dFlatTreeInterface->Fill_Fundamental<Double_t>("vanHove_omega",omega);
                // If we were to do Q-Values for the pi0pi0 system we will probably do the same thing and use one of the pi0 as the discriminator varible and check against the other
                dFlatTreeInterface->Fill_Fundamental<Double_t>("pi0_energy", locPi0E_Kin );
		//FILL FLAT TREE
                //dComboWrapper->Set_ComboIndex(loc_i);  // Combo succeeded // this might be redundant since this is after the continue command in the cut condition above, combo must have succeeded already. But Elton has it! 
                Fill_FlatTree(); //for the active combo
	} // end of combo loop
	++eventIdx;
        if(showOutput){cout << "\n\n **************** Finishing the combo loop ***************\n**********************************************************\n" << endl;}


	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
/*
	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEnergy = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/
	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
/*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}

	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
/*
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != "")
		Fill_OutputTree();
*/

        // Even though we have this Set_IsComboCut and continue, that only "continues" the combo loop. This is outside that loop so we have to loop over all the combos again and check to see if any of the combos have 
        // been cut. We only the fill the tree if at least one combo succeeded. If that happens it breaks the loop and begins to fill the output (if there is a file you want to fill it in. 
        Bool_t locIsEventCut = true;
        for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
                //Set branch array indices for combo and all combo particles
                dComboWrapper->Set_ComboIndex(loc_i);
                // Is used to indicate when combos have been cut
                if(dComboWrapper->Get_IsComboCut())
                  continue;
                locIsEventCut = false; // At least one combo succeeded                                                     
                break;
        }
        if(!locIsEventCut && dOutputTreeFileName != ""){ Fill_OutputTree(); }

	//}//closes the //if(itersToRun) condition
        return kTRUE; // this return should close the process loop to return false as the kTrue as the output.
}// end of process loop

void DSelector_ver20::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.
        if(showOutput){ std::cout << "Runs Included In Analysis:" << endl; }
        for (std::set<UInt_t>::iterator it=usedRuns.begin(); it!=usedRuns.end(); ++it){
                if (showOutput){ cout << *it << endl; }
        }

	//if(showOutput){cout << "Num passed UnusedEnergy: " << std::to_string(count_UnusedEnergy)<<endl;}
        //if(showOutput){cout << "Num passed CLKinFit: " << std::to_string(count_ChiSq)<<endl;}
        //if(showOutput){cout << "Num passed DeltaTRF: " << std::to_string(count_DeltaTRF)<<endl;}
        //if(showOutput){cout << "Num passed dij3pass: " << std::to_string(count_dij3pass)<<endl;}
        //if(showOutput){cout << "Num passed PhotonE: " << std::to_string(count_PhotonE)<<endl;}
        //if(showOutput){cout << "Num passed PhotonTheta: " << std::to_string(count_PhotonTheta)<<endl;}
        //if(showOutput){cout << "Num passed MagP3Proton: " << std::to_string(count_MagP3Proton)<<endl;}
        //if(showOutput){cout << "Num passed zCutmin: " << std::to_string(count_zCutmin)<<endl;}
        //if(showOutput){cout << "Num passed RProton: " << std::to_string(count_RProton)<<endl;}
        //if(showOutput){cout << "Num passed MissingMassSquared: " << std::to_string(count_MissingMassSquared)<<endl;}
        //if(showOutput){cout << "Num passed dEdxCDCProton: " << std::to_string(count_dEdxCDCProton)<<endl;}
        //if(showOutput){cout << "Num passed insideEllipse: " << std::to_string(count_insideEllipse)<<endl;}
	//if(showOutput){cout << "The next two show how the numbers are reduced when doing uniqueness tracking" << endl;}
        //if(showOutput){cout << "Num passed allGeneralCutsPassed: " << std::to_string(count_allGeneralCutsPassed)<<endl;}
        ////if(showOutput){cout << "Num passed allGeneralCutsPassedPlusTracked: " << std::to_string(count_allGeneralCutsPassedPlusTracked)<<endl;}

	if(true){cout << "Num passed ShowerQuality: " << std::to_string(count_ShowerQuality)<<endl; }
	if(true){cout << "Num passed BeamE8GeVPlus: " << std::to_string(count_BeamE8GeVPlus)<<endl; }
	if(true){cout << "Num passed UnusedEnergy: " << std::to_string(count_UnusedEnergy)<<endl;}
        if(true){cout << "Num passed CLKinFit: " << std::to_string(count_ChiSq)<<endl;}
        if(true){cout << "Num passed DeltaTRF: " << std::to_string(count_DeltaTRF)<<endl;}
        if(true){cout << "Num passed dij3pass: " << std::to_string(count_dij3pass)<<endl;}
        if(true){cout << "Num passed PhotonE: " << std::to_string(count_PhotonE)<<endl;}
        if(true){cout << "Num passed PhotonTheta: " << std::to_string(count_PhotonTheta)<<endl;}
        if(true){cout << "Num passed MagP3Proton: " << std::to_string(count_MagP3Proton)<<endl;}
        if(true){cout << "Num passed zCutmin: " << std::to_string(count_zCutmin)<<endl;}
        if(true){cout << "Num passed RProton: " << std::to_string(count_RProton)<<endl;}
        if(true){cout << "Num passed MissingMassSquared: " << std::to_string(count_MissingMassSquared)<<endl;}
        if(true){cout << "Num passed dEdxCDCProton: " << std::to_string(count_dEdxCDCProton)<<endl;}
        if(true){cout << "Num passed insideEllipse: " << std::to_string(count_insideEllipse)<<endl;}
	if(true){cout << "The next two show how the numbers are reduced when doing uniqueness tracking" << endl;}
        if(true){cout << "Num passed allGeneralCutsPassed: " << std::to_string(count_allGeneralCutsPassed)<<endl;}
	// we can only use the below code when we are using setupTest.sh. DOesnt work with proof since it will probably try to do this for every thread...
	//dHist_Cuts->SetStats(0);
	//dHist_Cuts->SetCanExtend(TH1::kAllAxes);
	//dHist_Cuts->LabelsDeflate("X");
	//dHist_Cuts->LabelsOption("v");

	//CALL THIS LAST
        if(showOutput){cout << "finalizing!" << endl;}
	DSelector::Finalize(); //Saves results to the output file
}



