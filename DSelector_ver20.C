#include "DSelector_ver20.h"
bool NoCut=0;
// degXXX where XXX = {000,045,090,135,All} where All is polarization independent. Actually anything other than the first 4 cases work but
// MUST BE ATLEAST 3 CHARACTERS LONG.
//string degAngle = "a0a2a2pi1_";
string degAngle="pi0eta_data";
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

string selectDetector="ALL";

void DSelector_ver20::Init(TTree *locTree)
{
        cout << "STARTING" << endl;
	targetCenter = {0,0,65};

        if (is_pi0eta) {
                lowMass = 0.7;
                upMass = 2.0;
                etaProtonBaryonCut = 1.65;
                pi0ProtonBaryonCut = 2;
		binScale = (upMass-lowMass)/numBinsMass;

		// using the meas data
		//ellipseX = 0.134285; ellipseY = 0.540383; ellipseXr = 0.02308917; ellipseYr = 0.0699516;  

                //ellipseX = 0.134547; ellipseY = 0.541950; ellipseXr = 0.025449; ellipseYr = 0.069267;

                //using the kin data
                ellipseX = 0.135881; ellipseY = 0.548625; ellipseXr = 3*0.0076; ellipseYr = 3*0.0191;
//                ellipseXBS1 = 0.135881; ellipseYBS1 = 0.548625; ellipseXrBS1 = 0.022; ellipseYrBS1 = 0.06;
//                ellipseXBS2 = 0.135881; ellipseYBS2 = 0.548625; ellipseXrBS2 = 0.045; ellipseYrBS2 = 0.165;
//		//ellipseXr_loose=0.0391; ellipseYr_loose=0.131;
//		//ellipseXr_loose=0.0391; ellipseYr_loose=0.14;
//		ellipseXr_loose=0.02; ellipseYr_loose=100;
//		areaRatio = 0.067; //double checked these values;

		ellipseYr += 0.02;
		double ellipseXYratio = ellipseYr/ellipseXr;
		double scaleUp = 3;
		ellipseXBS1 = 0.135881; ellipseYBS1 = 0.548625; ellipseXrBS1 = 3*0.0076; ellipseYrBS1 = 3*0.0076*ellipseXYratio;	
		double outerXVal = TMath::Sqrt( scaleUp*ellipseXr*ellipseXr + ellipseXrBS1*ellipseXrBS1 );
		ellipseXBS2 = 0.135881; ellipseYBS2 = 0.548625; ellipseXrBS2 = outerXVal ; ellipseYrBS2 = outerXVal*ellipseXYratio;	
		areaRatio=1/scaleUp;

		skipY = ellipseYrBS1 - ellipseYr;
		skipX = ellipseXrBS1 - ellipseXr;

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
        dOutputFileName = degAngle+"_DSelector_output.root"; //"" for none
        dOutputTreeFileName = degAngle+"_tree_DSelector.root"; //"" for none
        dFlatTreeFileName = degAngle+"_treeFlat_DSelector.root"; //output flat tree (one combo per tree entry), "" for none
        dFlatTreeName = degAngle+"_tree_flat"; //if blank, default name will be chosen

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

	dHistThrownTopologies = new TH1F("hThrownTopologies","hThrownTopologies", 10, -0.5, 9.5);

	vector<TString> locThrownTopologies;
	locThrownTopologies.push_back("4#gammap[#pi^{0},#eta]");
	locThrownTopologies.push_back("6#gammap[3#pi^{0}]");		
	locThrownTopologies.push_back("5#gammap[2#pi^{0},#omega]");
	locThrownTopologies.push_back("6#gammap[2#pi^{0},#eta]");
	locThrownTopologies.push_back("4#gammap[2#pi^{0}]");
	locThrownTopologies.push_back("8#gammap[4#pi^{0},#eta]");
	locThrownTopologies.push_back("6#gamma#pi^{#plus}#pi^{#minus}p[3#pi^{0}]");		
	locThrownTopologies.push_back("4#gamma#pi^{#plus}#pi^{#minus}p[2#pi^{0}]");
	locThrownTopologies.push_back("4#gamma#pi^{#plus}#pi^{#minus}p[#pi^{0},#eta]");		
	locThrownTopologies.push_back("8#gammap[3#pi^{0},#eta]");
	locThrownTopologies.push_back("3#gammap[#pi^{0},#omega]");
	for(uint i=0; i<locThrownTopologies.size(); i++) {
		dHistInvariantMass_ThrownTopology[locThrownTopologies[i]] = new TH1I(Form("hInvariantMass_ThrownTopology_%d", i),Form("Invariant Mass Topology: %s", locThrownTopologies[i].Data()), 1000, 0.5, 2.0);
	}

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/		

	// cutString will be used by alot of histograms when defining names 
	// Just defining alot of stuff before hand which will be used in the labeling of the histograms 
	std::string cutString;
	std::string cutsApplied="";
	std::string cutsBase="";
	std::vector<std::string> cutVariations;
	std::string massBin;   

        cout << "INITILIZED ACTIONS" << endl;

        dHist_BeamAngle = new TH1F("BeamAngle", "Beam Angle with no cuts applied;Beam Angle (GeV)", 180,0,180);
        dHist_BeamAngle->SetYTitle("Events / Degree");
	dHist_Cuts = new TH1F("CutsPassed", "Number of times a cut has been passed", 17,0,17);
	for (int i =0; i<3; ++i){
		if (is_pi0eta){
			dHist_checkEllipseBS[i] = new TH2F(("checkEllipseBS"+std::to_string(i)+"noCutOnlyRegionSelected").c_str(), ";#pi^{0} Mass (GeV) with Events / 0.001 GeV;#eta Mass (GeV) with Events / 0.0025 GeV", atof(pi0BinRange[0].c_str()), atof(pi0BinRange[1].c_str()), atof(pi0BinRange[2].c_str()), atof(etaBinRange[0].c_str()), atof(etaBinRange[1].c_str()), atof(etaBinRange[2].c_str()));
		}
		else {
			dHist_checkEllipseBS[i] = new TH2F(("checkEllipseBS"+std::to_string(i)+"noCutOnlyRegionSelected").c_str(), ";#pi^{0} Mass (GeV) with Events / 0.001 GeV;#pi^{0} Mass (GeV) with Events / 0.001 GeV", atof(pi0BinRange[0].c_str()), atof(pi0BinRange[1].c_str()), atof(pi0BinRange[2].c_str()), atof(pi0BinRange[0].c_str()), atof(pi0BinRange[1].c_str()), atof(pi0BinRange[2].c_str()));
		}
	}

// **************************************** START INITIALZING VARIABLES TO USE WITH HISTO BUILDING ********************************************//
// **************************************** START INITIALZING VARIABLES TO USE WITH HISTO BUILDING ********************************************//

	//	else {
	//		// ********** NOT SURE IF WE SHOULD INCLUDE mMPi0P14 IN PI0PI0 REACION
	//       		++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
	//		if (cutString=="") { histCuts[id] = cutsToApply; }
	//		else { histCuts[id] = mMPi0P14; cutsApplied="mMPi0P14";}
	//       		histList[id] = {("pi0proton1D_1"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi^{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
	//       		//histCuts[id] = cutsToApply;
	//		cutsApplied=cutsBase;
	//       		++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin};
	//		if (cutString=="") { histCuts[id] = cutsToApply; }
	//		else { histCuts[id] = mMPi0P14; cutsApplied="mMPi0P14";}
	//       		histList[id] = {("pi0proton1D_2"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi^{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
	//       		//histCuts[id] = cutsToApply;
	//		cutsApplied=cutsBase;
	//       		++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locPi0Proton_Kin};
	//       		histList[id] = {("pi0pi0Pi0Proton_1"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "90", "0.", "4.5", "M(#pi^{0}#pi^{0}) (GeV) with Events / 0.025 GeV", "M(#pi^{0}Proton) (GeV) with Events / 0.05 GeV"};
	//       		histCuts[id] = cutsToApply;
	//       		++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locEtaProton_Kin};
	//       		histList[id] = {("pi0pi0Pi0Proton_2"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "90", "0.", "4.5", "M(#pi^{0}#pi^{0}) (GeV) with Events / 0.025 GeV", "M(#pi^{0}Proton) (GeV) with Events / 0.05 GeV"};
	//       		histCuts[id] = cutsToApply;
	//		
	//       		++id; histVals[id] = {5, applyAccSub, locPi0E_Kin};
	//       		histList[id] = {("pi0E_1"+cutString).c_str(), ("Cuts="+cutsApplied+"*pSelectf2;E (GeV)").c_str(), "100", "0", "10", "Events / 0.1 GeV"};
	//       		histCuts[id] = cutsToApply;
	//       		++id; histVals[id] = {6, applyAccSub, locEtaE_Kin};
	//       		histList[id] = {("pi0E_2"+cutString).c_str(), ("Cuts="+cutsApplied+"*pSelectf2;E (GeV)").c_str(), "100", "0", "10", "Events / 0.1 GeV"};
	//       		histCuts[id] = cutsToApply;
	//	}

	//       // Kinematic Hists
	//       ++id; histVals[id] = {8, applyAccSub, locPhi}; 
	//       histList[id] = {("prodPlanePSphi"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi").c_str(), "180", "-540", "540", "Events / 6 degrees"};
	//       histCuts[id] = cutsToApply;
	//       //++id; histVals[id] = {8, applyAccSub, cosTheta_decayPlane_hel};
	//       //histList[id] = {("decayPlane_cosTheta_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";cos(#theta) of decay plane").c_str(), "100","-1","1", "Events / 0.02"};
	//       //histCuts[id] = cutsToApply;
	//       //++id; histVals[id] = {8, applyAccSub, phi_decayPlane_hel};
	//       //histList[id] = {("decayPlane_phi_hel"+cutString).c_str(), ("Cuts="+cutsApplied+";#phi of decay plane").c_str(), "160","-1.6","1.6", "Events / 0.02"};
	//       //histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {8, applyAccSub, omega};
	//       histList[id] = {("vanHove_omega"+cutString).c_str(), ("Cuts="+cutsApplied+";Van Hove Omega Plot").c_str(), "120","-360","360", "60", "Events / 6 degrees"};
	//       histCuts[id] = cutsToApply;


	//
	//       // Invariant Mass Hists EBeam > 8 GeV
	////       cutsApplied=cutsBase+"*pBeamE8GeVPlus";
	////       ++id; histVals[id] = {9, applyAccSub, locPi0Proton_Kin};
	////       histList[id] = {("pi0proton1D8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#pi^{0}proton) (GeV)").c_str(), "400", "0", "4", "Events / 0.01 GeV"};
	////       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	////       ++id; histVals[id] = {10, applyAccSub, locEtaProton_Kin}; 
	////       histList[id] = {("etaproton1D8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied+";M(#etaproton) (GeV)").c_str(), "450", "0", "4.5", "Events / 0.01 GeV"};
	////       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	////       ++id; histVals[id] = {11, applyAccSub, locEtaProton_Kin,locPi0Eta_Kin}; 
	////       histList[id] = {("pi0etaProton8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "90", "0.", "4.5", "160", "0.", "4", "M(#etaProton) (GeV) with Events / 0.05 GeV", "M(#pi^{0}#eta) (GeV) with Events / 0.025 GeV"};
	////       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	////       ++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locPi0Proton_Kin};
	////       histList[id] = {("pi0etaPi0Proton8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "90", "0.", "4.5", "M(#pi^{0}#eta) (GeV) with Events / 0.025 GeV", "M(#pi^{0}Proton) (GeV) with Events / 0.05 GeV"};
	////       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	////       ++id; histVals[id] = {11, applyAccSub, locPi0Eta_Kin, locEtaProton_Kin};
	////       histList[id] = {("pi0etaEtaProton8GeVPlus"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "160", "0.", "4", "100", "0.", "5", "M(#pi^{0}#eta) (GeV) with Events / 0.025 GeV", "M(#etaProton) (GeV) with Events / 0.05 GeV"};
	////       histCuts[id] = cutsToApply*pBeamE8GeVPlus;
	////	cutsApplied=cutsBase;
	//
	//       // Charged Track Hists
	//       ++id; histVals[id] = {2, applyAccSub, locRProton}; 
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mRProton; cutsApplied="mRProton";}
	//       histList[id] = {("RadiusProton"+cutString).c_str(), ("Cuts="+cutsApplied+";Radius(proton) (cm)").c_str(), "200", "0" , "10", "Events / 0.05 cm"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {2, applyAccSub, locdzProton};
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mZMin; cutsApplied="mZMin";}
	//       histList[id] = {("dzProton"+cutString).c_str(), ("Cuts="+cutsApplied+";z(Proton) (cm)").c_str(), "160", "0" , "160", "Events / 1 cm"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {2, applyAccSub, locdEdxCDCProton};
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mdEdxCDC; cutsApplied="mdEdxCDC";}
	//       histList[id] = {("dEdxProtonCDC"+cutString).c_str(), ("Cuts="+cutsApplied+";dEdx(proton) GeV/cm").c_str(), "200", "0." , "0.00003", "Events / 1.5E-7 GeV/cm"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {2, applyAccSub, locPzProton}; 
	//       histList[id] = {("PzProton"+cutString).c_str(), ("Cuts="+cutsApplied+";Pz GeV").c_str(), "200", "0" , "5", "Events / 0.025 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locXProton,locYProton}; 
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mRProton; cutsApplied="mRProton";}
	//       histList[id] = {("XYplaneProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "100", "-4", "4", "100", "-4", "4", "x(proton) (cm) with Events / 0.08 cm", "y(proton) (cm) with Events / 0.08 cm"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {2, applyAccSub, locRProton,locdzProton};
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mRProtonZMin; cutsApplied="mRProtonZMin";}
	//       histList[id] = {("RZplaneProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "100", "0", "4", "100", "0", "100", "R(proton) (cm) wtih Events / 0.04 cm","z(proton) (cm) with Events / 1 cm"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {2, applyAccSub, locMagP3Proton,locdEdxCDCProton};
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mdEdxCDC; cutsApplied="dzRP";} // recall this is no where the final value of histCuts is set!
	//       histList[id] = {("P3dEdxCDCProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "100", "0", "4", "100", "0", "0.00003", "Momentum(proton) CDC (GeV/c) with Events / 0.04 GeV/c", "dEdx(proton) CDC (GeV/cm) with Events / 3E-7"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {2, applyAccSub, locMagP3Proton,locdEdxFDCProton};
	//       histList[id] = {("P3dEdxFDCProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "100", "0", "4", "100", "0", "0.00003", "Momentum(proton) FDC (GeV/c) with Events / 0.04 GeV/c", "dEdx(proton) FDC (GeV/cm) with Events / 3E-7"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locPzProton,locPtProton}; 
	//       histList[id] = {("PzPtProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "125", "0", "5", "125", "0", "2.5", "Pz(proton) (GeV/c) with Events / 0.04 GeV/c", "Pt(proton) (GeV/c) with Events / 0.02 GeV/c"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locPzProton,locPolarAngleProton};
	//       histList[id] = {("PzThetaProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "125", "0", "5", "250", "0", "100", "Pz(proton) (GeV/c) with Events / 0.04 GeV/c", "Theta(proton) (degrees) with Events / 0.4 degrees"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locPtProton,locPolarAngleProton};
	//       histList[id] = {("PtThetaProton"+cutString).c_str(),("Cuts="+cutsApplied).c_str() , "125", "0", "5", "250", "0", "100", "Pz(proton) (GeV/c) with Events / 0.04 GeV/c", "Theta(proton) (degrees) with Events / 0.4 degrees"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	//       histList[id] = {("PolarAngleProton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(Proton) degrees").c_str()," 200", "0", "100", "Events / 0.5 degrees"};
	//       histCuts[id] = cutsToApply;
	//	if (outputThetaRegion){
	//       		++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	//       		histList[id] = {("ThetaRegion1Proton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(degrees)").c_str(), "125", "0", "100", "Events / 0.08 degrees"};
	//       		histCuts[id] = cutsToApply*pReg1; 
	//       		++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	//       		histList[id] = {("ThetaRegion2Proton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(degrees)").c_str(), "125", "0", "100", "Events / 0.08 degrees"};
	//       		histCuts[id] = cutsToApply*pReg2; 
	//       		++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	//       		histList[id] = {("ThetaRegion3Proton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(degrees)").c_str(), "125", "0", "100", "Events / 0.08 degrees"};
	//       		histCuts[id] = cutsToApply*pReg3; 
	//       		++id; histVals[id] = {2, applyAccSub, locPolarAngleProton}; 
	//       		histList[id] = {("ThetaRegion4Proton"+cutString).c_str(), ("Cuts="+cutsApplied+";#theta(degrees)").c_str(), "125", "0", "100", "Events / 0.08 degrees"};
	//       		histCuts[id] = cutsToApply*pReg4; 
	//	}
	//
	//
	//       // Neutral Track Hists
	//       ++id; histVals[id] = {13, applyAccSub, photonEnergies[0]};
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mPhotonE; cutsApplied="mPhotonE";}
	//       histList[id] = {("PhotonShowerE1"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy(#gamma) GeV").c_str(), "200", "0" , "10", "Events / 0.05 GeV"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {13, applyAccSub, photonThetas[0]}; 
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mPhotonTheta; cutsApplied="mPhotonTheta";}
	//       histList[id] = {("PhotonShowerTheta1"+cutString).c_str(), ("Cuts="+cutsApplied+";PolarAngle(#gamma) degrees").c_str(), "300", "0" , "150", "Events / 0.5 degrees"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {13, applyAccSub, photonEnergies[0], photonThetas[0]};
	//       histList[id] = {("thetaEPhotonShower"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "10", "400", "0", "150", "E(#gamma) (GeV) wtih Events / 0.05", "#theta(#gamma) (radians) with Events / 0.25"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, photonXs_Shower[0], photonYs_Shower[0]}; 
	//       histList[id] = {("XYPhotonShower"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "260", "-130", "130", "260", "-130", "130", "x(cm) with Events / 1 cm", "y(cm) with Events / 1 cm"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, photonXs_Shower[0], photonYs_Shower[0]}; 
	//       histList[id] = {("XYPhotonShowerBCAL"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "260", "-130", "130", "260", "-130", "130", "x(cm) with Events / 1 cm", "y(cm) with Events / 1 cm"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, photonPhis[0], photonThetas[0]};
	//       histList[id] = {("ThetaPhiPhotonShower"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "360", "-180", "180", "130", "0", "130", "#phi(degrees) with Events / 1 cm", "#theta(degrees) with Events / 1 cm"};
	//       histCuts[id] = cutsToApply;
	//       // These below are tracked differently than the charged tracks and neutral tracks so they have their own unique identifiers.
	//
	//       // Timing
	//       if(outputTimingHists){
	//       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton,RFtimeProton}; 
	//       		histList[id] = {("timeProtonRFvsTheta"+cutString).c_str(),("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees","T_{Proton} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply;
	//       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	//       		histList[id] = {("timeProtonRFvsP3"+cutString).c_str(),("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c", "T_{Proton} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply;
	//       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton,RFtimeProton}; 
	//       		histList[id] = {("timeBCALRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees","T_{BCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	//       		histList[id] = {("timeBCALRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120","-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{BCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
	//       		histList[id] = {("timeTOFRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees", "T_{TOF} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	//       		histList[id] = {("timeTOFRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{TOF} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
	//       		histList[id] = {("timeFCALRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200"," 0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 cm", "T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	//       		histList[id] = {("timeFCALRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c", "T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
	//       		histList[id] = {("timeSTARTRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees", "T_{START} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	//       		histList[id] = {("timeSTARTRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6", "6", "Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{START} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, locPolarAngleProton, RFtimeProton};
	//       		histList[id] = {("timeSYSNULLRFvsTheta"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "100", "120", "-6", "6", "#theta(Proton) degrees with Events / 0.5 degrees","T_{SYSNULL} - T_{RF} (ns) with Events / 0.1 ns "};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, locMagP3Proton,RFtimeProton}; 
	//       		histList[id] = {("timeSYSNULLRFvsP3"+cutString).c_str(), ("Cuts="+cutsApplied).c_str(), "200", "0", "5", "120", "-6" , "6","Proton track momentum (GeV/c) with Events / 0.025 GeV/c","T_{SYSNULL} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 



	//       		++id; histVals[id] = {8, applyAccSub, photonThetas[0], photonDeltaTs[0]}; 
	//       		histList[id] = {("timePhotonBCALFCALRFvsTheta"+cutString).c_str(),("Cuts="+cutsApplied).c_str(),"300","0","120","120","-6","6", "#theta(Proton) degrees with Events / 0.4 degrees","T_{BCAL/FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		++id; histVals[id] = {8, applyAccSub, photonDeltaTs[0]}; 
	//       		histList[id] = {("timePhotonBCALFCALRF"+cutString).c_str(),("Cuts="+cutsApplied+";T_{BCAL/FCAL} - T_{RF} (ns)").c_str(),"120","-6","6", "Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		histCuts[id] = cutsToApply;// recall this is not where we actually fill the histCuts, its just for checking sizes of things at this point... 
	//       		++id; histVals[id] = {8, applyAccSub, photonThetas[0], photonDeltaTs[0]}; 
	//       		histList[id] = {("timePhotonFCALRFvsTheta"+cutString).c_str(),("Cuts="+cutsApplied).c_str(),"60","0","12","120","-6","6", "#theta(Proton) degrees with Events / 0.2 degrees","T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, photonEnergies[0], photonDeltaTs[0]};
	//       		histList[id] = {("timePhotonFCALRFvsEnergy"+cutString).c_str(),("Cuts="+cutsApplied).c_str(),"200","0","12","120","-6","6", "Photon momentum (GeV/c) with Events / 0.06 GeV/c", "T_{FCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply*pPhotonInFCAL[0]; 
	//       		++id; histVals[id] = {8, applyAccSub, photonDeltaTs[0]};                                                                                    
	//       		histList[id] = {("timePhotonFCALRF"+cutString).c_str(),("Cuts="+cutsApplied+";T_{FCAL} - T_{RF} (ns)").c_str(),"120","-6","6", "Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, photonEnergies[0], photonDeltaTs[0]};
	//       		histList[id] = {("timePhotonBCALRFvsEnergy"+cutString).c_str(),("Cuts="+cutsApplied).c_str(),"200","0","6","120","-6","6", "Photon momentum (GeV/c) with Events / 0.03 GeV/c","T_{BCAL} - T_{RF} (ns) with Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply; 
	//       		++id; histVals[id] = {8, applyAccSub, photonDeltaTs[0]}; 
	//       		histList[id] = {("timePhotonBCALRF"+cutString).c_str(),("Cuts="+cutsApplied+";T_{BCAL} - T_{RF} (ns)").c_str(),"120","-6","6", "Events / 0.1 ns"};
	//       		histCuts[id] = cutsToApply;                                                                                                  
	//	}
	//

	//       // Non track related hists
	//       ++id; histVals[id] = {8, applyAccSub, locCLKinFit}; 
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mChiSq; cutsApplied="mChiSq";} 
	//       histList[id] = {("P4CLKinFit"+cutString).c_str(), ("Cuts="+cutsApplied+";Confidence Level").c_str() , "200", "0", ".30", "Entries / 0.0015"};
	//       cutsApplied=cutsBase;
	//       ++id; histVals[id] = {8, applyAccSub, locUnusedEnergy};
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mUE; cutsApplied="mUE";} 
	//       histList[id] = {("UnusedEnergy"+cutString).c_str(), ("Cuts="+cutsApplied+";Unused Energy (GeV)").c_str(), "100", "0", "1", "Entries / 0.01 GeV"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {8, applyAccSub, locChiSqKinFit}; 
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mChiSq; cutsApplied="mChiSq";} 
	//       histList[id] = {("P4ChiSqKinFit"+cutString).c_str(), ("Cuts="+cutsApplied+";Chi Squared").c_str(), "150", "0", "300", "Entries / 1"};
	//	cutsApplied=cutsBase;
	//       ++id; histVals[id] = {8, applyAccSub, locDOFKinFit}; 
	//       if (cutString=="") { histCuts[id] = cutsToApply; }
	//       else { histCuts[id] = mChiSq; cutsApplied="mChiSq";} 
	//       histList[id] = {("P4DOFKinFit"+cutString).c_str(), ("Cuts="+cutsApplied+";DOF").c_str(), "30", "0", "30", "Entries / 1"};
	//	cutsApplied=cutsBase;

	//	//Shower shape variables
	//	// neutral showers
	//       ++id; histVals[id] = {13, applyAccSub, E1E9_FCAL[0]}; 
	//       histList[id] = {("E1E9_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";E1E9_FCAL").c_str(), "100", "-1", "2", "Entries / 0.03"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, E9E25_FCAL[0]}; 
	//       histList[id] = {("E9E25_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";E9E25_FCAL").c_str(), "100", "-1", "2", "Entries / 0.03"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, SumU_FCAL[0]}; 
	//       histList[id] = {("SumU_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SumU_FCAL").c_str(), "200", "-1", "9", "Entries / 0.5"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, SumV_FCAL[0]}; 
	//       histList[id] = {("SumV_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SumV_FCAL").c_str(), "200", "-1", "9", "Entries / 0.5"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, DOCA_FCAL[0]}; 
	//       histList[id] = {("DOCA_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";DOCA_FCAL").c_str(), "160", "0", "16", "Entries / 0.1 cm"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, showerQuality_FCAL[0]}; 
	//       histList[id] = {("showerQuality_FCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";showerQuality_FCAL").c_str(), "50", "0", "1", "Entries / 0.02"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, Energy_BCALPreshower[0]}; 
	//       histList[id] = {("Energy_BCALPreshower"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy_BCALPreshower").c_str(), "200", "-1", "9", "Entries / 0.5 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, Energy_BCAL[0]}; 
	//       histList[id] = {("Energy_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, SigLong_BCAL[0]}; 
	//       histList[id] = {("sigLong_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SigLong_BCAL").c_str(), "200", "0", "10", "Entries / 0.05 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, SigTrans_BCAL[0]}; 
	//       histList[id] = {("sigTrans_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTrans_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, SigTheta_BCAL[0]}; 
	//       histList[id] = {("sigTheta_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTheta_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, DeltaPhi_BCAL[0]}; 
	//       histList[id] = {("DeltaPhi_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";DeltaPhi_BCAL").c_str(), "100", "-5", "5", "Entries / 0.1 degree"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {13, applyAccSub, DeltaZ_BCAL[0]}; 
	//       histList[id] = {("DeltaZ_BCAL"+cutString).c_str(), ("Cuts="+cutsApplied+";DeltaZ_BCAL").c_str(), "200", "-100", "100", "Entries / 1 cm"};
	//       histCuts[id] = cutsToApply;
	//	
	//	// charged showers
	//       ++id; histVals[id] = {2, applyAccSub, locE1E9_FCAL_proton}; 
	//       histList[id] = {("E1E9_FCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";E1E9_FCAL").c_str(), "100", "-1", "2", "Entries / 0.03"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locE9E25_FCAL_proton}; 
	//       histList[id] = {("E9E25_FCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";E9E25_FCAL").c_str(), "100", "-1", "2", "Entries / 0.03"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locSumU_FCAL_proton}; 
	//       histList[id] = {("SumU_FCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";SumU_FCAL").c_str(), "200", "-1", "9", "Entries / 0.05"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locSumV_FCAL_proton}; 
	//       histList[id] = {("SumV_FCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";SumV_FCAL").c_str(), "200", "-1", "9", "Entries / 0.05"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locEnergy_BCALPreshower_proton}; 
	//       histList[id] = {("Energy_BCALPreshowerproton"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy_BCALPreshower").c_str(), "100", "-1", "9", "Entries / 0.05 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locEnergy_BCAL_proton}; 
	//       histList[id] = {("Energy_BCALproton"+cutString).c_str(), ("Cuts="+cutsApplied+";Energy_BCAL").c_str(), "200", "-1", "9", "Entries / 0.05 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locSigLong_BCAL_proton}; 
	//       histList[id] = {("sigLong_BCAL_proton"+cutString).c_str(), ("Cuts="+cutsApplied+";SigLong_BCAL").c_str(), "200", "0", "10", "Entries / 0.05 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locSigTrans_BCAL_proton}; 
	//       histList[id] = {("sigTrans_BCAL_proton"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTrans_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	//       histCuts[id] = cutsToApply;
	//       ++id; histVals[id] = {2, applyAccSub, locSigTheta_BCAL_proton}; 
	//       histList[id] = {("sigTheta_BCAL_proton"+cutString).c_str(), ("Cuts="+cutsApplied+";SigTheta_BCAL").c_str(), "80", "0", "0.2", "Entries / 0.0025 GeV"};
	//       histCuts[id] = cutsToApply;




	//cutVariations = {"", "Cut", "mUE", "mUEChiSq100"};
	//for (int j=0; j<4; ++j){
	//       if (j==0) { cutsToApply = 1; cutsApplied="None"; applyAccSub=noAccSub;}
	//       if (j==1) { cutsToApply = allGeneralCutsPassed; cutsApplied="GeneralCuts"; applyAccSub=weight;}
	//       if (j==2) { cutsToApply = mUE; cutsApplied="mUE"; applyAccSub=weight;}
	//       if (j==3) { cutsToApply = mUEChiSq*chiSq100; cutsApplied="mUEChiSq+ChiSq100"; applyAccSub=weight;}
	//       cutString = cutVariations[j];
	//       ++id; histVals[id] = {8, applyAccSub, locNumExtraNeutralShowers}; 
	//       histList[id] = {("numExtraNeutralShowers"+cutString).c_str(),("Cuts="+cutsApplied+";Number of extra showers").c_str(),"10","0","10", "Events / 1"};
	//       histCuts[id] = cutsToApply;                                                                                                  
	//}




        string name;
        histDef_1D histdef;
        histDef_2D histdef2d;
	// ********************************** Kinematics related *********************************************
        histdef.clear();
        name="P4ChiSqKinFit_mChiSq";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mChiSq;Chi Squared; Entries / 2", 150, 0, 300);
        histdef.name = name; histdef.cut=&mChiSq; histdef.weights = &weightAS;
        histdef.values.push_back( &locChiSqKinFit );
        group_1234BP.insert(histdef); 

        histdef.clear();
        name="UnusedEnergy_noCut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=noCut;Chi Squared;Entries / 0.01 GeV", 100, 0, 1);
        histdef.name = name; histdef.cut=&noCut; histdef.weights = &weightAS;
        histdef.values.push_back( &locUnusedEnergy );
        group_1234BP.insert(histdef); 

	// ********************************** DETECTOR SPECIFIC RELATED PLOTS FOR PI0/ETA *****************************************
        histdef.clear();
        name="pi0DetectedIn";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tAll;0-FCAL 1-BCAL 2-SPLIT;", 3, 0, 3);
        histdef.name = name; histdef.cut=&mEllipse_pre_tAll; histdef.weights = &noWeight;
        histdef.values.push_back( &pi0DetectedIn );
        group_12B.insert(histdef); 
        histdef.clear();
        name="etaDetectedIn";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tAll;0-FCAL 1-BCAL 2-SPLIT;", 3, 0, 3);
        histdef.name = name; histdef.cut=&mEllipse_pre_tAll; histdef.weights = &noWeight;
        histdef.values.push_back( &etaDetectedIn );
        group_12B.insert(histdef); 

	// ********************************** DECK REALTED PLOTS *****************************************
        histdef.clear();
        name="mandelstam_t";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT;-t momentum transfer of #pi^{0}+#eta;Events / 0.06 GeV", 100, 0, 6);
        histdef.name = name; histdef.cut=&mMandelstamT; histdef.weights = &weightAS;
        histdef.values.push_back( &mandelstam_abst );
        group_1234B.insert(histdef); 

	// Will leave bin to contain all the bad regions
        histdef.clear();
        name="tetaVsMpi0eta_recCounts_tLT05";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tLT05;efficiency", numHists+1, -1, numHists);
        histdef.name = name; histdef.cut=&mEllipse_pre_tLT05; histdef.weights = &weightAS;
        histdef.values.push_back( &teta_recCounts );
	//    chose this uniqueness set since I use the 2D distribution to select out the regions. So technically I use this set
        group_34B_1234B.insert(histdef); 

        histdef.clear();
        name="tetaVsMpi0eta_recCounts_tGT05LT1";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tGT05LT1;efficiency", numHists+1, -1, numHists);
        histdef.name = name; histdef.cut=&mEllipse_pre_tGT05LT1; histdef.weights = &weightAS;
        histdef.values.push_back( &teta_recCounts );
        group_34B_1234B.insert(histdef); 

        histdef.clear();
        name="tetaVsMpi0eta_recCounts_tGT1";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tGT1;efficiency", numHists+1, -1, numHists);
        histdef.name = name; histdef.cut=&mEllipse_pre_tGT1; histdef.weights = &weightAS;
        histdef.values.push_back( &teta_recCounts );
        group_34B_1234B.insert(histdef); 

        histdef.clear();
        name="tetaVsMpi0eta_recCounts_tAll";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tAll;efficiency", numHists+1, -1, numHists);
        histdef.name = name; histdef.cut=&mEllipse_pre_tAll; histdef.weights = &weightAS;
        histdef.values.push_back( &teta_recCounts );
        group_34B_1234B.insert(histdef); 

        histdef.clear();
        name="tpi0VsMpi0eta_recCounts_tLT05";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tLT05;efficiency", numHists+1, -1, numHists);
        histdef.name = name; histdef.cut=&mEllipse_pre_tLT05; histdef.weights = &weightAS;
        histdef.values.push_back( &tpi0_recCounts );
        group_34B_1234B.insert(histdef); 

        histdef.clear();
        name="tpi0VsMpi0eta_recCounts_tGT05LT1";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tGT05LT1;efficiency", numHists+1, -1, numHists);
        histdef.name = name; histdef.cut=&mEllipse_pre_tGT05LT1; histdef.weights = &weightAS;
        histdef.values.push_back( &tpi0_recCounts );
        group_34B_1234B.insert(histdef); 

        histdef.clear();
        name="tpi0VsMpi0eta_recCounts_tGT1";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tGT1;efficiency", numHists+1, -1, numHists);
        histdef.name = name; histdef.cut=&mEllipse_pre_tGT1; histdef.weights = &weightAS;
        histdef.values.push_back( &tpi0_recCounts );
        group_34B_1234B.insert(histdef); 

        histdef.clear();
        name="tpi0VsMpi0eta_recCounts_tAll";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tAll;efficiency", numHists+1, -1, numHists);
        histdef.name = name; histdef.cut=&mEllipse_pre_tAll; histdef.weights = &weightAS;
        histdef.values.push_back( &tpi0_recCounts );
        group_34B_1234B.insert(histdef); 

	//for ( int iMass=0; iMass < num_massBins; ++iMass){
        //	histdef.clear();
        //	name="tetaMassBinned"+to_string(iMass);
        //	histdef.hist = new TH1F(name.c_str(), "Cuts=passMassBin_tetaIntegrated;t_{#eta} (GeV^2)", num_tBins,tMin,tMax);
        //	histdef.name = name; histdef.cut=&passMassBin_tetaIntegrated[iMass]; histdef.weights = &weightAS;
        //	histdef.values.push_back( &mandelstam_teta );
        //	group_34B_1234B.insert(histdef); 

        //	histdef.clear();
        //	name="tpi0MassBinned"+to_string(iMass);
        //	histdef.hist = new TH1F(name.c_str(), "Cuts=passMassBin_tpi0Integrated;t_{#eta} (GeV^2)", num_tBins,tMin,tMax);
        //	histdef.name = name; histdef.cut=&passMassBin_tpi0Integrated[iMass]; histdef.weights = &weightAS;
        //	histdef.values.push_back( &mandelstam_tpi0 );
        //	group_34B_1234B.insert(histdef); 
	//}

	// *************************** PHOTON VECTOR QUANTITIES ****************************

        histdef.clear();
        cutsApplied="true";
        histdef.hist = new TH1F("X_Kin","Cuts=GeneralCuts;position (cm);Events / 1 cm", 160, 0, 160);
        histdef.name = "X_Kin";
        histdef.cut = &allGeneralCutsPassed;
        histdef.weights = &noWeight;
        histdef.values.push_back( &(photonXs_Kin[0]) );
        histdef.values.push_back( &(photonXs_Kin[1]) );
        histdef.values.push_back( &(photonXs_Kin[2]) );
        histdef.values.push_back( &(photonXs_Kin[3]) );
        group_PhNB.insert(histdef);
       
        histdef.clear();
        histdef.hist = new TH1F("Y_Kin","Cuts=GeneralCuts;position (cm);Events / 1 cm", 160, 0, 160);
        histdef.name = "Y_Kin";
        histdef.cut = &allGeneralCutsPassed;
        histdef.weights = &noWeight;
        histdef.values.push_back( &(photonYs_Kin[0]) );
        histdef.values.push_back( &(photonYs_Kin[1]) );
        histdef.values.push_back( &(photonYs_Kin[2]) );
        histdef.values.push_back( &(photonYs_Kin[3]) );
        group_PhNB.insert(histdef);

        histdef.clear();
        cutsApplied="true";
        histdef.hist = new TH1F("Z_Kin","Cuts=GeneralCuts;position (cm);Events / 1 cm", 160, 0, 160);
        histdef.name = "Z_Kin";
        histdef.cut = &allGeneralCutsPassed;
        histdef.weights = &noWeight;
        histdef.values.push_back( &(photonZs_Kin[0]) );
        histdef.values.push_back( &(photonZs_Kin[1]) );
        histdef.values.push_back( &(photonZs_Kin[2]) );
        histdef.values.push_back( &(photonZs_Kin[3]) );
        group_PhNB.insert(histdef);

        histdef.clear();
        cutsApplied="true";
        histdef.hist = new TH1F("T_Kin","Cuts=GeneralCuts;position (cm);Events / 1 cm", 160, 0, 160);
        histdef.name = "T_Kin";
        histdef.cut = &allGeneralCutsPassed;
        histdef.weights = &noWeight;
        //histdef.values.insert(histdef.values.end(), &std::begin(photonXs_Kin), &std::end(photonXs_Kin)); 
        histdef.values.push_back( &(photonTs_Kin[0]) );
        histdef.values.push_back( &(photonTs_Kin[1]) );
        histdef.values.push_back( &(photonTs_Kin[2]) );
        histdef.values.push_back( &(photonTs_Kin[3]) );
        group_PhNB.insert(histdef);

        histdef.clear();
        histdef.hist = new TH1F("X_Shower","Cuts=GeneralCuts;position (cm);Events / 1 cm", 160, 0, 160);
        histdef.name = "X_Shower";
        histdef.cut = &allGeneralCutsPassed;
        histdef.weights = &noWeight;
        histdef.values.push_back( &(photonXs_Shower[0]) );
        histdef.values.push_back( &(photonXs_Shower[1]) );
        histdef.values.push_back( &(photonXs_Shower[2]) );
        histdef.values.push_back( &(photonXs_Shower[3]) );
        group_PhNB.insert(histdef);

        histdef.clear();
        histdef.hist = new TH1F("Y_Shower","Cuts=GeneralCuts;position (cm);Events / 1 cm", 160, 0, 160);
        histdef.name = "Y_Shower";
        histdef.cut = &allGeneralCutsPassed;
        histdef.weights = &noWeight;
        histdef.values.push_back( &(photonYs_Shower[0]) );
        histdef.values.push_back( &(photonYs_Shower[1]) );
        histdef.values.push_back( &(photonYs_Shower[2]) );
        histdef.values.push_back( &(photonYs_Shower[3]) );
        group_PhNB.insert(histdef);

        histdef.clear();
        histdef.hist = new TH1F("Z_Shower","Cuts=GeneralCuts;position (cm);Events / 1 cm", 160, 0, 160);
        histdef.name = "Z_Shower";
        histdef.cut = &allGeneralCutsPassed;
        histdef.weights = &noWeight;
        histdef.values.push_back( &(photonZs_Shower[0]) );
        histdef.values.push_back( &(photonZs_Shower[1]) );
        histdef.values.push_back( &(photonZs_Shower[2]) );
        histdef.values.push_back( &(photonZs_Shower[3]) );
        group_PhNB.insert(histdef);
 
        histdef.clear();
        histdef.hist = new TH1F("T_Shower","Cuts=GeneralCuts;position (cm);Events / 1 cm", 160, 0, 160);
        histdef.name = "T_Shower";
        histdef.cut = &allGeneralCutsPassed;
        histdef.weights = &noWeight;
        histdef.values.push_back( &(photonTs_Shower[0]) );
        histdef.values.push_back( &(photonTs_Shower[1]) );
        histdef.values.push_back( &(photonTs_Shower[2]) );
        histdef.values.push_back( &(photonTs_Shower[3]) );
        group_PhNB.insert(histdef);

        histdef.clear();
        histdef.hist = new TH1F("pi0gamma","Cuts=mEllipse_pre; Mass GeV;Events / 0.01 cm", 350, 0, 3.5);
        histdef.name = "pi0gamma";
        histdef.cut = &mEllipse_pre;
        histdef.weights = &weightAS;
        histdef.values.push_back( &( massGammaPi0[0]) );
        histdef.values.push_back( &( massGammaPi0[1]) );
        group_12PhNB.insert(histdef);

        histdef.clear();
        histdef.hist = new TH1F("etagamma","Cuts=mEllipse_pre; Mass GeV;Events / 0.01 cm", 350, 0, 3.5);
        histdef.name = "etagamma";
        histdef.cut = &mEllipse_pre;
        histdef.weights = &weightAS;
        histdef.values.push_back( &( massGammaEta[0]) );
        histdef.values.push_back( &( massGammaEta[1]) );
        group_34PhNB.insert(histdef);

        // ************************** PROTON RELATED HISTS ***********************
        histdef.clear();
        name="RadiusProton_mRProton";
        histdef.hist = new TH1F(name.c_str(),"Cuts=mRProton;Radius(proton) (cm);Events / 0.05 cm", 200, 0 , 10);
        histdef.name = name; histdef.cut = &mRProton; histdef.weights = &weightAS;
        histdef.values.push_back( &locRProton ); 
        group_PB.insert(histdef);

        histdef.clear();
        name="dzProton_ZMin";
        histdef.hist = new TH1F(name.c_str(), "Cuts=Zmin;z(Proton) (cm);Events / 1 cm", 160, 0 , 160);
        histdef.name = name; histdef.cut = &mZMin; histdef.weights = &weightAS;
        histdef.values.push_back( &locdzProton );
        group_PB.insert(histdef);

        histdef.clear();
        name="dEdxProtonCDC_mdEdxCDC";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mdEdxCDC;dEdx(proton) GeV/cm;Events / 1.5E-7 GeV/cm", 200, 0. , 0.00003);
        histdef.name = name; histdef.cut=&mdEdxCDC; histdef.weights = &weightAS;
        histdef.values.push_back( &locdEdxCDCProton );
        group_PB.insert(histdef);

        histdef.clear();
        name = "PzProton_Cut";
        histdef.hist = new TH1F(name.c_str(),  "Cuts=GeneralCuts;Pz GeV;Events / 0.025 GeV", 200, 0, 5 );
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &locPzProton );
        group_PB.insert(histdef); 

        histdef.clear();
        name="PolarAngleProton_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts;#theta(Proton) degrees;Events / 0.5 degrees", 200, 0, 100);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed ; histdef.weights = &weightAS;
        histdef.values.push_back( &locPolarAngleProton );
        group_PB.insert(histdef); 

        histdef2d.clear();
        name = "XYplaneProton_mRProton";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mRProton;x(proton) (cm) with Events / 0.08 cm;y(proton) (cm) with Events / 0.08 cm", 100,-4,4,100,-4,4);
        histdef2d.name = name; histdef2d.cut=&mRProton; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locXProton );
        histdef2d.valuesY.push_back( &locYProton );
        group_PB.insert_2D(histdef2d); 

        histdef2d.clear();
        name = "RZplaneProton_mRProtonZmin";
        histdef2d.hist = new TH2F(name.c_str(),"Cuts=mRProtonZmin;R(proton) (cm) wtih Events / 0.04 cm;z(proton) (cm) with Events / 1 cm",100,0,4,100,0,100);
        histdef2d.name = name; histdef2d.cut=&mRProtonZMin; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locRProton );
        histdef2d.valuesY.push_back( &locdzProton );
        group_PB.insert_2D(histdef2d); 
        
        histdef.clear();
        name = "P3Proton_dzR";
        histdef.hist = new TH1F(name.c_str(), "Cuts=dzR;Momentum(proton) CDC (GeV/c) with Events / 0.04 GeV/c",100,0,4);
        histdef.name = name; histdef.cut=&dzR; histdef.weights = &weightAS;
        histdef.values.push_back( &locMagP3Proton );
        group_PB.insert(histdef); 

        histdef2d.clear();
        name = "P3dEdxCDCProton_dzRP";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=dzRP;Momentum(proton) CDC (GeV/c) with Events / 0.04 GeV/c;dEdx(proton) CDC (GeV/cm) with Events / 3E-7",100,0,4,100,0,0.00003);
        histdef2d.name = name; histdef2d.cut=&dzRP; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locMagP3Proton );
        histdef2d.valuesY.push_back( &locdEdxCDCProton );
        group_PB.insert_2D(histdef2d); 

        histdef2d.clear();
        name = "P3dEdxFDCProton_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=GeneralCuts;Momentum(proton) FDC (GeV/c) with Events / 0.04 GeV/c;dEdx(proton) FDC (GeV/cm) with Events / 3E-7",100,0,4,100,0,0.00003);
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locMagP3Proton );
        histdef2d.valuesY.push_back( &locdEdxFDCProton );
        group_PB.insert_2D(histdef2d); 
        
        histdef2d.clear();
        name = "PzPtProton_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=GeneralCuts;Pz(proton) (GeV/c) with Events / 0.04 GeV/c;Pt(proton) (GeV/c) with Events / 0.04 GeV/c",125,0,5,125,0,5);
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPzProton );
        histdef2d.valuesY.push_back( &locPtProton );
        group_PB.insert_2D(histdef2d);

        histdef2d.clear();
        name = "PzThetaProton_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=GeneralCuts;Pz(proton) (GeV/c) with Events / 0.04 GeV/c;Theta(proton) (degrees) with Events / 0.4 degrees",125,0,5,250,0,100);
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPzProton );
        histdef2d.valuesY.push_back( &locPolarAngleProton );
        group_PB.insert_2D(histdef2d); 

        histdef2d.clear();
        name="PtThetaProton_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=GeneralCuts;Pt(proton) (GeV/c) with Events / 0.04 GeV/c;Theta(proton) (degrees) with Events / 0.4 degrees",125,0,5,250,0,100);
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPtProton );
        histdef2d.valuesY.push_back( &locPolarAngleProton );
        group_PB.insert_2D(histdef2d); 

        // ********************** ANGULAR PLOTS ********************
        
        histdef.clear();
        name="eta_cosTheta_hel_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts;cos(#theta) of #pi^{0};Events / 0.02", 100,-1,1);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &cosTheta_eta_hel );
        group_34B.insert(histdef); 
        
        histdef.clear();
        name = "pi0eta_cosTheta_hel_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts;cos(#theta) of #pi^{0}+#eta;Events / 0.02",100,-1,1);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &cosTheta_pi0eta_hel );
        group_1234B.insert(histdef); 

        histdef.clear();
        name = "eta_phi_hel_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts;#phi of #eta;Events / 4 degrees",90,-180,180);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &phi_eta_hel);
        group_34B.insert(histdef); 

        histdef.clear();
        name = "pi0eta_phi_hel_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cut=GeneralCuts;#phi of #pi^{0}+#eta;Events / 4 degrees",90,-180,180);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &phi_pi0eta_hel );
        group_1234B.insert(histdef); 

        histdef.clear();
        name = "eta_cosTheta_GJ_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cut=GeneralCuts;cos(#theta) of #eta;Events / 0.02",100,-1,1);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &cosTheta_eta_GJ );
        group_34B.insert(histdef); 

        histdef.clear();
        name="eta_phi_GJ_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cut=GeneralCuts;#phi of #eta;Events / 4 degrees", 90,-180,180);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &phi_eta_GJ );
        group_34B.insert(histdef); 

        histdef2d.clear();
        name = "eta_cosTheta_GJvsM_mMandelstamT";
        histdef2d.hist = new TH2F(name.c_str(), "Cut=mMandelstamT;M(#pi^{0}#eta)Events / 0.01 GeV;Cos(#theta) of #eta Events / 0.2",350,0,3.5,100,-1,1);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &cosTheta_eta_GJ );
        group_34B_1234B.insert_2D(histdef2d); 

        histdef2d.clear();
        name = "eta_cosTheta_GJvsMpi0p_mMandelstamT_mdelta";
        histdef2d.hist = new TH2F(name.c_str(), "Cut=mMandelstamT_mdelta;M(#pi^{0} p)Events / 0.01 GeV;Cos(#theta) of #eta Events / 0.2",350,0,3.5,100,-1,1);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT_mdelta; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Proton_Kin );
        histdef2d.valuesY.push_back( &cosTheta_eta_GJ );
        group_34B_1234B.insert_2D(histdef2d); 

        histdef2d.clear();
        name = "eta_cosTheta_GJvsMetap_mMandelstamT_mdelta";
        histdef2d.hist = new TH2F(name.c_str(), "Cut=mMandelstamT_mdelta;M(#eta p)Events / 0.01 GeV;Cos(#theta) of #eta Events / 0.2",400,0,4,100,-1,1);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT_mdelta; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locEtaProton_Kin );
        histdef2d.valuesY.push_back( &cosTheta_eta_GJ );
        group_34B_1234B.insert_2D(histdef2d); 

        name = "eta_phi_GJvsM_mMandelstamT";
        histdef2d.hist = new TH2F(name.c_str(), "Cut=mMandelstamT AS;M(#pi^{0}#eta)Events / 0.01 GeV;#phi of #eta Events / 4 degrees",350,0,3.5,100,-180,180);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &phi_eta_GJ );
        group_34B_1234B.insert_2D(histdef2d); 

        histdef2d.clear();
        name = "eta_cosTheta_GJvsM_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cut=GeneralCuts;M(#pi^{0}#eta)Events / 0.01 GeV;Cos(#theta) of #eta Events / 0.2",350,0,3.5,100,-1,1);
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &cosTheta_eta_GJ );
        group_34B_1234B.insert_2D(histdef2d); 
        histdef2d.clear();

        name = "eta_phi_GJvsM_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cut=GeneralCuts;M(#pi^{0}#eta)Events / 0.01 GeV;#phi of #eta Events / 4 degrees",350,0,3.5,90,-180,180);
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &phi_eta_GJ );
        group_34B_1234B.insert_2D(histdef2d); 

        histdef.clear();
        name = "pi0eta_cosTheta_CM_Cut";
        histdef.hist = new TH1F(name.c_str(),"Cut=GeneralCuts;cos(#theta) of #pi^{0}+#eta;Events / 0.02", 100,-1,1 );
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &cosTheta_pi0eta_CM );
        group_1234B.insert(histdef); 
        
        histdef.clear();
        name = "pi0eta_phi_CM_Cut";
        histdef.hist = new TH1F(name.c_str(),"Cut=GeneralCuts;#phi of #pi^{0}+#eta;Events / 4 degrees", 90,-180,180 );
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &phi_pi0eta_CM );
        group_PB.insert(histdef); 

        // *********************** MASS SHIFT PLOTS *************************
        histdef.clear();
        name="pi0Mass_chargedVertex";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts;M(#pi^{0}) GeV; Events / 0.001 GeV", 200, 0.05,0.25);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Mass_charged );
        group_12B.insert(histdef); 
        
        histdef.clear();
        name="etaMass_chargedVertex";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts;M(#eta) GeV; Events / 0.002 GeV", 300, 0.25, 0.85);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass_charged );
        group_34B.insert(histdef); 

        histdef2d.clear();
        name="pi0eta_chargedVertex";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=GeneralCuts;M(#pi^{0}) GeV Events / 0.001 GeV;#eta Mass (GeV) with Events / 0.002 GeV", 200,0.05,0.25,300,0.25,0.85 );
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Mass_charged );
        histdef2d.valuesY.push_back( &locEtaMass_charged );
        group_12B_34B.insert_2D(histdef2d); 
        
        histdef.clear();
        name="pi0Mass_targetVertex";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts;M(#pi^{0}) GeV; Events / 0.001 GeV", 200, 0.05,0.25);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Mass_target );
        group_12B.insert(histdef); 
        
        histdef.clear();
        name="etaMass_targetVertex";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts;M(#eta) GeV; Events / 0.002 GeV", 300, 0.25, 0.85);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass_target );
        group_34B.insert(histdef); 

        histdef2d.clear();
        name="pi0eta_targetVertex";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=GeneralCuts;M(#pi^{0}) GeV Events / 0.001 GeV;#eta Mass (GeV) with Events / 0.002 GeV", 200,0.05,0.25,300,0.25,0.85 );
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Mass_target );
        histdef2d.valuesY.push_back( &locEtaMass_target );
        group_12B_34B.insert_2D(histdef2d); 

        // ********************* MASS PLOTS IN DIFFERENT PARTS OF DETECTOR **********************
        histdef.clear();
        name="pi0Mass_Kin_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre;M(#pi^{0}) (GeV);Events / 0.001 GeV",200,0.05,0.25);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Mass_Kin );
        group_12B.insert(histdef); 

        histdef.clear();
        name="pi0Mass_Meas_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre;M(#pi^{0}) (GeV);Events / 0.001 GeV",200,0.05,0.25);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Mass );
        group_12B.insert(histdef); 
        histdef.clear();

        histdef.clear();
        name="pi0MassFCAL_Kin_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre*InFCAL;M(#pi^{0}) (GeV);Events / 0.001 GeV",200,0.05,0.25);
        histdef.name = name; histdef.cut=&mEllipse_pre_pi0FCAL; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Mass_Kin );
        group_12B.insert(histdef); 
        
        histdef.clear();
        name="pi0MassBCAL_Kin_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre*InBCAL;M(#pi^{0}) (GeV);Events / 0.001 GeV",200,0.05,0.25);
        histdef.name = name; histdef.cut=&mEllipse_pre_pi0BCAL; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Mass_Kin );
        group_12B.insert(histdef); 

        histdef.clear();
        name="pi0MassSPLIT_Kin_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre*SPLIT;M(#pi^{0}) (GeV);Events / 0.001 GeV",200,0.05,0.25);
        histdef.name = name; histdef.cut=&mEllipse_pre_pi0SPLIT; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Mass_Kin );
        group_12B.insert(histdef); 

        histdef.clear();
        name="etaMass_Kin_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre;M(#eta) (GeV);Events / 0.002 GeV",300,0.25,0.85);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass_Kin );
        group_34B.insert(histdef); 

        histdef.clear();
        name="etaMass_Meas_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre;M(#eta) (GeV);Events / 0.002 GeV",300,0.25,0.85);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass );
        group_34B.insert(histdef); 
        histdef.clear();

        histdef.clear();
        name="etaMassFCAL_Kin_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre*InFCAL;M(#eta) (GeV);Events / 0.002 GeV",300,0.25,0.85);
        histdef.name = name; histdef.cut=&mEllipse_pre_etaFCAL; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass_Kin );
        group_34B.insert(histdef); 
        
        histdef.clear();
        name="etaMassBCAL_Kin_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre*InBCAL;M(#eta) (GeV);Events / 0.002 GeV",300,0.25,0.85);
        histdef.name = name; histdef.cut=&mEllipse_pre_etaBCAL; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass_Kin );
        group_34B.insert(histdef); 

        histdef.clear();
        name="etaMassSPLIT_Kin_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre*SPLIT;M(#eta) (GeV);Events / 0.002 GeV",300,0.25,0.85);
        histdef.name = name; histdef.cut=&mEllipse_pre_etaSPLIT; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass_Kin );
        group_34B.insert(histdef); 

        histdef.clear();
        name="etaMassFCAL_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre*InFCAL;M(#eta) (GeV);Events / 0.002 GeV",300,0.25,0.85);
        histdef.name = name; histdef.cut=&mEllipse_pre_etaFCAL; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass );
        group_34B.insert(histdef); 
        
        histdef.clear();
        name="etaMassBCAL_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre*InBCAL;M(#eta) (GeV);Events / 0.002 GeV",300,0.25,0.85);
        histdef.name = name; histdef.cut=&mEllipse_pre_etaBCAL; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass );
        group_34B.insert(histdef); 

        histdef.clear();
        name="etaMassSPLIT_mEllipsePre";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre*SPLIT;M(#eta) (GeV);Events / 0.002 GeV",300,0.25,0.85);
        histdef.name = name; histdef.cut=&mEllipse_pre_etaSPLIT; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaMass );
        group_34B.insert(histdef); 
        histdef2d.clear();

        histdef2d.clear();
        name="pi0eta_mEllipsePre";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mEllipse_pre;#pi^{0} Mass (GeV) with Events / 0.001 GeV;#eta Mass (GeV) with Events / 0.0025 GeV",200,0.05,0.25,300,0.25,0.85);
        histdef2d.name = name; histdef2d.cut=&mEllipse_pre; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Mass_Kin );
        histdef2d.valuesY.push_back( &locEtaMass_Kin );
        group_12B_34B.insert_2D(histdef2d); 

        histdef2d.clear();
        name="pi0eta_Meas_mEllipsePre";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mEllipse_pre;#pi^{0} Mass (GeV) with Events / 0.001 GeV;#eta Mass (GeV) with Events / 0.0025 GeV",200,0.05,0.25,300,0.25,0.85);
        histdef2d.name = name; histdef2d.cut=&mEllipse_pre; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Mass );
        histdef2d.valuesY.push_back( &locEtaMass );
        group_12B_34B.insert_2D(histdef2d); 

	//for (int i=0; i<10; ++i){
        histdef2d.clear();
        name="pi0eta_mEllipsePre_rectSBregions_kin";
        //name="pi0eta_mEllipsePre_rectSBregion"+to_string(i);
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mEllipse_pre;#pi^{0} Mass (GeV) with Events / 0.001 GeV;#eta Mass (GeV) with Events / 0.0025 GeV",200,0.05,0.25,300,0.25,0.85);
        histdef2d.name = name; histdef2d.cut=&(inBox[9]); histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Mass_Kin );
        histdef2d.valuesY.push_back( &locEtaMass_Kin );
        group_12B_34B.insert_2D(histdef2d); 
	//}
        // *********************** PHOTON PAIR PLOTS ************************
        histdef.clear();
        name="3DistanceBetweenPhotonsFCAL_mdij3";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mdij3;FCAL - Distance between Pairs of Photons (cm);Events / 0.5 cm", 500,0,250);
        histdef.name = name; histdef.cut=&mdij3; histdef.weights = &weightAS;
        group_pairFCAL.insert(histdef); 
        
        histdef.clear();
        name="distanceOnCylinderBCAL_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts;BCAL - Distance on cylinder between Pairs of Photons (cm);Events / 1 cm",400,0,400);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        group_pairBCAL.insert(histdef); 

        // *********************** PI0ETA MASS PLOTS ******************************
        histdef.clear();
        name="pi0eta1D_mMandelstamT_mBeamE8GeVPlus";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mBeamE8GeVPlus;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mMandelstamT_mBeamE8GeVPlus; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_thrown_mMandelstamT_mBeamE8GeVPlus";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mBeamE8GeVPlus;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mMandelstamT_mBeamE8GeVPlus; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_thrown );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_res_mMandelstamT_mBeamE8GeVPlus";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mBeamE8GeVPlus;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 100, -0.5, 0.5);
        histdef.name = name; histdef.cut=&mMandelstamT_mBeamE8GeVPlus; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_resolution );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=allGeneralCutsPassed;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&allGeneralCutsPassed; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_ASBS";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS_BS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_ASBS_meas";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS_BS;
        histdef.values.push_back( &locPi0Eta );
        group_1234B.insert(histdef); 

	for (int iTBin=0; iTBin<10; ++iTBin){
        	histdef.clear();
        	name="pi0eta1D_Cut_tBin"+to_string(iTBin);
        	histdef.hist = new TH1F(name.c_str(), ("Cuts=mMandelstamT*p_massTBinned["+to_string(iTBin)+"];M(#pi^{0}#eta) (GeV);Events / 0.02 GeV").c_str(), 125, 0, 3.5);
        	histdef.name = name; histdef.cut=&p_massTBinned[iTBin]; histdef.weights = &weightAS_BS;
        	histdef.values.push_back( &locPi0Eta_Kin );
        	group_1234B.insert(histdef); 
	}

        histdef.clear();
        name="pi0eta1D_Cut_tAll";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre_tAll;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mEllipse_pre_tAll; histdef.weights = &weightAS_BS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
        histdef.clear();

        name="pi0eta1D_BaryBkg";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts*BaryBkg;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&allGen_barybkg; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_vanHove";
        histdef.hist = new TH1F(name.c_str(), "Cuts=GeneralCuts*vanHove;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&allGen_vanHove; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
        
        histdef.clear();
        // so in side band subraction we have 4 regions. 1=circular signal. 2=disc skip region. 3=disc bkg region. 4=large disc reject region. 
        name="pi0eta1D_0_1_1pR_1";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre -- 0_1_1+R_1;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS_B;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_1_0_mR_0";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre -- 1_0_-R_0;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS_BS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_1_1_1_1";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre -- 1_1_1_1;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_RectSBSubRegion4";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre -- RectSBSubRregion4;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&(inBox[4]); histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
        histdef.clear();
        name="pi0eta1D_RectSBSubRegion0268";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre -- RectSBSubRregion0268;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&(inBox[10]); histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
        histdef.clear();
        name="pi0eta1D_RectSBSubRegion17";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre -- RectSBSubRregion17;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&(inBox[11]); histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
        histdef.clear();
        name="pi0eta1D_RectSBSubRegion35";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre -- RectSBSubRregion35;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&(inBox[12]); histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
	


        // ************************** GENERAL KINEMATIC QUANTITIES ***************************
        histdef2d.clear();
        name="tetaVsMpi0eta";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mMandelstamT;M(#pi^{0}#eta) (GeV);t_{#eta} (GeV^2)", 260, 0.6, 3.2, 80,0,8);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &mandelstam_teta_Kin );
        group_34B_1234B.insert_2D(histdef2d); 

        histdef2d.clear();
        name="tpi0VsMpi0eta";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mMandelstamT;M(#pi^{0}#eta) (GeV);t_{#pi^{0}} (GeV^2)", 260, 0.6, 3.2, 80,0,8);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &mandelstam_tpi0_Kin );
        group_34B_1234B.insert_2D(histdef2d); 

        histdef.clear();
        name="mandelstam_tp";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT;-t' momentum transfer of #pi^{0}+#eta;Events / 0.06 GeV", 100, 0, 6);
        histdef.name = name; histdef.cut=&mMandelstamT; histdef.weights = &weightAS_BS;
        histdef.values.push_back( &mandelstam_tp );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="MissingMassSquared_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMMSq;Missing Mass Squared (GeV/c^{2})^{2};Events / 0.002 GeV/c^{2}",200,-0.2,0.2);
        histdef.name = name; histdef.cut=&mMMSq; histdef.weights = &weightAS;
        histdef.values.push_back( &locMissingMassSquared );
        group_1234BP.insert(histdef); 

        histdef.clear();
        name="BeamEnergy_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=noCut;Beam Energy (GeV) ;Events / 0.1 GeV",120,0,12);
        histdef.name = name; histdef.cut=&noCut; histdef.weights = &weightAS;
        histdef.values.push_back( &locBeamE );
        group_B.insert(histdef); 

        histdef.clear();
        name="pi0proton1D_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMPi0P14_ellipse;M(#pi^{0}proton) (GeV);Events / 0.01 GeV", 400,0,4);
        histdef.name = name; histdef.cut=&mMPi0P14_ellipse; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Proton_Kin);
        group_12PB.insert(histdef); 

	// need this one for the graph, just delete afterwards
        histdef.clear();
        name="pi0proton1D_Cut_ASBS";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMPi0P14_ellipse;M(#pi^{0}proton) (GeV);Events / 0.01 GeV", 400,0,4);
        histdef.name = name; histdef.cut=&mMPi0P14_ellipse; histdef.weights = &weightAS_BS;
        histdef.values.push_back( &locPi0Proton_Kin);
        group_12PB.insert(histdef); 
        histdef.clear();

        name="etaproton1D_Cut";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mEllipse_pre;M(#etaproton) (GeV);Events / 0.01 GeV",450,0,4.5);
        histdef.name = name; histdef.cut=&mEllipse_pre; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaProton_Kin );
        group_34PB.insert(histdef); 

        histdef2d.clear();
        name="pi0etaPi0Proton_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=GeneralCuts;M(#pi^{0}#eta) (GeV) with Events / 0.025 GeV;M(#pi^{0}Proton) (GeV) with Events / 0.05 GeV",160,0,4,90,0,4.5);
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &locPi0Proton_Kin );
        group_1234B_12PB.insert_2D(histdef2d); 

        histdef2d.clear();
        name="pi0etaEtaProton_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=GeneralCuts;M(#pi^{0}#eta) (GeV) with Events / 0.025 GeV;M(#etaProton) (GeV) with Events / 0.05 GeV",160,0,4,100,0,5);
        histdef2d.name = name; histdef2d.cut=&allGeneralCutsPassed; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &locEtaProton_Kin );
        group_1234B_34PB.insert_2D(histdef2d); 
        
        // *************** BARYON BACKGROUND HISTOGRAMS ***********************
	// The goal of this section is to show how the 4 bkg rejection histograms work. { M(pi0proton), M(etaproton), vanHove, t } 
	// So we will always keep track of M(pi0eta) to see how these cuts affect it, so 5 histograms to keep track of
	// We will in turn select one of the 4 bkg rejection histograms, cut on them, display the remaining 3+M(pi0eta) histogram on a 2x2 canvas
	// mMandelstamT_mdelta = pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        // mMandelstamT = !pMPi0P14*mMandelstamT_mdelta; // with delta cut applied
	// mMandelstamT_mdelta_petaProton = pEtaProtonBaryonCut*mMandelstamT_mdelta; 
	// mMandelstamT_mdelta_pVanHove = pVanHove*mMandelstamT_mdelta; 
	// mDelta = ptpLT1*mMandelstamT_mdelta; // t Cut applied

	// -------------------- SHOW THESE FIRST BEFORE ANY CUTTING ------------------------------
        histdef2d.clear();
        name="vanHove_mMandelstamT_mdelta";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mMandelstamT_mdelta;Events / 0.1 degrees;Events / 0.1 degrees",60,-3,3,60,-3,3);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT_mdelta; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &vanHove_x );
        histdef2d.valuesY.push_back( &vanHove_y );
        group_1234BP.insert_2D(histdef2d); 

        histdef.clear();
        name="mandelstam_tp_mMandelstamT_mdelta";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta;-t' momentum transfer of #pi^{0}+#eta;Events / 0.06 GeV", 100, 0, 6);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta; histdef.weights = &weightAS_BS;
        histdef.values.push_back( &mandelstam_tp );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0proton1D_mMandelstamT_mdelta";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta;M(#pi^{0}proton) (GeV);Events / 0.01 GeV", 400,0,4);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Proton_Kin);
        group_12PB.insert(histdef); 

        histdef.clear();
        name="etaproton1D_mMandelstamT_mdelta";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta;M(#etaproton) (GeV);Events / 0.01 GeV",450,0,4.5);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaProton_Kin );
        group_34PB.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_mMandelstamT_mdelta";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
	
	// -------------------- NOW WE USE VANHOVE AS THE PRIMARY ---------------------
        histdef.clear();
        name="mandelstam_tp_mMandelstamT_mdelta_pVanHove";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta_pVanHove;-t' momentum transfer of #pi^{0}+#eta;Events / 0.06 GeV", 100, 0, 6);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta_pVanHove; histdef.weights = &weightAS;
        histdef.values.push_back( &mandelstam_tp );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0proton1D_mMandelstamT_mdelta_pVanHove";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta_pVanHove;M(#pi^{0}proton) (GeV);Events / 0.01 GeV", 400,0,4);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta_pVanHove; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Proton_Kin);
        group_12PB.insert(histdef); 

        histdef.clear();
        name="etaproton1D_mMandelstamT_mdelta_pVanHove";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta_pVanHove;M(#etaproton) (GeV);Events / 0.01 GeV",450,0,4.5);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta_pVanHove; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaProton_Kin );
        group_34PB.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_mMandelstamT_mdelta_pVanHove";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta_pVanHove;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta_pVanHove; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
	// -------------------- NOW WE USE MANDELSTAM_T AS THE PRIMARY ---------------------
        histdef2d.clear();
        name="vanHove_mDelta";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mDelta;Events / 0.1 degrees;Events / 0.1 degrees",60,-3,3,60,-3,3);
        histdef2d.name = name; histdef2d.cut=&mDelta; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &vanHove_x );
        histdef2d.valuesY.push_back( &vanHove_y );
        group_1234BP.insert_2D(histdef2d); 

        histdef.clear();
        name="pi0proton1D_mDelta";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mDelta;M(#pi^{0}proton) (GeV);Events / 0.01 GeV", 400,0,4);
        histdef.name = name; histdef.cut=&mDelta; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Proton_Kin);
        group_12PB.insert(histdef); 

        histdef.clear();
        name="etaproton1D_mDelta";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mDelta;M(#etaproton) (GeV);Events / 0.01 GeV",450,0,4.5);
        histdef.name = name; histdef.cut=&mDelta; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaProton_Kin );
        group_34PB.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_mDelta";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mDelta;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mDelta; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
	// -------------------- NOW WE USE M(PI0PROTON) AS THE PRIMARY ---------------------
        histdef2d.clear();
        name="vanHove_mMandelstamT";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mMandelstamT;Events / 0.1 degrees;Events / 0.1 degrees",60,-3,3,60,-3,3);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &vanHove_x );
        histdef2d.valuesY.push_back( &vanHove_y );
        group_1234BP.insert_2D(histdef2d); 

        histdef.clear();
        name="mandelstam_tp_mMandelstamT";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT;-t' momentum transfer of #pi^{0}+#eta;Events / 0.06 GeV", 100, 0, 6);
        histdef.name = name; histdef.cut=&mMandelstamT; histdef.weights = &weightAS;
        histdef.values.push_back( &mandelstam_tp );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="etaproton1D_mMandelstamT";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT;M(#etaproton) (GeV);Events / 0.01 GeV",450,0,4.5);
        histdef.name = name; histdef.cut=&mMandelstamT; histdef.weights = &weightAS;
        histdef.values.push_back( &locEtaProton_Kin );
        group_34PB.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_mMandelstamT";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mMandelstamT; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
	// -------------------- NOW WE USE M(ETAPROTON) AS THE PRIMARY ---------------------
        histdef2d.clear();
        name="vanHove_mMandelstamT_mdelta_petaProton";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mMandelstamT_mdelta_petaProton;Events / 0.1 degrees;Events / 0.1 degrees",60,-3,3,60,-3,3);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT_mdelta_petaProton; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &vanHove_x );
        histdef2d.valuesY.push_back( &vanHove_y );
        group_1234BP.insert_2D(histdef2d); 

        histdef.clear();
        name="mandelstam_tp_mMandelstamT_mdelta_petaProton";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta_petaProton;-t' momentum transfer of #pi^{0}+#eta;Events / 0.06 GeV", 100, 0, 6);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta_petaProton; histdef.weights = &weightAS;
        histdef.values.push_back( &mandelstam_tp );
        group_1234B.insert(histdef); 

        histdef.clear();
        name="pi0proton1D_mMandelstamT_mdelta_petaProton";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta_petaProton;M(#pi^{0}proton) (GeV);Events / 0.01 GeV", 400,0,4);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta_petaProton; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Proton_Kin);
        group_12PB.insert(histdef); 

        histdef.clear();
        name="pi0eta1D_mMandelstamT_mdelta_petaProton";
        histdef.hist = new TH1F(name.c_str(), "Cuts=mMandelstamT_mdelta_petaProton;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 350, 0, 3.5);
        histdef.name = name; histdef.cut=&mMandelstamT_mdelta_petaProton; histdef.weights = &weightAS;
        histdef.values.push_back( &locPi0Eta_Kin );
        group_1234B.insert(histdef); 
	
	// ------------------ END OF BARYON CUT HISTGRAMS ----------------------
	// ---------------------------------------------------------------------

	// To show the t' distribution as a function of M(pi0eta)
        histdef2d.clear();
        name="mandelstam_tpVspi0etaMass_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mMandelstamT;M(#pi^{0}#eta) with Events / 0.01  GeV;-t' momentum transfer of #pi^{0}+#eta with Events / 0.06 GeV",100,0.8,1.8,100,0,6);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &mandelstam_tp );
		// I THINK WE CAN JUST USE GROUP_1234B SINCE BOTH OF THE VARAIBLES DEPENDS ON THIS TRACKING
        group_1234B.insert_2D(histdef2d); 

	// To show the t distribution as a function of M(pi0eta)
        histdef2d.clear();
        name="mandelstam_tVspi0etaMass_Cut";
        histdef2d.hist = new TH2F(name.c_str(), "Cuts=mMandelstamT;M(#pi^{0}#eta) with Events / 0.01  GeV;-t' momentum transfer of #pi^{0}+#eta with Events / 0.06 GeV",100,0.8,1.8,100,0,6);
        histdef2d.name = name; histdef2d.cut=&mMandelstamT; histdef2d.weights = &weightAS;
        histdef2d.valuesX.push_back( &locPi0Eta_Kin );
        histdef2d.valuesY.push_back( &mandelstam_abst );
		// I THINK WE CAN JUST USE GROUP_1234B SINCE BOTH OF THE VARAIBLES DEPENDS ON THIS TRACKING
        group_1234B.insert_2D(histdef2d); 







        
        //name;
        //histdef2d.hist = new TH2F(name.c_str(), );
        //histdef2d.name = name; histdef2d.cut=; histdef2d.weights = &weightAS;
        //histdef2d.valuesX.push_back( & );
        //histdef2d.valuesY.push_back( & );
        //group_PB.insert_2D(histdef2d); 

        //histdef.clear();
        //name;
        //histdef.hist = new TH1F(name.c_str(), );
        //histdef.name = name; histdef.cut=; histdef.weights = &weightAS;
        //histdef.values.push_back( & );
        //group_PB.insert(histdef); 
        
        

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
        //for(int groupNum=1; groupNum<totalGroups+1; ++groupNum){
        //	group_ids.clear();
        //	// loop through all the histograms
        //	for(int id_idx=0; id_idx<id+1; ++id_idx){
        //		if(histVals[id_idx][0]==groupNum){
        //			group_ids.push_back(id_idx);
        //			//if(showOutput){ cout << histList[id_idx][0] << " is under the group: " << std::to_string(groupNum) << endl; }
        //		}
        //	}
        //	vec_group_ids.push_back(group_ids);
        //}


        //if(onlyNamesPi0_1){ cout << "ONLY HISTOGRAMS WITH AFFIX _1 and _2 FOR MERGING PI0_1 AND PI0_2 HISTOGRAMS !!!" << endl; }
        //int hist_dim;
        //for(std::size_t i = 0; i<vec_group_ids.size(); ++i){
        //	if((!onlyNamesPi0_1)*showOutput) {cout << " hist Ids for a group: " << std::to_string(i+1) << " - " << groupNames[i] << endl << "=============================" << endl;}
        //	for(std::size_t j = 0; j<vec_group_ids[i].size(); ++j){
        //		if(onlyNamesPi0_1){
        //			if(histList[vec_group_ids[i][j]][0].find("_1") != std::string::npos) {
        //				if ( histList[vec_group_ids[i][j]].size() > 8 ){ hist_dim = 2; }
        //				else { hist_dim = 1;}
        //				if(showOutput){cout << histList[vec_group_ids[i][j]][0] << "," << hist_dim << endl;}
        //			}
        //		}
        //		else {
        //			if(showOutput){cout << histList[vec_group_ids[i][j]][0] << " with ID:" << histVals[vec_group_ids[i][j]][0] << endl;}
        //		}
        //	}
        //}

        //if(showOutput) { cout << endl << "-------------------------------" << endl << "Total Number of histograms: " << std::to_string(id+1) << endl << "----------------------------------" << endl;}

        // *****//*********************************** MAKE THE HISTOGRAMS ********************************************//
        //for(int i=0; i<id+1; ++i){// need + 1 since we have to include the last histogram.... 
        //	if(showOutput) {cout << "Building histogram: " << histList[i][0] << " with id = " << std::to_string(i) << " and group = " << std::to_string(histVals[i][0]) << endl; }
        //	int histListSize = (int)histList[i].size();
        //	if(histListSize == 6 || histListSize == 7 ){
        //		dHist_all1DHists[i] = new TH1F(histList[i][0].c_str(), histList[i][1].c_str(), atof(histList[i][2].c_str()), atof(histList[i][3].c_str()), atof(histList[i][4].c_str()));
        //		dHist_all1DHists[i]->SetYTitle(histList[i][5].c_str());
        //		// some have a title making the size 1 larger
        //		if(histListSize==7){
        //			dHist_all1DHists[i]->SetTitle(histList[i][6].c_str());
        //		}
        //		
        //	}
        //	else if (histListSize==10 || histListSize==11){
        //		dHist_all2DHists[i] = new TH2F(histList[i][0].c_str(), histList[i][1].c_str(), atof(histList[i][2].c_str()), atof(histList[i][3].c_str()), atof(histList[i][4].c_str()), 
        //			atof(histList[i][5].c_str()), atof(histList[i][6].c_str()), atof(histList[i][7].c_str()));
        //		dHist_all2DHists[i]->SetXTitle(histList[i][8].c_str());
        //		dHist_all2DHists[i]->SetYTitle(histList[i][9].c_str());
        //		// it could exist that a 2D hist has a title making total num of params = 11
        //		if((int)histList[i].size()==11){
        //			dHist_all2DHists[i]->SetTitle(histList[i][10].c_str());
        //		}
        //	}
        //	else { cout << "ERROR!" << endl  << "ERROR!" << endl << "ERROR!" << endl; }
        //}


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
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("weightASBS"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("uniqueComboID"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mandelstam_tp"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("ptGT1");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("ptLT05");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("ptGT05LT1");
        dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("finalStateComboID"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("chiSq"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("unusedEnergy"); //fundamental = char, int, float, double, etc.
        if (is_pi0eta){
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0g1"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0g2"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0g1");
        	dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0g2");

        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Meta"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0eta"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0_meas"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Meta_meas"); //fundamental = char, int, float, double, etc.
        	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0eta_meas"); //fundamental = char, int, float, double, etc.
		dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mandelstam_teta_meas");	
		dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mandelstam_teta");	
		dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mandelstam_tpi0_meas");	
		dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("mandelstam_tpi0");	
        }
        else{
            dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0"); //fundamental = char, int, float, double, etc.
            dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0pi0"); //fundamental = char, int, float, double, etc.
        }
        // introducing some variables to keep track of how much unique particles we have
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_eta");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0eta");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_eta_pi0eta");
        dFlatTreeInterface->Create_Branch_Fundamental<Bool_t>("isNotRepeated_pi0_pi0eta");
        dFlatTreeInterface->Create_Branch_Fundamental<ULong64_t>("spectroscopicComboID"); //fundamental = char, int, float, double, etc.

	// introducing variables to show the pi0/eta was detected in
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("pi0DetectedIn"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("etaDetectedIn"); //fundamental = char, int, float, double, etc.
	
        // Introduce some angles to use in the phase space distance calculation
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("cosTheta_X_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("cosTheta_X_cm_meas"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phi_X_cm"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("cosTheta_eta_gj"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("cosTheta_eta_gj_meas"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phi_eta_gj"); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("phi_eta_gj_meas"); //fundamental = char, int, float, double, etc.
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
	//
	if(showOutput) { cout << "Finished creating branch funamentals" <<endl; }




} // end of initialization

Bool_t DSelector_ver20::Process(Long64_t locEntry)
{
	++count_totEvents;
    	group_PB.clear_tracking();
    	group_12B_1234B.clear_tracking();
    	group_34B_1234B.clear_tracking();
    	group_12B.clear_tracking();  
    	group_34B.clear_tracking();  
    	group_12B_34B.clear_tracking();
    	group_1234B.clear_tracking();
    	group_PhNB.clear_tracking();
    	group_12PhNB.clear_tracking();
    	group_34PhNB.clear_tracking();
    	group_pairFCAL.clear_tracking();
    	group_pairBCAL.clear_tracking();
    	group_1234BP.clear_tracking();
    	group_12PB.clear_tracking();       
    	group_34PB.clear_tracking();       
    	group_B.clear_tracking();
    	group_1234B_12PB.clear_tracking();
    	group_1234B_34PB.clear_tracking();

    	// everytime we start a new event we have to reset the tracking
    	set< map<Particle_t, set<Int_t> > > used1234B;
    	set< pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > > used12B_1234B;
    	set< pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > > used34B_1234B;
    	set< map<Particle_t, set<Int_t> > > uniqueSpectroscopicID;
    	map< map<Particle_t, set<Int_t> >, Int_t > map_uniqueSpectroscopicID;
    	set< map<Particle_t, set<Int_t> > > used12B;
    	set< map<Particle_t, set<Int_t> > > used34B;

    	set< map<Particle_t, set<Int_t> > > used123B;
    	set< map<Particle_t, set<Int_t> > > used124B;

    	//if(itersToRun<1000000){ ++itersToRun; //so we can just try to show the outut of one event 
	++count_events;
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
    	//
	paddedRun=to_string(Get_RunNumber());
	digitsInRun=paddedRun.length();
	for ( int idigit=0; idigit < maxDigitsInRun-digitsInRun; ++idigit){
		paddedRun = "0" + paddedRun;
	}
	paddedEvent=to_string(Get_EventNumber());
	digitsInEvent=paddedEvent.length();
	for ( int idigit=0; idigit < maxDigitsInEvent-digitsInEvent; ++idigit){
		paddedEvent = "0" + paddedEvent;
	}

	std::vector<int> parentArray;
	TLorentzVector etaP4;
	TLorentzVector pi0P4;
	TLorentzVector pi0etaP4;
	std::vector<int> pids;
	int locNumThrown = Get_NumThrown();

	/************************************************* PARSE THROWN TOPOLOGY ***************************************/
	TString locThrownTopology = Get_ThrownTopologyString();

	 // WE HAVE TO CHECK IF THERE IS THROWN DATA FIRST. I USE THIS CONDITION TO DETERMINE IF THE TREE IS MC OR DATA
	if (Get_NumThrown() != 0 ) {
		for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
		{	
			//Set branch array indices corresponding to this particle
			dThrownWrapper->Set_ArrayIndex(loc_i);

			//Do stuff with the wrapper here ...
			Particle_t locPID = dThrownWrapper->Get_PID();
			Int_t locParentPID = dThrownWrapper->Get_ParentIndex();

			parentArray.push_back(locParentPID);
			pids.push_back(locPID);
			if (locPID==7 && locParentPID==-1) { pi0P4 = dThrownWrapper->Get_P4(); } 
			if (locPID==17 && locParentPID==-1) { etaP4 = dThrownWrapper->Get_P4(); } 
		}

		std::vector<int> parents;
		findParents(parentArray,parents);
		//cout << "\nTHESE ARE THE PARENTS" << endl;
		for ( auto parent=0; parent < (int)parents.size(); ++parent){
			//cout << parents[parent] << endl;
		}

		int pi0ToNGamma=0;
		int etaToNGamma=0;
		bool correctFinalState=false;
		for (auto parent : parents){
			std::vector<int> daughters;
			cout << "Parent: " << parent << " which has PID=" << pids[parent] << " has children:" << endl;
			findDaughters( parentArray, daughters, parent );
			for ( auto daughter=0; daughter < (int)daughters.size(); ++daughter) {
				if ( pids[parent]==7 && pids[daughters[daughter]] == 1 && daughters.size()==2){ ++pi0ToNGamma; }
				if ( pids[parent]==17 && pids[daughters[daughter]] == 1 && daughters.size()==2){ ++etaToNGamma; }
			}
		}
		if ( pi0ToNGamma==2 && etaToNGamma==2) {
			//correctFinalState=true;
			//cout << "THIS EVENT HAS 4 GAMMA FINAL STATE!" << endl;
			pi0etaP4 = pi0P4+etaP4;
			locPi0Eta_thrown = pi0etaP4.M();
			++count_correctTopology;
		}
		else { 
        		if(showOutput) { cout << "\n\n\n*********************************************************\n**************************************************\n########    EventIdx: " << (eventIdx) << "    #############" << endl; }
			++eventIdx;
			cout << "SKIPPING THIS COMBO SINCE IT DOESN'T HAVE THE 4 GAMMA TOPOLOGY" << endl;
			cout << "THESE ARE THE PARTICLES:" << endl;
			for ( auto pid : pids ) {
				cout << pid << " ";
			}
			cout << endl;
			if (!showThrownTopology) { 
        			dComboWrapper->Set_IsComboCut(true); return kTRUE;
			}
		}
	}




    /******************************************** GET POLARIZATION ORIENTATION ******************************************/

    //Only if the run number change
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
    //So, for each quantity you histogram, keep track of what particles you used (for a givMpi0g2//Then for each combo, just compare to what you used before, and make sure it's unique


    //EXmap<Particle_t, AMPLE 2: Combo-specific info:
    //In general: Could have multiple particles with the same PID: Use a set of Int_t's
    //In general: Multiple PIDs, so multiple sets: Contain within a map
    //Multiple combos: Contain maps within a set (easier, faster to search)
    //use set<map .... > > even for reactions with the only one type of final state particles as in the case of 4 gammas
    // MAKING SURE ALL THE SETS ARE EMPTY SO WE DONT START A LOOP WITH ELEMENTS IN IT ALREADY...

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


        std::vector< TLorentzVector > allPhotonP4Vectors_Kin;
        allPhotonP4Vectors_Kin.push_back(locPhoton1P4_Kin);
        allPhotonP4Vectors_Kin.push_back(locPhoton2P4_Kin);
        allPhotonP4Vectors_Kin.push_back(locPhoton3P4_Kin);
        allPhotonP4Vectors_Kin.push_back(locPhoton4P4_Kin);
        std::vector< TLorentzVector > allPhotonX4Vectors_Kin;
        allPhotonX4Vectors_Kin.push_back(locPhoton1X4_Kin);
        allPhotonX4Vectors_Kin.push_back(locPhoton2X4_Kin);
        allPhotonX4Vectors_Kin.push_back(locPhoton3X4_Kin);
        allPhotonX4Vectors_Kin.push_back(locPhoton4X4_Kin);
        std::vector< TLorentzVector > allPhoton4Vectors_Shower;
        allPhoton4Vectors_Shower.push_back(locPhoton1X4_Shower);
        allPhoton4Vectors_Shower.push_back(locPhoton2X4_Shower);
        allPhoton4Vectors_Shower.push_back(locPhoton3X4_Shower);
        allPhoton4Vectors_Shower.push_back(locPhoton4X4_Shower);
        std::vector< DNeutralParticleHypothesis* > allPhotonWrappers;
        allPhotonWrappers.push_back(dPhoton1Wrapper);
        allPhotonWrappers.push_back(dPhoton2Wrapper);
        allPhotonWrappers.push_back(dPhoton3Wrapper);
        allPhotonWrappers.push_back(dPhoton4Wrapper);
        std::vector<Int_t> photonIds;
        photonIds.push_back(locPhoton1NeutralID);
        photonIds.push_back(locPhoton2NeutralID);
        photonIds.push_back(locPhoton3NeutralID);
        photonIds.push_back(locPhoton4NeutralID);

	int nPhotons = (int) (allPhotonP4Vectors_Kin.size());
	if(showOutput){cout << "There are " << nPhotons << " photons" << endl;  }

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

	massGammaPi0[0] = (locPhoton1P4_Kin + locPhoton2P4_Kin + locPhoton3P4_Kin).M();
	massGammaPi0[1] = (locPhoton1P4_Kin + locPhoton2P4_Kin + locPhoton4P4_Kin).M();
	massGammaEta[0] = (locPhoton3P4_Kin + locPhoton4P4_Kin + locPhoton1P4_Kin).M();
	massGammaEta[1] = (locPhoton3P4_Kin + locPhoton4P4_Kin + locPhoton2P4_Kin).M();


        if(showOutput){cout << "Combined 4-vectors" << endl;}

        /******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

        // Loop through the analysis actions, executing them in order for the active particle combo

        // Perform_Action gives an error  with our current code probably because it does not have any actions to perform
        //dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
        if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
            continue;


        //if you manually execute any actions, and it fails a cut, be sure to call:
        //dComboWrapper->Set_IsComboCut(true);

	TLorentzVector locPi0EtaP4 = locPhoton1P4 + locPhoton2P4 + locPhoton3P4 + locPhoton4P4;
	TLorentzVector locPi0EtaP4_Kin = locPhoton1P4_Kin + locPhoton2P4_Kin + locPhoton3P4_Kin + locPhoton4P4_Kin;

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
	locPi0Eta_resolution = locPi0Eta_Kin-locPi0Eta_thrown;
	cout << "locPi0Eta_Kin: " << locPi0Eta_Kin << endl;
	cout << "locPi0Eta_thrown: " << locPi0Eta_thrown << endl;
	cout << "locPi0Eta_resolution: " << locPi0Eta_resolution << endl;
        locPi0Eta = mixingPi0Eta.M();

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


	if (showOutput){ cout << "Filling photon variables" << endl; } 

        TLorentzVector newPhotonP4Vector;
        TLorentzVector newPhotonX4Vector;
        TLorentzVector newPhotonVector_Shower;
        DNeutralParticleHypothesis* newPhotonWrapper;
        for (auto iPhoton=0; iPhoton<nPhotons ; ++iPhoton){
              newPhotonP4Vector = allPhotonP4Vectors_Kin[iPhoton];
              newPhotonX4Vector = allPhotonX4Vectors_Kin[iPhoton];
              newPhotonVector_Shower = allPhoton4Vectors_Shower[iPhoton];
              newPhotonWrapper = allPhotonWrappers[iPhoton];
              photonThetas[iPhoton] = newPhotonP4Vector.Theta()*radToDeg;
              photonPhis[iPhoton] = newPhotonP4Vector.Phi()*radToDeg;
              photonEnergies[iPhoton] = newPhotonP4Vector.E();
              photonXs_Kin[iPhoton] = newPhotonX4Vector.X();
              photonYs_Kin[iPhoton] = newPhotonX4Vector.Y();
              photonZs_Kin[iPhoton] = newPhotonX4Vector.Z();
              photonTs_Kin[iPhoton] = newPhotonX4Vector.T();
              photonXs_Shower[iPhoton] = newPhotonVector_Shower.X();
              photonYs_Shower[iPhoton] = newPhotonVector_Shower.Y();
              photonZs_Shower[iPhoton] = newPhotonVector_Shower.Z();
              photonTs_Shower[iPhoton] = newPhotonVector_Shower.T();
              photonDeltaTs[iPhoton] = newPhotonX4Vector.T()-(dComboWrapper->Get_RFTime() + (newPhotonX4Vector.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
              photonDetectedSyss[iPhoton] = newPhotonWrapper->Get_Detector_System_Timing();
              E1E9_FCAL[iPhoton] = newPhotonWrapper->Get_E1E9_FCAL();
              E9E25_FCAL[iPhoton] = newPhotonWrapper->Get_E9E25_FCAL();
              SumU_FCAL[iPhoton] = newPhotonWrapper->Get_SumU_FCAL();
              SumV_FCAL[iPhoton] = newPhotonWrapper->Get_SumV_FCAL();
              Energy_BCALPreshower[iPhoton] = newPhotonWrapper->Get_Energy_BCALPreshower();
              Energy_BCAL[iPhoton] = newPhotonWrapper->Get_Energy_BCAL();
              SigLong_BCAL[iPhoton] = newPhotonWrapper->Get_SigLong_BCAL();
              SigTheta_BCAL[iPhoton] = newPhotonWrapper->Get_SigTheta_BCAL();
              showerQuality_FCAL[iPhoton] = newPhotonWrapper->Get_Shower_Quality();
              SigTrans_BCAL[iPhoton] = newPhotonWrapper->Get_SigTrans_BCAL();
              DeltaPhi_BCAL[iPhoton] = newPhotonWrapper->Get_TrackBCAL_DeltaPhi();
              DeltaZ_BCAL[iPhoton] = newPhotonWrapper->Get_TrackBCAL_DeltaZ();
              DOCA_FCAL[iPhoton] = newPhotonWrapper->Get_TrackFCAL_DOCA();
        }

        if(showOutput) { cout << "Getting some combo specific variables like CL" << endl; }
        //Kinematic fit variables
        locCLKinFit = dComboWrapper->Get_ConfidenceLevel_KinFit("");
        locUnusedEnergy = dComboWrapper->Get_Energy_UnusedShowers();
        locNumExtraNeutralShowers = Get_NumNeutralHypos()-4;
        locDOFKinFit = dComboWrapper->Get_NDF_KinFit("");
        locChiSqKinFit = dComboWrapper->Get_ChiSq_KinFit("");
        //TString DOF; DOF.Form("%d", dComboWrapper->Get_NDF_KinFit(""));
        //TString title = "Num DOF: "; title+=DOF.Data();

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
        //std::vector<double> angle_ijVec;
        //std::vector<double> deltaZ_ijVec;
        //std::vector<double> deltaPhi_ijVec;
        std::vector< double> dij3VecFCAL;
        std::vector< double> dijVecBCAL;
        std::vector< map<Particle_t, set<Int_t> > > usedPairIdsFCAL;// will be used to do our uniqueness trackingo over the entire event
        std::vector< map<Particle_t, set<Int_t> > > usedPairIdsBCAL;

        map<Particle_t, set<Int_t> > usingPairIdsFCAL;
        map<Particle_t, set<Int_t> > usingPairIdsBCAL; 



        if(showOutput){ cout << "\nCalculating pairwise quantities between photons, doing uniqueness tracking now instead of later since there will be a vector of data to fill\n-----------------------------------------" << endl;}
        std::vector<TLorentzVector> photonX4 = {locPhoton1X4_Shower, locPhoton2X4_Shower, locPhoton3X4_Shower, locPhoton4X4_Shower};
        if(showOutput){ cout << "Set up photonX4 and photonIds to use for pairwise quantity calculating!" << endl;}
        for (int i=0;i<4;++i){
        	for(int j=i+1; j<4;j++){
        		if(showOutput){ cout << "-- ith,jth photon: " << std::to_string(i) << ", " << std::to_string(j) << endl;}
        		// if both the photons are not in the resppective detector there is no point on continuing to check the histograms + cuts + uniquenss tracking
        		if(photonDetectedSyss[i]==SYS_FCAL && photonDetectedSyss[j]==SYS_FCAL){
        			++countBothInFCAL;
        			usingPairIdsFCAL.clear();
                                usingPairIdsFCAL[Unknown].insert(locBeamID);
        			usingPairIdsFCAL[Gamma].insert(photonIds[i]);
        			usingPairIdsFCAL[Gamma].insert(photonIds[j]);
                        	dij3 = (photonX4[i]-photonX4[j]).Vect().Mag();
                                //if (usedPairIdsFCAL_test.find(usingPairIdsFCAL)==usedPairIdsFCAL_test.end()){
                                //        usedPairIdsFCAL_test.insert(usingPairIdsFCAL);
        			//        dij3VecFCAL.push_back(dij3);
                                //}
        			dij3VecFCAL.push_back(dij3);
                                usedPairIdsFCAL.push_back(usingPairIdsFCAL);
        		}
        		else if(photonDetectedSyss[i]==SYS_BCAL && photonDetectedSyss[j]==SYS_BCAL){
        			++countBothInBCAL;
        			usingPairIdsBCAL.clear();
                                usingPairIdsBCAL[Unknown].insert(locBeamID);
        			usingPairIdsBCAL[Gamma].insert(photonIds[i]);
        			usingPairIdsBCAL[Gamma].insert(photonIds[j]);
                                usedPairIdsBCAL.push_back(usingPairIdsBCAL);
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

                                cout << "BCAL dij: " << dij << endl;
        			dijVecBCAL.push_back(dij);
        			//angle_ijVec[hist_id].push_back(angle_ij);
        			//deltaZ_ijVec[hist_id].push_back(deltaZ_ij);
        		}
        		else if((photonDetectedSyss[i]==SYS_BCAL && photonDetectedSyss[j]==SYS_FCAL) || (photonDetectedSyss[i]==SYS_FCAL && photonDetectedSyss[j]==SYS_BCAL)){
                                cout << "SPLIT" << endl;
        			++countNotInEither;
        			if(showOutput){ cout << "ONE IN FCAL ONE IN BCAL" << endl;}
        		}
        		else { cout << "\n\n*************************\nERROR - Pairwise both in FCAL, both in FCAL, in either, is not complete!\n***********************************" << endl;}
        	}
        }
        if((countBothInBCAL+countNotInEither+countBothInFCAL)==6){ if(showOutput) {cout << "\n\n*********************\nPairwise photons counted correctly!\n*******************\nNum both in FCAL: " + std::to_string(countBothInFCAL) +"\nNum both in BCAL: " + std::to_string(countBothInBCAL) + "\nNum one in either: " + std::to_string(countNotInEither)<< endl; }}
        else { cout << "\n\n*********************\nERROR - Pairwise photons not counted correctly!\n*******************\n" << endl;}

        // *********************** PHOTON PAIR PLOTS ************************
        for (UInt_t iPair=0; iPair<dij3VecFCAL.size(); ++iPair){
           //cout << "FCAL dij3: " << dij3 << endl;
           //for ( auto elem : usedPairIdsFCAL[iPair] ){
           //    cout << elem.first << " ";
           //    for (auto it=elem.second.begin(); it != elem.second.end(); ++it){
           //        cout << *it << " ";
           //    }
           //    cout << endl;
           //}
           group_pairFCAL.allHists_1D[0].values.push_back( &dij3VecFCAL[iPair]);
        }
        for (UInt_t iPair=0; iPair<dijVecBCAL.size(); ++iPair){
            group_pairBCAL.allHists_1D[0].values.push_back( &dijVecBCAL[iPair]);
        }

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
	mandelstam_abst = abs(mandelstam_t);
        mandelstam_t_pe = (locBeamP4_Kin-mixingPi0Eta).M2();
        mandelstam_teta = -(locBeamP4-locEtaP4).M2();
        mandelstam_tpi0 = -(locBeamP4-locPi0P4).M2();
        mandelstam_teta_Kin = -(locBeamP4_Kin-locEtaP4_Kin).M2();
        mandelstam_tpi0_Kin = -(locBeamP4_Kin-locPi0P4_Kin).M2();
	idx_t_eta = (int)( (mandelstam_teta-tMin)/tStep ); 
	idx_t_pi0 = (int)( (mandelstam_tpi0-tMin)/tStep ); 
	idx_m = (int)( (locPi0Eta-mMin)/mStep );
	if ( mandelstam_teta < tMin || locPi0Eta < mMin || mandelstam_teta>tMax || locPi0Eta>mMax ) {
		teta_recCounts = -1;
	}
	else {
		teta_recCounts = idx_t_eta+num_tBins*idx_m;
	}
	if ( mandelstam_tpi0 < tMin || locPi0Eta < mMin || mandelstam_tpi0>tMax || locPi0Eta>mMax ) {
		tpi0_recCounts = -1;
	}
	else {
		tpi0_recCounts = idx_t_pi0+num_tBins*idx_m;
	}
	cout << "teta_recCounts: " << teta_recCounts << endl;
	cout << "tpi0_recCounts: " << tpi0_recCounts << endl;

        // We will also calculate the cos theta in the lab frame. This is used to see if we are actually getting the right amount of events in the detector
        theta_pi0_lab = locPi0P4.Theta()*radToDeg;
        theta_eta_lab = locEtaP4.Theta()*radToDeg;

        // Determine the CM vectors to be used to calculate the angles in the Hel frame
        // also calculating the cosTheta of pi0eta system, pi0, eta  in the CM frame
        TLorentzVector cm_vec = locBeamP4_Kin+dTargetP4;
        TLorentzVector cm_vec_meas = locBeamP4+dTargetP4;
        TLorentzVector mixingPi0Eta_cm = mixingPi0Eta_Kin;
        TLorentzVector mixingPi0Eta_cm_meas = mixingPi0Eta;
        TLorentzVector pi0_cm = locPi0P4_Kin;
        TLorentzVector pi0_cm_meas = locPi0P4;
        TLorentzVector eta_cm = locEtaP4_Kin;
        TLorentzVector eta_cm_meas = locEtaP4;
        TLorentzVector beam_cm = locBeamP4_Kin;
        TLorentzVector beam_cm_meas = locBeamP4;
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
        mixingPi0Eta_cm_meas.Boost(-cm_vec_meas.BoostVector());
        pi0_cm.Boost(-cm_vec.BoostVector());
        eta_cm.Boost(-cm_vec.BoostVector());
        eta_cm_meas.Boost(-cm_vec_meas.BoostVector());
        beam_cm.Boost(-cm_vec.BoostVector());
        beam_cm_meas.Boost(-cm_vec_meas.BoostVector());
        recoil_cm.Boost(-cm_vec.BoostVector());

        cosTheta_pi0eta_CM = mixingPi0Eta_cm.CosTheta();
        cosTheta_pi0eta_CM_meas = mixingPi0Eta_cm_meas.CosTheta();
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
	cout << "cosTheta_pi0eta_hel: " << cosTheta_pi0eta_hel << endl;
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
        TLorentzVector beam_res_meas = beam_cm_meas;
        TLorentzVector recoil_res = recoil_cm;
        TLorentzVector eta_res = eta_cm;
        TLorentzVector eta_res_meas = eta_cm_meas;
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
        phi_pi0_GJ = angles_pi0.Phi()*radToDeg;
        phi_eta_GJ = angles_eta.Phi()*radToDeg;	

	// some extra variables from the measured 4 vectors
        beam_res_meas.Boost(-mixingPi0Eta_cm_meas.BoostVector());
       	eta_res_meas.Boost(-mixingPi0Eta_cm_meas.BoostVector());
        z = beam_res_meas.Vect().Unit();
        y = mixingPi0Eta_cm_meas.Vect().Cross(beam_cm_meas.Vect()).Unit();
        locYDotZ_GJ=y.Dot(z);
        x = y.Cross(z).Unit();
        TVector3 eta_res_unit_meas = eta_res_meas.Vect().Unit();	
        angles_eta.SetXYZ ( eta_res_unit_meas.Dot(x), eta_res_unit_meas.Dot(y), eta_res_unit_meas.Dot(z) );
        cosTheta_eta_GJ_meas = angles_eta.CosTheta();
        phi_eta_GJ_meas = angles_eta.Phi()*radToDeg;	

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

        map<Particle_t, set<Int_t> > using12B;
        using12B[Unknown].insert(locBeamID); //beam
        using12B[Gamma].insert(locPhoton1NeutralID);
        using12B[Gamma].insert(locPhoton2NeutralID);

        map<Particle_t, set<Int_t> > using34B;
        using34B[Unknown].insert(locBeamID); //beam
        using34B[Gamma].insert(locPhoton3NeutralID);
        using34B[Gamma].insert(locPhoton4NeutralID);


        map<Particle_t, set<Int_t> > using13B;
        using13B[Unknown].insert(locBeamID); //beam
        using13B[Gamma].insert(locPhoton1NeutralID);
        using13B[Gamma].insert(locPhoton3NeutralID);

        map<Particle_t, set<Int_t> > using12;
        using12[Gamma].insert(locPhoton1NeutralID);
        using12[Gamma].insert(locPhoton2NeutralID);

        map<Particle_t, set<Int_t> > using24B;
        using24B[Unknown].insert(locBeamID); //beam
        using24B[Gamma].insert(locPhoton2NeutralID);
        using24B[Gamma].insert(locPhoton4NeutralID);

        map<Particle_t, set<Int_t> > using34;
        using34[Gamma].insert(locPhoton3NeutralID);
        using34[Gamma].insert(locPhoton4NeutralID);

        map<Particle_t, set<Int_t> > using1234B;
        using1234B[Unknown].insert(locBeamID); //beam
        using1234B[Gamma].insert(locPhoton1NeutralID);
        using1234B[Gamma].insert(locPhoton2NeutralID);
        using1234B[Gamma].insert(locPhoton3NeutralID);
        using1234B[Gamma].insert(locPhoton4NeutralID);

        map<Particle_t, set<Int_t> > using1234;
        using1234[Gamma].insert(locPhoton1NeutralID);
        using1234[Gamma].insert(locPhoton2NeutralID);
        using1234[Gamma].insert(locPhoton3NeutralID);
        using1234[Gamma].insert(locPhoton4NeutralID);

        map<Particle_t, set<Int_t> > using12PB;
        using12PB[Unknown].insert(locBeamID); //beam
        using12PB[Proton].insert(locProtonTrackID);
        using12PB[Gamma].insert(locPhoton1NeutralID);
        using12PB[Gamma].insert(locPhoton2NeutralID);

        map<Particle_t, set<Int_t> > using34PB;
        using34PB[Unknown].insert(locBeamID); //beam
        using34PB[Proton].insert(locProtonTrackID);
        using34PB[Gamma].insert(locPhoton3NeutralID);
        using34PB[Gamma].insert(locPhoton4NeutralID);

        map<Particle_t, set<Int_t> > using1234PB;
        using1234PB[Unknown].insert(locBeamID); //beam
        using1234PB[Proton].insert(locProtonTrackID);
        using1234PB[Gamma].insert(locPhoton1NeutralID);
        using1234PB[Gamma].insert(locPhoton2NeutralID);
        using1234PB[Gamma].insert(locPhoton3NeutralID);
        using1234PB[Gamma].insert(locPhoton4NeutralID);

        std::vector< map<Particle_t, set<Int_t>> > beingUsedNeutralIds; 
        map<Particle_t, set<Int_t>> usingPhB;
        usingPhB[Unknown].insert(locBeamID);
        usingPhB[Gamma].insert(locPhoton1NeutralID);
        beingUsedNeutralIds.push_back( usingPhB );
        usingPhB.clear();
        usingPhB[Unknown].insert(locBeamID);
        usingPhB[Gamma].insert(locPhoton2NeutralID);
        beingUsedNeutralIds.push_back( usingPhB );
        usingPhB.clear();
        usingPhB[Unknown].insert(locBeamID);
        usingPhB[Gamma].insert(locPhoton3NeutralID);
        beingUsedNeutralIds.push_back( usingPhB );
        usingPhB.clear();
        usingPhB[Unknown].insert(locBeamID);
        usingPhB[Gamma].insert(locPhoton4NeutralID);
        beingUsedNeutralIds.push_back( usingPhB );

	std::vector< map<Particle_t, set<Int_t>> > using12PhNBs;
	map<Particle_t, set<Int_t>> using12PhNB;
	using12PhNB[Unknown].insert(locBeamID);
	using12PhNB[Gamma].insert(locPhoton1NeutralID);
	using12PhNB[Gamma].insert(locPhoton2NeutralID);
	using12PhNB[Gamma].insert(locPhoton3NeutralID);
	using12PhNBs.push_back(using12PhNB);
	using12PhNB.clear();
	using12PhNB[Unknown].insert(locBeamID);
	using12PhNB[Gamma].insert(locPhoton1NeutralID);
	using12PhNB[Gamma].insert(locPhoton2NeutralID);
	using12PhNB[Gamma].insert(locPhoton4NeutralID);
	using12PhNBs.push_back(using12PhNB);

	std::vector< map<Particle_t, set<Int_t>> > using34PhNBs;
	map<Particle_t, set<Int_t>> using34PhNB;
	using34PhNB[Unknown].insert(locBeamID);
	using34PhNB[Gamma].insert(locPhoton3NeutralID);
	using34PhNB[Gamma].insert(locPhoton4NeutralID);
	using34PhNB[Gamma].insert(locPhoton1NeutralID);
	using34PhNBs.push_back(using34PhNB);
	using34PhNB.clear();
	using34PhNB[Unknown].insert(locBeamID);
	using34PhNB[Gamma].insert(locPhoton3NeutralID);
	using34PhNB[Gamma].insert(locPhoton4NeutralID);
	using34PhNB[Gamma].insert(locPhoton2NeutralID);
	using34PhNBs.push_back(using34PhNB);

        pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > using12B_1234B = make_pair(using12B,using1234B);
        pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > using34B_1234B = make_pair(using34B,using1234B);
        pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > using1234B_12PB = make_pair(using1234B,using12PB);
        pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > using1234B_34PB = make_pair(using1234B,using34PB);
        pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > using12B_34B = make_pair(using12B,using34B);

        map<Particle_t, set<Int_t> > usingPB;
        usingPB[Unknown].insert(locBeamID);
        usingPB[Proton].insert(locProtonTrackID);
        

        map<Particle_t, set<Int_t> > usingB;
        usingB[Unknown].insert(locBeamID);

 
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
        pBeamE8to9 = locBeamE > 8.0 && locBeamE< 9.0;

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
        if (locPi0Mass_Kin<=(ellipseX-ellipseXr_loose) || locPi0Mass_Kin>=(ellipseX+ellipseXr_loose)){ outsideEllipse_loose = 1;}
        double ellipseDeltaY = TMath::Sqrt((1 - TMath::Power((locPi0Mass_Kin-ellipseX)/ellipseXr_loose,2))*TMath::Power(ellipseYr_loose,2));
        if (locEtaMass_Kin>=(ellipseDeltaY+ellipseY) || locEtaMass_Kin<=(-ellipseDeltaY+ellipseY)){ outsideEllipse_loose = 1; }
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
            weightBS = -areaRatio; 
            weightB = 1+areaRatio;
        } 
        else if (pinsideEllipse){
            dHist_checkEllipseBS[1]->Fill(locPi0Mass_Kin, locEtaMass_Kin);
            weightBS=1;
            weightB=0;
        }
        else { 
            weightBS=0;
            weightB=1;
            dHist_checkEllipseBS[2]->Fill(locPi0Mass_Kin, locEtaMass_Kin);
        }

        // now that we have defined both the weights we can multiply them together
        //weight = weightAS*weightBS;
        weightAS_BS = weightAS*weightBS;
	if ( abs(weightAS_BS)>1){
		cout << "weightBS,weightAS: " << weightBS << ", " << weightAS << endl;
	}
	// weightB would be the conjugate of the weight such that weightB+weightBS=(1,1,1,1)
        weightAS_B = weightAS*weightB;

        // General Cuts
        pUnusedEnergy = locUnusedEnergy <= unusedEnergyCut;
	pLooseUnusedEnergy = locUnusedEnergy <= 0.5;
        //pCLKinFit1 = locCLKinFit >= CLCut1;
        //pCLKinFit = locChiSqKinFit <= ChiSqCut;	
        pChiSq = locChiSqKinFit <= ChiSqCut;	
	pLooseChiSq = locChiSqKinFit <= 500;
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
        for (std::size_t i=0; i<dij3VecFCAL.size();i++){
            // since we want all the dij3FCAL to be <= 12.5 cm we can check if any of them is > 12.5 since its simpler. 
            if ( dij3VecFCAL[i] <= dijCut){
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

	if ( pPi0InFCAL ) { pi0DetectedIn = 0; } 
	if ( pPi0InBCAL ) { pi0DetectedIn = 1; } 
	if ( pPi0InSplit ) { pi0DetectedIn = 2; } 
	if ( pEtaInFCAL ) { etaDetectedIn = 0; } 
	if ( pEtaInBCAL ) { etaDetectedIn = 1; } 
	if ( pEtaInSplit ) { etaDetectedIn = 2; } 

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

        ptpLT1 = mandelstam_tp<1; 
        ptGT1 = mandelstam_abst>1; 
        ptLT05 = mandelstam_abst<0.5; 
        ptGT05LT1 = mandelstam_abst<1 && mandelstam_abst>0.5; 



        // cut to accept the delta peak in M(pi0proton)
        // we will actually use it to reject the delta peak in the cuts. pi0pi0 should also have this resonance right?
        pMPi0P14 = locPi0Proton_Kin<1.4;

        // Various combinations of cuts, the majority of them will be used just for a few histograms like when showing unused energy graph we will use mUE which
        // removes the UE cut from allGeneralCutsPassed. m prefix basically stands for minus
        //
        dzRP = pMagP3Proton*pzCutmin*pRProton;
        dzR = pzCutmin*pRProton;

        //pShowerQuality=pShowerQuality0*pShowerQuality1*pShowerQuality2*pShowerQuality3;
	pShowerQuality=true;
	pdij3pass=true;
	pUnusedEnergy=true;


	
	baseCuts = pShowerQuality*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
        allGeneralCutsPassed = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
	mMandelstamT_mBeamE8GeVPlus = !pMPi0P14*pShowerQuality*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;

	// -----------------------------------------------------------------------
	// This will be the basis of the cuts for the baryon rejection histograms
	// -----------------------------------------------------------------------
	mMandelstamT_mdelta = pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mMandelstamT = !pMPi0P14*mMandelstamT_mdelta; // with delta cut applied
	mMandelstamT_mdelta_petaProton = pEtaProtonBaryonCut*mMandelstamT_mdelta; 
	mMandelstamT_mdelta_pVanHove = pVanHove*mMandelstamT_mdelta; 
	mDelta = ptpLT1*mMandelstamT_mdelta; // t Cut applied
	// -----------------------------------------------------------------------
	// -----------------------------------------------------------------------
	
        mMPi0P14_ellipse = ptpLT1*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;

        mBeamE = ptpLT1*!pMPi0P14*pShowerQuality*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mMMSq = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pdEdxCDCProton*pinsideEllipse;
        //pDiffCL = pBeamE8GeVPlus*pUnusedEnergy*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pdEdxCDCProton*pinsideEllipse*pMissingMassSquared;
        mRProton = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mRProtonZMin = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mdEdxCDC = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pinsideEllipse;
        mZMin = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mMagP3 = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mPhotonE = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mPhotonTheta = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mdij3 = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mUE = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mUEChiSq = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton*pinsideEllipse;
        mChiSq = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pdEdxCDCProton*pinsideEllipse*pMissingMassSquared;
        // ------ 
        mEllipse_pre = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
        mEllipse_pre_tAll = !pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
	mEllipse_pre_tGT1 = mEllipse_pre_tAll*ptGT1;
	mEllipse_pre_tLT05 = mEllipse_pre_tAll*ptLT05;
	mEllipse_pre_tGT05LT1 = mEllipse_pre_tAll*ptGT05LT1;

        mEllipse_pre_tAll_delta = pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
	for ( int iTBin=0; iTBin< 10; ++iTBin ) {
		p_massTBinned[iTBin] = ( (mandelstam_tp < 0.1*(iTBin+1))  && (mandelstam_tp > 0.1*(iTBin)) ) * mMandelstamT;	
	}
        mEllipseUE_pre = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pChiSq*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
        mEllipseUEChiSq_pre = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
        mEllipseLooseUEChiSq_pre = pLooseUnusedEnergy*pLooseChiSq*ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
        mEllipseChiSq_pre = ptpLT1*!pMPi0P14*pShowerQuality*pBeamE8GeVPlus*pUnusedEnergy*pdij3pass*pPhotonE*pPhotonTheta*pMagP3Proton*pzCutmin*pRProton*pMissingMassSquared*pdEdxCDCProton;
        // ------ Rejects the 0-weight region since that would mess with the uniqueness tracking. We use the OR logic to sum the subsets of the yellow and the red region.
        mEllipse = mEllipse_pre*pYellowBKG || allGeneralCutsPassed;
        mEllipseUE = mEllipseUE_pre*pYellowBKG || mEllipseUE_pre*pinsideEllipse;
        mEllipseUEChiSq = mEllipseUEChiSq_pre*pYellowBKG || mEllipseUEChiSq_pre*pinsideEllipse;
        mEllipseChiSq = mEllipseChiSq_pre*pYellowBKG || mEllipseChiSq_pre*pinsideEllipse;

        mEllipse_pre_pi0BCAL = mEllipse_pre*pPi0InBCAL;
        mEllipse_pre_pi0FCAL = mEllipse_pre*pPi0InFCAL;
        mEllipse_pre_pi0SPLIT = mEllipse_pre*pPi0InSplit;
        mEllipse_pre_etaBCAL = mEllipse_pre*pEtaInBCAL;
        mEllipse_pre_etaFCAL = mEllipse_pre*pEtaInFCAL;
        mEllipse_pre_etaSPLIT = mEllipse_pre*pEtaInSplit;

        allGen_barybkg=allGeneralCutsPassed*pEtaProtonBaryonCut*ppi0ProtonBaryonCut;
        allGen_vanHove=allGeneralCutsPassed*pVanHove;



	// SOME MORE BOOLS FOR THE DECK PLOTS. 
	for ( int iMass=0; iMass < num_massBins; ++iMass){
		if ( mandelstam_teta > tMin && mandelstam_teta < tMax && mEllipse_pre_tAll) {
			if ( locPi0Eta > mMin+iMass*mStep && locPi0Eta < mMin+(iMass+1)*mStep ) {
				passMassBin_tetaIntegrated[iMass] = true;
			}
			else { passMassBin_tetaIntegrated[iMass] = false; } 
		}
		else { passMassBin_tetaIntegrated[iMass] = false; } 

		if ( mandelstam_tpi0 > tMin && mandelstam_tpi0 < tMax && mEllipse_pre_tAll) {
			if ( locPi0Eta > mMin+iMass*mStep && locPi0Eta < mMin+(iMass+1)*mStep ) {
				passMassBin_tpi0Integrated[iMass] = true;
			}
			else { passMassBin_tpi0Integrated[iMass] = false; } 
		}
		else { passMassBin_tpi0Integrated[iMass] = false; } 
	}

		// Full rect
	//withinBox(inBox, inBox_noOtherCuts, mEllipse_pre, locPi0Mass_Kin,locEtaMass_Kin,ellipseX-ellipseXr, ellipseX+ellipseXr, ellipseY-0.03, ellipseY+0.03, skipX-0.005, 0.030);
		// Diff rect size
		// 	1st: half size
		// 	2nd: half size closer
	//withinBox(inBox, inBox_noOtherCuts, mEllipse_pre, locPi0Mass_Kin,locEtaMass_Kin, 0.115, 0.155, 0.49, 0.59, 0.01, 0.02);
	//withinBox(inBox, inBox_noOtherCuts, mEllipse_pre, locPi0Mass_Kin,locEtaMass_Kin, 0.12, 0.15, 0.51, 0.57, 0.01, 0.035);
		// Using meas values, full rect
		//  ETA
		// par3 0.540383
		// par4 0.0233172
		//  PI0
		// par3 0.134285
		// par4 0.00769639
	withinBox(inBox, inBox_noOtherCuts, mEllipse_pre, locPi0Mass_Kin,locEtaMass_Kin, 0.135881-2*0.0076, 0.135881+2*0.0076, 0.548625-2*0.0191, 0.548625+2*0.0191, 0.01, 0.02);
	if ( inBox_noOtherCuts[10] ) { weightBS = 0.25; weightB = -1.25; }
	else if ( inBox_noOtherCuts[11] ) { weightBS = -0.5; weightB = 1.5; }
	else if ( inBox_noOtherCuts[12] ) { weightBS = -0.5; weightB = 1.5; }
	else if ( inBox_noOtherCuts[4] ) { weightBS = 1; weightB = 0; }
	else { weightBS=0; weightB=1; }

        //weightAS_BS = weightBS;//weightAS*weightBS;
        weightAS_B = weightAS*weightB;

	++count_combos;
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
	if(!pMPi0P14) { ++count_MPi0P14; dHist_Cuts->Fill(cutNames[15],1);}
	if(mMandelstamT_mBeamE8GeVPlus){ ++count_seanResTest; dHist_Cuts->Fill(cutNames[16],1);}


        if(showOutput) {cout << "Start Filling histVals and histCuts" << endl;}


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

        // Cut out this combo from the output tree.
        //if (!allGeneralCutsPassed) { dComboWrapper->Set_IsComboCut(true); continue; }

        // METHOD 1: We want to remove the elliptical cut (which doesn't affect our histograms since we filled then already) so that we have a flat tree that we can pass to another program
        // 	that will do the bkg subtraction. 
        // METHOD 2: Here we will apply a loose elliptical cut so that the data that we pass into SPlot has a better form that we can fit to
        //
        //

        group_PB.fillHistograms_Map( usingPB );
        group_12B_1234B.fillHistograms_pairMap( using12B_1234B );
        group_34B_1234B.fillHistograms_pairMap( using34B_1234B );
        group_12B.fillHistograms_Map( using12B );
        group_34B.fillHistograms_Map( using34B );
        group_12B_34B.fillHistograms_pairMap( using12B_34B );
        group_1234B.fillHistograms_Map( using1234B );
        group_PhNB.fillHistograms_vectorMap(beingUsedNeutralIds);
        group_12PhNB.fillHistograms_vectorMap( using12PhNBs );
        group_34PhNB.fillHistograms_vectorMap( using34PhNBs );
        group_pairFCAL.fillHistograms_vectorMap(usedPairIdsFCAL);
        group_pairBCAL.fillHistograms_vectorMap(usedPairIdsBCAL);
        group_1234BP.fillHistograms_Map( usingCombo );
        group_12PB.fillHistograms_Map(using12PB);       
        group_34PB.fillHistograms_Map(using34PB);       
        group_B.fillHistograms_Map(usingB);
        group_1234B_12PB.fillHistograms_pairMap(using1234B_12PB);
        group_1234B_34PB.fillHistograms_pairMap(using1234B_34PB);

	if(showOutput){ cout << "Filling histogram's uniqueness elements" << endl; }

	if ( selectDetector == "FCAL" ) { detectorCut=pEtaInFCAL; }
	else if ( selectDetector == "BCAL" ) { detectorCut=pEtaInBCAL; }
	else if ( selectDetector == "SPLIT" ) { detectorCut=pEtaInSplit; }
	else { detectorCut=true; }

        /****************************************** DO NOT CUT - FILL THROWN TOPLOGY (IF DESIRED) ******************************************/
						///     Obviously dont use this section if you are going to cut....
	// Fill histogram of thrown topologies
	//if (showThrownTopology){
	//	if(dHistInvariantMass_ThrownTopology.find(locThrownTopology) != dHistInvariantMass_ThrownTopology.end()) {
	//		dHistInvariantMass_ThrownTopology[locThrownTopology]->Fill(locPi0EtaP4_Kin.M());
	//		dHistThrownTopologies->Fill(locThrownTopology.Data(),1);
	//	}
	//}

	/******************************************* CUT ON THE COMBINATION *********************************************************/

        //if (!mEllipse_pre_tAll || !detectorCut) {
        if (!mEllipse_pre || !detectorCut) {
        //if (!mEllipseLooseUEChiSq_pre || !detectorCut)
	//if (!mMandelstamT_mBeamE8GeVPlus){
	    if (showOutput) { cout << "Did not pass cut, moving on.... " << endl; }  
            dComboWrapper->Set_IsComboCut(true); continue; 
        }

        else { 
		if (showOutput) { cout << "Passed cut, continuing.... " << endl; }  
	    	if (used1234B.find(using1234B)==used1234B.end()){
            	    used1234B.insert(using1234B);
            	    isNotRepeated_pi0eta=true;
            	}
            	else { isNotRepeated_pi0eta=false; } 

	    	if (used12B_1234B.find(using12B_1234B)==used12B_1234B.end()){
            	    used12B_1234B.insert(using12B_1234B);
            	    isNotRepeated_pi0_pi0eta=true;
            	}
            	else { isNotRepeated_pi0_pi0eta=false; } 

	    	if (used34B_1234B.find(using34B_1234B)==used34B_1234B.end()){
            	    used34B_1234B.insert(using34B_1234B);
            	    isNotRepeated_eta_pi0eta=true;
            	}
            	else { isNotRepeated_eta_pi0eta=false; } 

	    	if (used12B.find(using12B)==used12B.end()){
            	    used12B.insert(using12B);
            	    isNotRepeated_pi0=true;
            	}
            	else { isNotRepeated_pi0=false; } 

	    	if (used34B.find(using34B)==used34B.end()){
            	    used34B.insert(using34B);
            	    isNotRepeated_eta=true;
            	}
            	else { isNotRepeated_eta=false; } 

	    	if (used123B.find(using12PhNBs[0])==used123B.end()){
            	    used123B.insert(using12PhNBs[0]);
            	    isNotRepeated_pi0g1=true;
            	}
            	else { isNotRepeated_pi0g1=false; } 
	    	if (used124B.find(using34PhNBs[1])==used124B.end()){
            	    used124B.insert(using34PhNBs[1]);
            	    isNotRepeated_pi0g2=true;
            	}
            	else { isNotRepeated_pi0g2=false; } 


            	if ( uniqueSpectroscopicID.find(using34B)==uniqueSpectroscopicID.end() ){
            		uniqueSpectroscopicID.insert(using34B);
			// might be able to get rid of the use of both a set (to find) and a map to track. I think map would implement a find key also
			// problem for another time
			map_uniqueSpectroscopicID[using34B] = loc_i;
            	}
    			
		
		paddedCombo=to_string(map_uniqueSpectroscopicID[using34B]);
		if(showOutput){ cout << "Before paddedCombo: " << paddedCombo << endl; }
		digitsInCombo=paddedCombo.length();
		if(showOutput){ cout << "digitsInCombo: " << digitsInCombo << endl; }
		if(showOutput){ cout << "maxDigitsInCombo: " << maxDigitsInCombo << endl; }
		if(showOutput){ cout << "maxDigitsInCombo-digitsInCombo: " << maxDigitsInCombo-digitsInCombo << endl; }
		for ( int idigit=0; idigit < maxDigitsInCombo-digitsInCombo; ++idigit){
			paddedCombo = "0" + paddedCombo;
		}
		if(showOutput){ cout << "Final paddedCombo: " << paddedCombo << endl; }
		if(showOutput){ cout << "Final paddedEvent: " << paddedEvent << endl; }
		if(showOutput){ cout << "Final paddedRun: " << paddedRun << endl; }
		string s_spectroscopicComboID = "1"+paddedRun+paddedEvent+paddedCombo;
		if(showOutput){ cout << "1+paddedRun+paddedEvent+paddedCombo: " << s_spectroscopicComboID << endl; }
		spectroscopicComboID = (ULong64_t)stoull(s_spectroscopicComboID); 
		//digitsInSpectroscopicComboID = (int)log10(spectroscopicComboID);
		if(showOutput){ cout << "spectroscopicComboID: " << spectroscopicComboID << endl; }

        }

	if (showOutput){ cout << "Calculated uniqueness booleans" << endl; } 




        /****************************************** CUT - FILL THROWN TOPLOGY (IF DESIRED) ******************************************/
	// Fill histogram of thrown topologies
	if(showThrownTopology){
		if(dHistInvariantMass_ThrownTopology.find(locThrownTopology) != dHistInvariantMass_ThrownTopology.end()) {
			dHistInvariantMass_ThrownTopology[locThrownTopology]->Fill(locPi0EtaP4_Kin.M());
			dHistThrownTopologies->Fill(locThrownTopology.Data(),1);
		}
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

        dFlatTreeInterface->Fill_Fundamental<Double_t>("AccWeight", weightAS);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("weightASBS", weightAS_BS);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("uniqueComboID", uniqueComboID);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("finalStateComboID", finalStateComboID);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("chiSq", locChiSqKinFit);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("unusedEnergy", locUnusedEnergy);

	if(showOutput){ cout << "Filled some fundamental branches" << endl; } 

        ++uniqueComboID;
        if (is_pi0eta){
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0g1", massGammaPi0[0]);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0g2", massGammaPi0[1]);
        	dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0g1",isNotRepeated_pi0g1);
        	dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0g2",isNotRepeated_pi0g2);

        	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0", locPi0Mass_Kin);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("Meta", locEtaMass_Kin);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0eta", locPi0Eta_Kin);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0_meas", locPi0Mass);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("Meta_meas", locEtaMass);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0eta_meas", locPi0Eta);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("mandelstam_teta_meas", mandelstam_teta);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("mandelstam_teta", mandelstam_teta_Kin);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("mandelstam_tpi0_meas", mandelstam_tpi0);
        	dFlatTreeInterface->Fill_Fundamental<Double_t>("mandelstam_tpi0", mandelstam_tpi0_Kin);
        }
        else{
            dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0", locPi0Mass_Kin);
            dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0", locEtaMass_Kin);
            dFlatTreeInterface->Fill_Fundamental<Double_t>("Mpi0pi0", locPi0Eta_Kin);
        }
        dFlatTreeInterface->Fill_Fundamental<Double_t>("mandelstam_tp", mandelstam_tp);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("ptGT1", ptGT1);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("ptLT05", ptLT05);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("ptGT05LT1", ptGT05LT1);
	if(showOutput){ cout << "Filled some more fundamental branches" << endl; } 
	
	// some variables to track uniqueness
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0",isNotRepeated_pi0);
	if(showOutput){ cout << "Filled even more fundamental branches" << endl; } 
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_eta",isNotRepeated_eta);
	if(showOutput){ cout << "Filled even more fundamental branches" << endl; } 
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0eta",isNotRepeated_pi0eta);
	if(showOutput){ cout << "Filled even more fundamental branches" << endl; } 
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_eta_pi0eta",isNotRepeated_eta_pi0eta);
        dFlatTreeInterface->Fill_Fundamental<Bool_t>("isNotRepeated_pi0_pi0eta",isNotRepeated_pi0_pi0eta);
        dFlatTreeInterface->Fill_Fundamental<ULong64_t>("spectroscopicComboID",spectroscopicComboID);
	if(showOutput){ cout << "Filled even more fundamental branches" << endl; } 

	// some variables to show where the pi0/eta detected
        dFlatTreeInterface->Fill_Fundamental<Double_t>("pi0DetectedIn", pi0DetectedIn); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("etaDetectedIn", etaDetectedIn); //fundamental = char, int, float, double, etc.
	
        // Introduce some angles to use in the phase space distance calculation
        dFlatTreeInterface->Fill_Fundamental<Double_t>("cosTheta_X_cm", cosTheta_pi0eta_CM); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("cosTheta_X_cm_meas", cosTheta_pi0eta_CM_meas); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("phi_X_cm", phi_pi0eta_CM); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("cosTheta_eta_gj", cosTheta_eta_GJ); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("cosTheta_eta_gj_meas", cosTheta_eta_GJ_meas); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("phi_eta_gj", phi_eta_GJ); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("phi_eta_gj_meas", phi_eta_GJ_meas); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("cosThetaHighestEphotonIneta_gj", cosTheta_largestEinEta_GJ); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("cosThetaHighestEphotonInpi0_cm", cosTheta_largestEinPi0_CM); //fundamental = char, int, float, double, etc.
        dFlatTreeInterface->Fill_Fundamental<Double_t>("vanHove_x",vanHove_x);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("vanHove_y",vanHove_y);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("vanHove_omega",omega);
        // If we were to do Q-Values for the pi0pi0 system we will probably do the same thing and use one of the pi0 as the discriminator varible and check against the other
        dFlatTreeInterface->Fill_Fundamental<Double_t>("pi0_energy", locPi0E_Kin );

	if(showOutput){ cout << "Filled fundamental branches" <<endl; }
        //FILL FLAT TREE
        //dComboWrapper->Set_ComboIndex(loc_i);  // Combo succeeded // this might be redundant since this is after the continue command in the cut condition above, combo must have succeeded already. But Elton has it! 
        Fill_FlatTree(); //for the active combo
    } // end of combo loop
    ++eventIdx;
    if(showOutput){cout << "\n\n **************** Finishing the combo loop ***************\n**********************************************************\n" << endl;}


    //FILL HISTOGRAMS: Num combos / events surviving actions
    Fill_NumCombosSurvivedHists();
    if(showOutput){ cout << "Fillied NumComboSurvivedHists" << endl; } 

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
    if(showOutput) { cout << "Looping over through the comobs to check passed or not" << endl; } 
    for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
        //Set branch array indices for combo and all combo particles
        dComboWrapper->Set_ComboIndex(loc_i);
        // Is used to indicate when combos have been cut
        if(dComboWrapper->Get_IsComboCut()){
    	    if(showOutput) { cout << "Combo did not pass cuts" << endl; } 
            continue;
	}
        locIsEventCut = false; // At least one combo succeeded                                                     
        if(showOutput) { cout << "Combo passed cuts!"  << endl; } 
        break;
    }
    if(!locIsEventCut && dOutputTreeFileName != ""){ 
	    if (showOutput) {cout<<"Filling outputt tree" << endl; }
	    Fill_OutputTree(); 
    }

    //}//closes the //if(itersToRun) condition
    return kTRUE; // this return should close the process loop to return false as the kTrue as the output.
}// end of process loop

void DSelector_ver20::Finalize(void)
{
    //group_PhNB.drawHistograms();
    //group_PB.drawHistograms();
    //group_12B_1234B.drawHistograms();
    //group_34B_1234B.drawHistograms();
    //group_12B.drawHistograms();
    //group_34B.drawHistograms();

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

    if(true){cout << "Count combos: " << std::to_string(count_combos) << endl;} 
    if(true){cout << "Count events: " << std::to_string(count_events) << endl;} 
    if(true){cout << "Count tot events: " << std::to_string(count_totEvents) << endl;} 
    if(true){cout << "Percent passed ShowerQuality: " << std::to_string((double)count_ShowerQuality/count_combos)<<endl; }
    if(true){cout << "Percent passed BeamE8GeVPlus: " << std::to_string((double)count_BeamE8GeVPlus/count_combos)<<endl; }
    if(true){cout << "Percent passed UnusedEnergy: " << std::to_string((double)count_UnusedEnergy/count_combos)<<endl;}
    if(true){cout << "Percent passed CLKinFit: " << std::to_string((double)count_ChiSq/count_combos)<<endl;}
    if(true){cout << "Percent passed DeltaTRF: " << std::to_string((double)count_DeltaTRF/count_combos)<<endl;}
    if(true){cout << "Percent passed dij3pass: " << std::to_string((double)count_dij3pass/count_combos)<<endl;}
    if(true){cout << "Percent passed PhotonE: " << std::to_string((double)count_PhotonE/count_combos)<<endl;}
    if(true){cout << "Percent passed PhotonTheta: " << std::to_string((double)count_PhotonTheta/count_combos)<<endl;}
    if(true){cout << "Percent passed MagP3Proton: " << std::to_string((double)count_MagP3Proton/count_combos)<<endl;}
    if(true){cout << "Percent passed zCutmin: " << std::to_string((double)count_zCutmin/count_combos)<<endl;}
    if(true){cout << "Percent passed RProton: " << std::to_string((double)count_RProton/count_combos)<<endl;}
    if(true){cout << "Percent passed MissingMassSquared: " << std::to_string((double)count_MissingMassSquared/count_combos)<<endl;}
    if(true){cout << "Percent passed dEdxCDCProton: " << std::to_string((double)count_dEdxCDCProton/count_combos)<<endl;}
    if(true){cout << "Percent passed insideEllipse: " << std::to_string((double)count_insideEllipse/count_combos)<<endl;}
    if(true){cout << "Percent passed allGeneralCutsPassed: " << std::to_string((double)count_allGeneralCutsPassed/count_combos)<<endl;}
    if(true){cout << "Percent passed mMPi0P14: " << std::to_string((double)count_MPi0P14/count_combos)<<endl;}
    if(true){cout << "Percent pass resolutionTest: " << std::to_string((double)count_seanResTest/count_combos) << endl; }
    if(true){cout << "Percent correct topology: " << std::to_string((double)count_correctTopology/count_events) << endl; }
    // we can only use the below code when we are using setupTest.sh. DOesnt work with proof since it will probably try to do this for every thread...
    //dHist_Cuts->SetStats(0);
    //dHist_Cuts->SetCanExtend(TH1::kAllAxes);
    //dHist_Cuts->LabelsDeflate("X");
    //dHist_Cuts->LabelsOption("v");

    //CALL THIS LAST
    if(showOutput){cout << "finalizing!" << endl;}
    DSelector::Finalize(); //Saves results to the output file
};

