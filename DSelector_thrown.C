#include "DSelector_thrown.h"
#include "TRandom.h"
string polarization="deg000";
string tag="_flat_gen";

void DSelector_thrown::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "v20_pi0eta_Thrown.root"; //"" for none
	//USERS: SET OUTPUT TREE FILES/NAMES //e.g. binning into separate files for AmpTools
        dOutputTreeFileNameMap["selected"] = polarization+tag+"_allEtaDecays_DSelector.root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["selected_tLT06"] = "thrownNotAmptoolsReady_a0a2_tLT06.root"; //key is user-defined, value is output file name
	//dOutputTreeFileNameMap["selected_tGT05LT1"] = "thrownNotAmptoolsReady_a0a2_tGT05LT1.root"; //key is user-defined, value is output file name

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	dPreviousRunNumber = 0;
        countThrownEvents = new TH1F("count thrown events", "Cuts=noCut;counts", 1, 0, 1);
        dHist_prodPlanePS_000 = new TH1F("prodPlanePS_000", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        dHist_prodPlanePS_045 = new TH1F("prodPlanePS_045", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        dHist_prodPlanePS_090 = new TH1F("prodPlanePS_090", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        dHist_prodPlanePS_135 = new TH1F("prodPlanePS_135", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        dHist_prodPlanePS_AMO = new TH1F("prodPlanePS_AM0", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //dHist_prodPlanePS_000_rejSamp = new TH1F("prodPlanePS_000_rejSamp", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //dHist_prodPlanePS_045_rejSamp = new TH1F("prodPlanePS_045_rejSamp", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //dHist_prodPlanePS_090_rejSamp = new TH1F("prodPlanePS_090_rejSamp", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //dHist_prodPlanePS_135_rejSamp = new TH1F("prodPlanePS_135_rejSamp", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        //dHist_prodPlanePS_AMO_rejSamp = new TH1F("prodPlanePS_AM0_rejSamp", "Cuts=ptEtaBeamAsym;#phi; Entries / 9 degrees", 40, -180, 180);
        dHist_BeamAngle = new TH1F("BeamAngle", "Beam Angle with no cuts applied;Beam Angle (GeV);Events / 2 Degree", 180,-180,180);
        dHist_SelectedBeamAngle = new TH1F("SelectedBeamAngle", "Beam Angle with no cuts applied;Beam Angle (GeV);Events / 2 Degree", 180,-180,180);
	dHist_PID = new TH1I("PID","",50,0,50);
	dHist_NumThrown = new TH1I("NumThrown","",10,0,10);
	dHist_beamE = new TH1F("beamE","Beam Energy", 100,0,15);
	dHist_beamECut = new TH1F("beamECut","Beam Energy", 100,0,15);
	mandelstam_tpAll = new TH1F("mandelstam_tpAll","tprime no cut",100,0,1.1);
	//mandelstam_tpAll = new TH1F("mandelstam_tpAll","tprime no cut",100,0,6);
	mandelstam_tpAll_selected = new TH1F("mandelstam_tpAll_selected","tprime selected",100,0,6);
	mandelstam_tAll = new TH1F("mandelstam_tAll","tprime",100,0,6);
	mandelstam_tpLT1 = new TH1F("mandelstam_tpLT1","tprime<1",100,0,6);
	mandelstam_tpLT06 = new TH1F("mandelstam_tpLT06","tprime<0.6",100,0,6);
	mandelstam_tpGT05LT1 = new TH1F("mandelstam_tpGT05LT1","0.5<tprime<1",100,0,6);
	dHist_cosThetaVsMass_tpAll = new TH2F("cosThetaVsMass_tpAll","cosTheta vs Mass",175,0,3.5,100,-1,1);
	dHist_cosThetaVsMass_tpLT1 = new TH2F("cosThetaVsMass_tpLT1","cosTheta vs Mass",175,0,3.5,100,-1,1);
	dHist_cosThetaVsMass_tpLT06 = new TH2F("cosThetaVsMass_tpLT06","cosTheta vs Mass",175,0,3.5,100,-1,1);
	dHist_cosThetaVsMass_tpGT05LT1 = new TH2F("cosThetaVsMasstpGT05LT1","cosTheta vs Mass",175,0,3.5,100,-1,1);
	dHist_phiVsMass = new TH2F("phiVsMass","phi vs mass", 150,0,3.5,60,-180,180);
	dHist_phi = new TH1F("phi","phi GJ",60,-180,180);
	dHist_cosTheta = new TH1F("cosTheta","cosTheta GJ",100,-1,1);

	dHist_numEventsOnePi0OneEta = new TH1F( "onePi0oneEta", "", 2,0,1);
	// included the -1 for underflow maybe?
        dHist_genCounts_eta_tLT05 = new TH1F("tetaVsMpi0eta_genCounts_tLT05", "genCounts", numHists+1, -1, numHists);
        dHist_genCounts_pi0_tLT05 = new TH1F("tpi0VsMpi0eta_genCounts_tLT05", "genCounts", numHists+1, -1, numHists);
        dHist_genCounts_eta_tGT05LT1 = new TH1F("tetaVsMpi0eta_genCounts_tGT05LT1", "genCounts", numHists+1, -1, numHists);
        dHist_genCounts_pi0_tGT05LT1 = new TH1F("tpi0VsMpi0eta_genCounts_tGT05LT1", "genCounts", numHists+1, -1, numHists);
        dHist_genCounts_eta_tGT1 = new TH1F("tetaVsMpi0eta_genCounts_tGT1", "genCounts", numHists+1, -1, numHists);
        dHist_genCounts_pi0_tGT1 = new TH1F("tpi0VsMpi0eta_genCounts_tGT1", "genCounts", numHists+1, -1, numHists);
        dHist_genCounts_eta_tAll = new TH1F("tetaVsMpi0eta_genCounts_tAll", "genCounts", numHists+1, -1, numHists);
        dHist_genCounts_pi0_tAll = new TH1F("tpi0VsMpi0eta_genCounts_tAll", "genCounts", numHists+1, -1, numHists);
	
	for (int beamE=0; beamE<12; ++beamE){
		dHist_pi0eta1DBeam[beamE] = new TH1F(("pi0eta1DbeamE"+std::to_string(beamE)).c_str(),"M(pi0eta)",150,0,3);
	}
	dHist_pi0eta1D = new TH1F("pi0eta1D","M(pi0eta)",350,0,3.5);
	dHist_phi8GeVPlus = new TH1F("phi8GeVPlus","phi GJ",60,-180,180);
	dHist_cosTheta8GeVPlus = new TH1F("cosTheta8GeVPlus","cosTheta GJ",60,-1,1);

        dHist_Mpi0eta = new TH1F("Mpi0eta", "Mpi0eta;M(#pi^{0}#eta) (GeV);Events / 0.01 GeV", 60,1.6,2.8);
	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want
	//
	ievent=0;
	maxevent=100;

        dTreeInterface->Create_Branch_Fundamental<Double_t>("mandelstam_tp"); //fundamental = char, int, float, double, etc.
        dTreeInterface->Create_Branch_Fundamental<Double_t>("Mpi0eta"); //fundamental = char, int, float, double, etc.
        dTreeInterface->Create_Branch_Fundamental<Double_t>("cosTheta_eta_gj"); //fundamental = char, int, float, double, etc.
}

Bool_t DSelector_thrown::Process(Long64_t locEntry)
{
	//if (ievent > maxevent) { return true; } 
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
	//
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//
	//

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	cout << "RunNumber = " << locRunNumber << endl;
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
        	hasPolarizationAngle = dAnalysisUtilities.Get_PolarizationAngle(locRunNumber, locPolarizationAngle);
		dPreviousRunNumber = locRunNumber;
        	cout << "Getting beam polarization and filling used runs" << endl;
        	if (hasPolarizationAngle) {
        	    dHist_BeamAngle->Fill(locPolarizationAngle);
        	}
		else {
		    dHist_BeamAngle->Fill(-1);
		}
	}

	countThrownEvents->Fill(1);

    	keepPolarization=false;
    	if (polarization=="degALL"){ keepPolarization=true; }
    	else if (locPolarizationAngle==0 && polarization=="deg000") { keepPolarization=true; }
    	else if (locPolarizationAngle==45 && polarization=="deg045") { keepPolarization=true; }
    	else if (locPolarizationAngle==90 && polarization=="deg090") { keepPolarization=true; }
    	else if (locPolarizationAngle==135 && polarization=="deg135") { keepPolarization=true; }
    	else if (!hasPolarizationAngle && polarization=="degAMO") { keepPolarization=true; }
    	else {  cout << "FUDGE! THERE IS AN UNEXPECTED OUTCOME FROM THE POLARIZATION CHECK" << endl; } 

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/
	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/******************************************* LOOP OVER THROWN DATA ***************************************/

	//Thrown beam: just use directly
	bool beamEisFinite = 0;
	double locBeamEnergyUsedForBinning = 0.0;
	if(dThrownBeam != NULL)
		locBeamEnergyUsedForBinning = dThrownBeam->Get_P4().E();
	if(TMath::Finite(locBeamEnergyUsedForBinning)==1){
		++beamEfin;
		beamEisFinite = 1;	
	}
	else { ++beamEinf; }	

        double cosTheta_eta_GJ; 
	double locPi0EtaMass;
	double mandelstam_t;

	std::vector<TLorentzVector> allP4;
	std::vector<int> parentArray;
	std::vector<int> pids;
	TLorentzVector locEtaP4;
	TLorentzVector locPi0P4;
	TLorentzVector locProtonP4;
	TLorentzVector locTargetP4 = {0,0,0,0.938};
	TLorentzVector locBeamP4=dThrownBeam->Get_P4();

	//Loop over throwns
	// we use is_in as a way to select the events we want to show instead of dumping everything
	const bool is_in = showOutput.find(eventIdx) != showOutput.end();
	if (is_in){
		cout << "########    EventIdx: " << eventIdx << "    #############" << endl;
	}

	int countPrimaryEta=0;
	int countPrimaryPi0=0;
	bool notCorrectTopology=0;
	int locNumThrown = Get_NumThrown();
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{	
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
		Particle_t locPID = dThrownWrapper->Get_PID();
		Int_t locParentPID = dThrownWrapper->Get_ParentIndex();
		allP4.push_back(dThrownWrapper->Get_P4());

		if (locParentPID==-1 && locPID==7) { ++countPrimaryPi0; locPi0P4 = dThrownWrapper->Get_P4(); } 
		if (locParentPID==-1 && locPID==17) { ++countPrimaryEta; locEtaP4 = dThrownWrapper->Get_P4(); } 
		if(locParentPID==-1 && locPID==14) { locProtonP4 = dThrownWrapper->Get_P4(); } 

		//Some verbal confirmation checking
		if (is_in){
			cout << "Thrown " << loc_i << " with (PID,ParentArrayLoc) = (" << locPID << " ," << locParentPID << ")"  << endl;
		}
		dHist_PID->Fill(locPID);	

		parentArray.push_back(locParentPID);
		pids.push_back(locPID);
	}
	if ( countPrimaryPi0!=1 || countPrimaryEta!=1 ) {
		// Interesting, that there is an event that starts with one eta -> 3pi0 -> 6 gamma
		cout << "SHOOT I HAVE TO SHOULD BREAK PROGRAM TO REMIND MYSELF THAT ASKING FOR A 7 AND 17 DOESNT GUARANTEE IT" << endl; 
		cout << " These are the ids and their parents: " << endl;
		for ( int idx=0; idx<(int)pids.size(); ++idx ){
			cout << "  " << pids[idx] << ", " << parentArray[idx] << endl;
		}
		notCorrectTopology=true;
		//exit(0);
	}

	if ( !notCorrectTopology ) { // if it has 1 pi0 and eta parent then we see if the daughters are both photons
		std::vector<int> parents;
		findParents(parentArray,parents);
		cout << "\nTHESE ARE THE PARENTS" << endl;
		for ( auto parent=0; parent < (int)parents.size(); ++parent){
			cout << parents[parent] << endl;
		}

		int pi0ToNGamma=0;
		int etaToNGamma=0;
		bool correctFinalState=false;
		for (auto parent : parents){ // each one of these parents are already a primary particle
			std::vector<int> daughters;
			cout << "Parent: " << parent << " which has PID=" << pids[parent] << " has children:" << endl;
			findDaughters( parentArray, daughters, parent );
			for ( auto daughter=0; daughter < (int)daughters.size(); ++daughter) {
				if ( pids[parent]==7 && pids[daughters[daughter]] == 1 && daughters.size()==2){ ++pi0ToNGamma; }
				if ( pids[parent]==17 && pids[daughters[daughter]] == 1 && daughters.size()==2){ ++etaToNGamma; }
				//cout << "-  " << daughters[daughter] << " which has PID=" << pids[daughters[daughter]] << endl;
				//
				// THESE CODE BELOW WILL GET THE SECONDARY DAUGHTERS
				//std::vector<int> secDaughters;
				//findDaughters( parentArray, secDaughters, daughters[daughter] ); 
				//for ( int secDaughter=0; secDaughter < secDaughters.size(); ++secDaughter) {
				//	cout << "--" << secDaughters[secDaughter] << endl;
				//}
			}
		}
		if ( pi0ToNGamma==2 && etaToNGamma==2) {
			correctFinalState=true;
                        cout << "(Correct topology)" << endl;
                }
                else {
                        cout << "(Incorrect topology) Some histograms will still be filled without the correct final state! FYI" << endl;
                }
		cout << "Found all the daughters" << endl;

	        // check for etas in idxInitial
	        Int_t eta = 17;
	        Int_t pi0 = 7;
	        Int_t proton =14;

	        double lowE = 0;
	        double uppE = 1;
	        for (int beamE=0; beamE<12; ++beamE){
	        	//cout << "lowE, uppE: "<<lowE<<", "<<uppE<<endl;
	        	pBeamE[beamE] = lowE<locBeamP4.E() && locBeamP4.E()<uppE;
	        	lowE+=1;
	        	uppE+=1;
	        }
	        cout << "Calculated EBeam bins" << endl;

                //double prodPlanePhi = dAnalysisUtilities.Calc_ProdPlanePhi_Pseudoscalar(locBeamP4.E(), Proton, locEtaP4);
	        double prodPlanePhi = locEtaP4.Phi()*180/3.14159;
	        
	        TLorentzVector locPi0EtaP4 = locPi0P4+locEtaP4;

                TLorentzVector cm_vec = locBeamP4+locTargetP4;
                TLorentzVector locPi0EtaP4_cm = locPi0EtaP4;
                TLorentzVector locPi0P4_cm = locPi0P4;
                TLorentzVector locEtaP4_cm = locEtaP4;
                TLorentzVector locBeamP4_cm = locBeamP4;
                TLorentzVector locProtonP4_cm = locProtonP4;
                locPi0EtaP4_cm.Boost(-cm_vec.BoostVector());
                locPi0P4_cm.Boost(-cm_vec.BoostVector());
                locEtaP4_cm.Boost(-cm_vec.BoostVector());
                locBeamP4_cm.Boost(-cm_vec.BoostVector());
                locProtonP4_cm.Boost(-cm_vec.BoostVector());

	        TLorentzVector locPi0Eta_gj = locPi0EtaP4_cm;	
	        TLorentzVector locPi0P4_gj = locPi0P4_cm;
	        TLorentzVector locEtaP4_gj = locEtaP4_cm;
	        TLorentzVector locBeamP4_gj = locBeamP4_cm;
	        TLorentzVector locProtonP4_gj = locProtonP4_cm;
	        locPi0Eta_gj.Boost(-locPi0EtaP4_cm.BoostVector());
	        locPi0P4_gj.Boost(-locPi0EtaP4_cm.BoostVector());
	        locEtaP4_gj.Boost(-locPi0EtaP4_cm.BoostVector());
	        locBeamP4_gj.Boost(-locPi0EtaP4_cm.BoostVector());
	        locProtonP4_gj.Boost(-locPi0EtaP4_cm.BoostVector());

	        double radToDeg = 57.3;
                TVector3 locPi0P4_gj_unit = locPi0P4_gj.Vect().Unit();
                TVector3 locEtaP4_gj_unit = locEtaP4_gj.Vect().Unit();
                TVector3 locPi0EtaP4_gj_unit = locPi0Eta_gj.Vect().Unit();
                // Calculate cosTheta, phi in maybe the GJ axes.
                // since we already defined the x,y,z as TVector3 we don't have to do it again.
                TVector3 z = locBeamP4_gj.Vect().Unit();
                // this y should be the normal of the production plane. If we do a boost in a direction in the production plane the perp direction doesn't change. We could use the beam and the recoiled proton to define the
                // production plane in this new frame. Let us define it in the CM frame. 
                TVector3 y = locPi0EtaP4_cm.Vect().Cross(locBeamP4_cm.Vect()).Unit();
                TVector3 x = y.Cross(z).Unit();

	        TVector3 angles_pi0;
	        TVector3 angles_eta;
                angles_pi0.SetXYZ ( locPi0P4_gj_unit.Dot(x), locPi0P4_gj_unit.Dot(y), locPi0P4_gj_unit.Dot(z) );
                angles_eta.SetXYZ ( locEtaP4_gj_unit.Dot(x), locEtaP4_gj_unit.Dot(y), locEtaP4_gj_unit.Dot(z) );

                double cosTheta_pi0_GJ = angles_pi0.CosTheta();
                cosTheta_eta_GJ= angles_eta.CosTheta();
                double phi_pi0_GJ = angles_pi0.Phi()*radToDeg;
                double phi_eta_GJ = angles_eta.Phi()*radToDeg;

	        cout << "Calculated the kinematic angles" << endl;

                locPi0EtaMass = locPi0EtaP4.M();
                mandelstam_t = -(locProtonP4-locTargetP4).M2();
	        double mandelstam_abst = abs(mandelstam_t);
	        //double mandelstam_t0 = -((locProtonP4.M2()-locPi0EtaP4.M2()-locTargetP4.M2())/(2*(locBeamP4+locTargetP4).M())-(locBeamP4_cm-locPi0EtaP4_cm).M2());
	        //above formulation is wrong. the last term uses magnitue, not p4
	        double mandelstam_t0 = -(TMath::Power(-locPi0EtaP4.M2()/(2*(locBeamP4+locTargetP4).M()),2)-TMath::Power(locBeamP4_cm.Vect().Mag()-locPi0EtaP4_cm.Vect().Mag(),2));
	        mandelstam_tp = mandelstam_t-mandelstam_t0;
	        
	        dHist_NumThrown->Fill(locNumThrown);
	        dHist_phiVsMass->Fill(locPi0EtaMass,phi_pi0_GJ);
	        dHist_phi->Fill(phi_pi0_GJ);
	        dHist_cosTheta->Fill(cosTheta_eta_GJ);
	        dHist_beamE->Fill(locBeamP4.E());
	        for (int beamE=0; beamE<12; ++beamE){
	        	if (pBeamE[beamE]){
	        		dHist_pi0eta1DBeam[beamE]->Fill(locPi0EtaMass);
	        	}
	        }
	        pBeamE8to9GeV=locBeamP4.E()>8 && locBeamP4.E()<9;	
	        cout << "Filled some histograms" << endl;

	        //if ( passCutsOneEtaPi0 && mandelstam_tp < 1) {
	        //}
                mandelstam_teta = -(locBeamP4-locEtaP4).M2();
                mandelstam_tpi0 = -(locBeamP4-locPi0P4).M2();
	        idx_t_eta = (int)( (mandelstam_teta-tMin)/tStep ); 
	        idx_t_pi0 = (int)( (mandelstam_tpi0-tMin)/tStep ); 
	        idx_m = (int)( (locPi0EtaP4.M()-mMin)/mStep );
	        if ( mandelstam_teta < tMin || locPi0EtaP4.M() < mMin || mandelstam_teta>tMax || locPi0EtaP4.M() >mMax ) {
	        	teta_genCounts = -1;
	        }
	        else {
	        	teta_genCounts = idx_t_eta+num_tBins*idx_m;
	        }
	        if ( mandelstam_tpi0 < tMin || locPi0EtaP4.M() < mMin || mandelstam_tpi0>tMax || locPi0EtaP4.M() >mMax ) {
	        	tpi0_genCounts = -1;
	        }
	        else {
	        	tpi0_genCounts = idx_t_pi0+num_tBins*idx_m;
	        }
	        cout << "Determined teta and tpi0 bins" << endl;


	        bool pBeamE8GeV = locBeamP4.E() > 8;
	        bool pBeamE8288 = 8.2 < locBeamP4.E() &&  locBeamP4.E() < 8.8;

	        if(correctFinalState*keepPolarization){
                        dHist_Mpi0eta->Fill(locPi0EtaMass);

	        	mandelstam_tpAll->Fill(mandelstam_tp);	
	        	mandelstam_tAll->Fill(mandelstam_abst);

	        	dHist_SelectedBeamAngle->Fill(locPolarizationAngle);

	        	dHist_beamECut->Fill(locBeamP4.E());
	        	dHist_cosThetaVsMass_tpAll->Fill(locPi0EtaMass,cosTheta_eta_GJ);
	        	dHist_genCounts_eta_tAll->Fill(teta_genCounts);
	        	dHist_genCounts_pi0_tAll->Fill(tpi0_genCounts);

	        	//Fill_OutputTree must be run with proof!
	        	dHist_pi0eta1D->Fill(locPi0EtaMass);
	        	dHist_phi8GeVPlus->Fill(phi_pi0_GJ);
	        	dHist_cosTheta8GeVPlus->Fill(cosTheta_eta_GJ);

	        	if (mandelstam_t<0.5){
	        		dHist_genCounts_eta_tLT05->Fill(teta_genCounts);
	        		dHist_genCounts_pi0_tLT05->Fill(tpi0_genCounts);
	        	}
	        	if (mandelstam_t>0.5 && mandelstam_t<1){
	        		dHist_genCounts_eta_tGT05LT1->Fill(teta_genCounts);
	        		dHist_genCounts_pi0_tGT05LT1->Fill(tpi0_genCounts);
	        	}
	        	if (mandelstam_t>1){
	        		dHist_genCounts_eta_tGT1->Fill(teta_genCounts);
	        		dHist_genCounts_pi0_tGT1->Fill(tpi0_genCounts);
	        	}

	        	if(mandelstam_t < 1){
	        		mandelstam_tpLT1->Fill(mandelstam_tp);
	        		dHist_cosThetaVsMass_tpLT1->Fill(locPi0EtaMass,cosTheta_eta_GJ);
	        	}

	        	mandelstam_tpAll_selected->Fill(mandelstam_tp);
	        	if ( locPolarizationAngle == 0 ) { 
	        		dHist_prodPlanePS_000->Fill(prodPlanePhi);
	        	}
	        	if ( locPolarizationAngle == 45 ) { 
	        		dHist_prodPlanePS_045->Fill(prodPlanePhi);
	        	}
	        	if ( locPolarizationAngle == 90 ) { 
	        		dHist_prodPlanePS_090->Fill(prodPlanePhi);
	        	}
	        	if ( locPolarizationAngle == 135 ) { 
	        		dHist_prodPlanePS_135->Fill(prodPlanePhi);
	        	}
	        	if ( !hasPolarizationAngle ) { 
	        		dHist_prodPlanePS_AMO->Fill(prodPlanePhi);
	        	}
                        dTreeInterface->Fill_Fundamental<Double_t>("mandelstam_tp", mandelstam_tp); //fundamental = char, int, float, double, etc.
                        dTreeInterface->Fill_Fundamental<Double_t>("cosTheta_eta_gj",cosTheta_eta_GJ); //fundamental = char, int, float, double, etc.
                        dTreeInterface->Fill_Fundamental<Double_t>("Mpi0eta", locPi0EtaMass); //fundamental = char, int, float, double, etc.
	                Fill_OutputTree();
	        } // if cuts not passed
	        Fill_OutputTree("selected"); //your user-defined key
        } // Cut for checking if pi0 and eta is there
	// If the following selections did not pass then this would have never executed and thus things are never filled

	++eventIdx;

	//OR Manually:
	//BEWARE: Do not expect the particles to be at the same array indices from one event to the next!!!!
	//Why? Because while your channel may be the same, the pions/kaons/etc. will decay differently each event.
	
	
	
	
	//cout << matchFOM << endl;

	//BRANCHES: https://halldweb.jlab.org/wiki/index.php/Analysis_TTreeFormat#TTree_Format:_Simulated_Data
/*
	Particle_t locThrown1PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[0]);
	TLorentzVector locThrown1P4 = *((TLorentzVector*)(*locP4Array)->At(0));
	cout << "Particle 1: " << locThrown1PID << ", " << locThrown1P4.Px() << ", " << locThrown1P4.Py() << ", " << locThrown1P4.Pz() << ", " << locThrown1P4.E() << endl;
	Particle_t locThrown2PID = PDGtoPType(((Int_t*)locPIDBranch->GetAddress())[1]);
	TLorentzVector locThrown2P4 = *((TLorentzVector*)(*locP4Array)->At(1));
	cout << "Particle 2: " << locThrown2PID << ", " << locThrown2P4.Px() << ", " << locThrown2P4.Py() << ", " << locThrown2P4.Pz() << ", " << locThrown2P4.E() << endl;
*/


	/******************************************* BIN THROWN DATA INTO SEPARATE TREES FOR AMPTOOLS ***************************************/

/*
	//THESE KEYS MUST BE DEFINED IN THE INIT SECTION (along with the output file names)
	if((locBeamEnergyUsedForBinning >= 8.0) && (locBeamEnergyUsedForBinning < 9.0))
		Fill_OutputTree("Bin1"); //your user-defined key
	else if((locBeamEnergyUsedForBinning >= 9.0) && (locBeamEnergyUsedForBinning < 10.0))
		Fill_OutputTree("Bin2"); //your user-defined key
	else if((locBeamEnergyUsedForBinning >= 10.0) && (locBeamEnergyUsedForBinning < 11.0))
		Fill_OutputTree("Bin3"); //your user-defined key
*/
	ievent++;
	return kTRUE;
}

void DSelector_thrown::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
