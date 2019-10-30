#include "DSelector_thrown_7_17_14.h"

void DSelector_thrown_7_17_14::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "v20_pi0eta_Thrown.root"; //"" for none
	//USERS: SET OUTPUT TREE FILES/NAMES //e.g. binning into separate files for AmpTools
	dOutputTreeFileNameMap["selected_tLT1"] = "thrownNotAmptoolsReady_a0a2_tLT1.root"; //key is user-defined, value is output file name
	dOutputTreeFileNameMap["selected_tLT06"] = "thrownNotAmptoolsReady_a0a2_tLT06.root"; //key is user-defined, value is output file name
	dOutputTreeFileNameMap["selected_tGT05LT1"] = "thrownNotAmptoolsReady_a0a2_tGT05LT1.root"; //key is user-defined, value is output file name

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	dPreviousRunNumber = 0;
	dHist_PID = new TH1I("PID","",50,0,50);
	dHist_NumThrown = new TH1I("NumThrown","",10,0,10);
	dHist_beamE = new TH1F("beamE","Beam Energy", 100,0,15);
	mandelstam_tpAll = new TH1F("mandelstam_tpAll","tprime",100,0,6);
	mandelstam_tpLT1 = new TH1F("mandelstam_tpLT1","tprime<1",100,0,6);
	mandelstam_tpLT06 = new TH1F("mandelstam_tpLT06","tprime<0.6",100,0,6);
	mandelstam_tpGT05LT1 = new TH1F("mandelstam_tpGT05LT1","0.5<tprime<1",100,0,6);
	dHist_cosThetaVsMass_tpAll = new TH2F("cosThetaVsMass_tpAll","cosTheta vs Mass",150,0,3.5,60,-1,1);
	dHist_cosThetaVsMass_tpLT1 = new TH2F("cosThetaVsMass_tpLT1","cosTheta vs Mass",150,0,3.5,60,-1,1);
	dHist_cosThetaVsMass_tpLT06 = new TH2F("cosThetaVsMass_tpLT06","cosTheta vs Mass",150,0,3.5,60,-1,1);
	dHist_cosThetaVsMass_tpGT05LT1 = new TH2F("cosThetaVsMasstpGT05LT1","cosTheta vs Mass",150,0,3.5,60,-1,1);
	dHist_phiVsMass = new TH2F("phiVsMass","phi vs mass", 150,0,3.5,60,-180,180);
	dHist_phi = new TH1F("phi","phi GJ",60,-180,180);
	dHist_cosTheta = new TH1F("cosTheta","cosTheta GJ",60,-1,1);
	
	for (int beamE=0; beamE<12; ++beamE){
		dHist_pi0eta1DBeam[beamE] = new TH1F(("pi0eta1DbeamE"+std::to_string(beamE)).c_str(),"M(pi0eta)",150,0,3);
	}
	dHist_pi0eta1D = new TH1F("pi0eta1D","M(pi0eta)",350,0,3.5);
	dHist_phi8GeVPlus = new TH1F("phi8GeVPlus","phi GJ",60,-180,180);
	dHist_cosTheta8GeVPlus = new TH1F("cosTheta8GeVPlus","cosTheta GJ",60,-1,1);
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
	maxevent=10;
}

Bool_t DSelector_thrown_7_17_14::Process(Long64_t locEntry)
{
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
	//if ( ievent<maxevent ){
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//
	//

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

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

	std::vector<TLorentzVector> allP4;
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
	vector<Int_t> vecParent;
	vector<Particle_t> vecParticle;
	Int_t numTotalEta = 0;
	Int_t numTotalPi0 = 0;
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{	
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
		Particle_t locPID = dThrownWrapper->Get_PID();
		Int_t locParentPID = dThrownWrapper->Get_ParentIndex();
		allP4.push_back(dThrownWrapper->Get_P4());
		//Some verbal confirmation checking
		if (is_in){
		// The loc_i'th particle would have a particle ID saying the particle type and a ParentPID that matches it with the loc_j'th particle to say its parent is the loc_j'th particle.
		cout << "Thrown " << loc_i << " with (PID,ParentPID) = (" << locPID << " ," << locParentPID << ")"  << endl;
		//cout << dThrownWrapper->Get_ParentIndex() << endl;
		}
		// seeding the decay chain..
		vecParticle.push_back(locPID);
		vecParent.push_back(locParentPID);
		dHist_PID->Fill(locPID);	
		if (locPID==17){ numTotalEta++;	}
		if (locPID==7){ numTotalPi0++;}	
	}
	

	///////////////// BOOLS to control which combos to cut ///////////////////
	bool etaTwoDaughters = 0;
	bool pi0TwoDaughters = 0;
	bool etaPhotonDaughters = 1;
	bool pi0PhotonDaughters = 1;
	bool oneInitialEta = 0;
	bool oneInitialPi0 = 0;
	bool oneEta = 0;
	bool onePi0 = 0;

	if (numTotalEta==1){ oneEta=1; }
	if (numTotalPi0==1){ onePi0=1; }

	//bool etaPi0NotTwo = 0; // the eta OR pi0 does not decay into two particles.
	//bool etaPi0NotPhoton = 0; // the eta OR pi0
	// the above two conditions are not enough to selec our event, since if there were no pi0 it would still pass, we introduce
	//bool notOneEtaOnePi0 = 0; 
	//////////////////////////////////////////////////////////////////////////
	
	vector<Int_t> idxInitial;
	set<Int_t> complete; // will be used to check if we use all the particles
	Int_t vecSize = static_cast<Int_t>(vecParent.size());
	// complete is a set complete set of intergers from 0 to vecSize
	for (Int_t i=0; i < vecSize; ++i){complete.insert(i);}
	set<Int_t> completing;
	Int_t matchId = 0;
        for (auto it = vecParent.begin(); it != vecParent.end(); ++it){
		// it is an iterator that iterates vecParent. *it would dereference to give the underlying value the pointer was pointing at. -1 = no parents
                if (*it == -1){
			// *it is = 1 when this particle has no parent
                        idxInitial.push_back(matchId);
			// say we have used this particle up.
			completing.insert(matchId);
                }
        	++matchId;
        }
	if (is_in) {
		for (int const& element : completing) { cout << element << " "; } 
		cout << "primaries" << endl;
	} 
	// since matchId would also correspond to the array index we can use it in looking for its daughters
	// searching for daughters of initial particles
	vector<vector<Int_t>> secondaries;
	vector<Int_t> secondary;
	vector<Int_t> vecSizeSec;
	// loop through the idxInital vector we just created with the parent particles AND the vector of parent particles
	for (auto it2 = idxInitial.begin(); it2 != idxInitial.end(); ++it2){
		matchId = 0;
		secondary.clear();
		if (is_in) {
			cout << *it2 << " has daughters : ";
		}
        	for (auto it = vecParent.begin(); it != vecParent.end(); ++it){
			if (*it == *it2){
				if (is_in){
					cout << matchId << ", ";
				}
				// basically this will be a vector of indices of particles that would have the same parent
				secondary.push_back(matchId);
				completing.insert(matchId);	
			}
			++matchId;
		}
		Int_t sizeSec = static_cast<Int_t>(secondary.size());
		if (is_in){ 
			cout << endl << "sizeSec: " << sizeSec << endl;
		}
		vecSizeSec.push_back(sizeSec);
		secondaries.push_back(secondary);
	}
	
	// here we check if we have used up all the particles
	if (completing != complete) {
		vecHasSecondaryDaughters.push_back(eventIdx);
	}
	
	//turns out there are decays with 3 fragementations. So we need to implement something for that.
	//secondaries_2 should have the same size and the number of daughters to the initial particles
	vector<vector<Int_t>> secondaries_2;
	vector<Int_t> secondary_2;
	for (auto it2 = secondaries.begin(); it2 != secondaries.end(); ++it2){
		//it2 is a pointer
		//here we iterate through all the secondaries where secondaries is a vector of vectors where the subvectors share the same parent
		for (auto it = it2->begin(); it != it2->end(); ++it){
			matchId = 0;
			secondary_2.clear();
			if (is_in){
				cout << *it << " has daughters : ";
			}
			// here we loop through the full vector of parent particles again to do our matching
			for (auto it3 = vecParent.begin(); it3 != vecParent.end(); ++it3){
			      //between vecParent and it2 which represents the daughters of the initial
				if (*it3 == *it){
			      		if (is_in){
						cout << matchId << ", ";		
					}
					secondary_2.push_back(matchId);
					completing.insert(matchId);
				}
				++matchId;
			}
			if (is_in){
				cout << endl;
			}
			secondaries_2.push_back(secondary_2);
		}
	}


	if (is_in) { 
		for (Int_t const& element : completing) { cout << element << " "; }; 
		cout << "completing" << endl;
		for (Int_t const& element : complete) { cout << element << " ";}
		cout << "complete" << endl;
	} 
	if (completing == complete) {
		vecEqual.push_back(eventIdx);
		++equal;
		}
	else {
		++notEqual;
		vecNotEqual.push_back(eventIdx);
	} 

	// check for etas in idxInitial
	Int_t eta = 17;
	Int_t pi0 = 7;
	Int_t proton =14;
	//itIter is indexing idxInitial, the initial particles, we have to keep an increment that mactches with itIter to use for indexing secondaries
	matchId = 0; //matchId will be the index of the initial particles that have a matching pi0 OR eta
	// checking that we have only one eta and one pi0
	Int_t numEta = 0;
	Int_t numPi0 = 0;
	Int_t startIdx;
	for (auto itIter = idxInitial.begin(); itIter != idxInitial.end(); ++itIter){
		numEta = 0; numPi0 = 0;
		// Only looks for an inital eta or pi0 that have no parents
		if (vecParticle[*itIter] == eta || vecParticle[*itIter] == pi0){
			vector<Int_t> etaPi0Daughters = secondaries[matchId]; // this will be a the daughters vector of the pi0 i.e. as a vector of two photons	
			Int_t sizeEtaPi0Daughters = vecSizeSec[matchId];
			if(vecParticle[*itIter] == proton){
				locProtonP4 = allP4[*itIter];
			}
			if (vecParticle[*itIter] == eta) {
				numEta += 1; 
				locEtaP4 = allP4[*itIter];
				if (sizeEtaPi0Daughters == 2){
					// checking eta has only two daughters....
					etaTwoDaughters = 1;	
				}
			}
			if (vecParticle[*itIter] == pi0) { 
				numPi0 += 1; 
				locPi0P4 = allP4[*itIter];
				if (sizeEtaPi0Daughters == 2){
					pi0TwoDaughters = 1;	
				}
			}
			if (numEta > 1 || numPi0 > 1){ throw std::invalid_argument( "numEta or numPi0 > 1!" ); }
			if (is_in){
				cout << "Found particle with PID " << vecParticle[*itIter] << " having daughters: " << endl;
			}
			startIdx = 0; 
			for (auto itSizes = vecSizeSec.begin(); itSizes != vecSizeSec.begin()+matchId;++itSizes){
				startIdx += *itSizes;
				//cout << "startIdx/itSizes: " << startIdx << "/" << *itSizes << endl;
			}
			//cout << "final startIdx: " << startIdx << endl;
			
			Int_t incrementStartIdx = 0;
			Int_t secondariesSize = static_cast<Int_t>(secondaries.size());
			for (auto itDaughter = etaPi0Daughters.begin(); itDaughter != etaPi0Daughters.end(); ++itDaughter){
				// first check if our decays are to photons. We initialize _PhotonDaughters to be 1. If we find any of the daughters to not be a gamma we set it to 0
				if (vecParticle[*itIter] == eta) {
					if (vecParticle[*itDaughter]!=1){ // photons have GEANT pid = 1
						etaPhotonDaughters = 0;
					}
				}
				if (vecParticle[*itIter] == pi0) {
					if (vecParticle[*itDaughter]!=1){ // photons have GEANT pid = 1
						etaPhotonDaughters = 0;
					}
				}
				// to grab the secondary daughters, where secondaries_2 has the same length as the total number of all daughters for all the initial
				// we first have to grab where the secondary daughters of the eta's begin which is denoted by startIdx and runs for sizeEtaDaughters
				if (is_in){
					cout <<  *itDaughter << " => ";
				}
				vector<Int_t> etaPi0Daughters2 = secondaries_2[startIdx+incrementStartIdx];
					if (etaPi0Daughters2.size() != 0){
						for (auto itDaughters2 = etaPi0Daughters2.begin(); itDaughters2 != etaPi0Daughters2.end(); ++itDaughters2){
							if (is_in){
								cout << "(" << *itDaughters2 << "," << vecParticle[*itDaughters2] << ") ";
							}
						}
					}
					else { if (is_in){cout << "(These are final states)";} }
					if (is_in){
						cout << " with (array_index, PID)" << endl;
					}
				++incrementStartIdx;
			}
			if (is_in){	
				cout << endl;
			}
		}
		++matchId;	
		if ( numEta == 1 ) { oneInitialEta = 1; } 
		if ( numPi0 == 1 ) { oneInitialPi0 = 1; } 
	}

	double lowE = 0;
	double uppE = 1;
	for (int beamE=0; beamE<12; ++beamE){
		//cout << "lowE, uppE: "<<lowE<<", "<<uppE<<endl;
		pBeamE[beamE] = lowE<locBeamP4.E() && locBeamP4.E()<uppE;
		lowE+=1;
		uppE+=1;
	}

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
        double cosTheta_eta_GJ = angles_eta.CosTheta();
        double phi_pi0_GJ = angles_pi0.Phi()*radToDeg;
        double phi_eta_GJ = angles_eta.Phi()*radToDeg;

	double locPi0EtaMass = locPi0EtaP4.M();
	double mandelstam_t = (locProtonP4-locTargetP4).M2();
	double mandelstam_t0 = TMath::Power((locBeamP4.M2()-locPi0EtaP4.M2()-locTargetP4.M2()+locProtonP4.M2())/(2*(locBeamP4+locTargetP4).M()),2)-(locBeamP4_cm-locPi0EtaP4_cm).M2();
	double mandelstam_tp = abs(mandelstam_t-mandelstam_t0);

	// the passCuts we use would depend on what we are looking for. Since our overall goal is to look at the pi0eta -> 4 photon final states. If we want to look at the full
	// distribution of eta and pi0 to see if there is an a0 or a2 resonance we should probably not restrict to 4 photon final states. 
	// We can probably use passCuts if we are trying to compare the q-value weighting of the mass distributions to the thrown distributions. Since we have not restricted on EBeam>8GeV
	// 	and we have used the fact that the eta and pi0 decay into two photons I think using this cut variable is appropriate.
	bool passCuts = beamEisFinite*etaTwoDaughters*pi0TwoDaughters*etaPhotonDaughters*pi0PhotonDaughters*oneInitialEta*oneInitialPi0;
	//bool passCutsOneEtaPi0 = beamEisFinite*oneInitialEta*oneInitialPi0;
	bool passCutsOneEtaPi0 = beamEisFinite*oneEta*onePi0;
	
	auto locNumThrown = Get_NumThrown();
	// There are some events with NumThrown=8; we have to check them out;
	// apparently there are states with NumThrown=4 also
	if (locNumThrown != 7){
		vecNumThrownNot7.push_back(eventIdx);
	}
	dHist_NumThrown->Fill(locNumThrown);
	dHist_cosThetaVsMass_tpAll->Fill(locPi0EtaMass,cosTheta_pi0_GJ);
	dHist_phiVsMass->Fill(locPi0EtaMass,phi_pi0_GJ);
	dHist_phi->Fill(phi_pi0_GJ);
	dHist_cosTheta->Fill(cosTheta_pi0_GJ);
	dHist_beamE->Fill(locBeamP4.E());
	for (int beamE=0; beamE<12; ++beamE){
		if (pBeamE[beamE]){
			dHist_pi0eta1DBeam[beamE]->Fill(locPi0EtaMass);
		}
	}
	pBeamE8to9GeV=locBeamP4.E()>8 && locBeamP4.E()<9;	
	if ( passCutsOneEtaPi0 && mandelstam_tp < 1) {
		//Fill_OutputTree("selected_tLT1"); //your user-defined key
	}
	if (passCutsOneEtaPi0*pBeamE8to9GeV){
		//Fill_OutputTree must be run with proof!
		mandelstam_tpAll->Fill(mandelstam_tp);
		dHist_pi0eta1D->Fill(locPi0EtaMass);
		//cout << locPi0P4.M() << endl;
		//cout << locEtaP4.M() << endl;
		dHist_cosThetaVsMass_tpAll->Fill(locPi0EtaMass,cosTheta_pi0_GJ);
		dHist_phi8GeVPlus->Fill(phi_pi0_GJ);
		dHist_cosTheta8GeVPlus->Fill(cosTheta_pi0_GJ);
		if(mandelstam_tp < 1){
			Fill_OutputTree("selected_tLT1"); //your user-defined key
			mandelstam_tpLT1->Fill(mandelstam_tp);
			dHist_cosThetaVsMass_tpLT1->Fill(locPi0EtaMass,cosTheta_pi0_GJ);
		}
		if(mandelstam_tp < 0.6){
			Fill_OutputTree("selected_tLT06"); //your user-defined key
			mandelstam_tpLT06->Fill(mandelstam_tp);
			dHist_cosThetaVsMass_tpLT06->Fill(locPi0EtaMass,cosTheta_pi0_GJ);
		}
		if((mandelstam_tp >= 0.5) && (mandelstam_tp < 1)) {
			Fill_OutputTree("selected_tGT05LT1"); //your user-defined key
			mandelstam_tpGT05LT1->Fill(mandelstam_tp);
			dHist_cosThetaVsMass_tpGT05LT1->Fill(locPi0EtaMass,cosTheta_pi0_GJ);
		}
		pass += 1;
		Fill_OutputTree();
	}
	else {	
		vecNotPassCuts.push_back(eventIdx);
		notPass += 1;
	}
	++eventIdx;

	//Fill_OutputTree();


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
	//}
	ievent++;
	return kTRUE;
}

void DSelector_thrown_7_17_14::Finalize(void)
{

	cout << "sets are equal: " << equal << " not equal: " << notEqual << endl;
	cout << "Some of the eventIdx which the sets are not equal: " << endl;
	Int_t N = 0;
	for (auto itVNE = vecNotEqual.begin(); itVNE != vecNotEqual.end(); ++itVNE){
		if (N < 10){
			cout << *itVNE << " ";
		} 
		++N;
	}
	N = 0;
	cout << "Pass cuts: " << pass << " Not pass cuts: " << notPass << endl;
	cout << "Some of the eventIdx which they are not pass cuts: " << endl;
	for (auto itVNPC = vecNotPassCuts.begin(); itVNPC != vecNotPassCuts.end(); ++itVNPC){
		if (N < 10){
			cout << *itVNPC << " ";
		} 
		++N;
	}
	N = 0;
	cout << endl << "Some of the eventIdx which has secondary Daughters: " << endl;
	for (auto itVHSD = vecHasSecondaryDaughters.begin(); itVHSD != vecHasSecondaryDaughters.end(); ++itVHSD){
		if (N < 10){
			cout << *itVHSD << " ";
		} 
		++N;
	}
	N = 0;
	cout << endl << "Some of the eventIdx which has NumThrown != 7: " << endl;
	for (auto itVNTN7 = vecNumThrownNot7.begin(); itVNTN7 != vecNumThrownNot7.end(); ++itVNTN7){
		if (N < 10){
			cout << *itVNTN7 << " ";
		} 
		++N;
	}

	cout << endl << "beamEinf: " << beamEinf << " beamEfin: " << beamEfin << endl;

	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
