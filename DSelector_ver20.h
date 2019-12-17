#ifndef DSelector_ver20_h
#define DSelector_ver20_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector3.h"
#include "THStack.h"
#include "TCanvas.h"

bool is_pi0eta=true;

void findDaughters( std::vector<int> parentArray, std::vector<int> &daughters, int selfLoc ) {
        for ( auto parent=0; parent<(int)parentArray.size(); ++parent ){
                if ( parentArray[parent] == selfLoc ){
                        daughters.push_back(parent);
                }
        }
}

void findParents( std::vector<int> parentArray, std::vector<int> &parents) {
	for (auto parent=0; parent<(int)parentArray.size(); ++parent) {
		if ( parentArray[parent] == -1 ){
			parents.push_back(parent);
		}
	}
}

void withinBox(bool inBox[], bool inBox_noOtherCuts[],bool additionalCut, double x, double y, double xmin, double xmax, double ymin, double ymax, double xskip, double yskip){
	// regions are:
	//  0 1 2
	//  3 4 5 
	//  6 7 8
	//  where xmin, xmax, ymin,ymax all belong to region 5
	//  the 10th element is the interesection of the negation of all regions 
	//  11th element is 1345 and 12th element is 0268
	double xlength = xmax-xmin;
	double ylength = ymax-ymin;
	inBox[4] = x<xmax && x>xmin && y<ymax && y>ymin;
	inBox[3] = x<(xmin-xskip) && x>(xmin-xskip-xlength/2) && y<ymax && y>ymin;
	inBox[5] = x<(xmax+xskip+xlength/2) && x>(xmax+xskip) && y<ymax && y>ymin;
	inBox[1] = x<xmax && x>xmin && y<(ymax+yskip+ylength) && y>(ymax+yskip);
	inBox[7] = x<xmax && x>xmin && y>(ymin-yskip-ylength) && y<(ymin-yskip);
	inBox[0] = x<(xmin-xskip) && x>(xmin-xskip-xlength/2) && y<(ymax+yskip+ylength) && y>(ymax+yskip);
	inBox[2] = x<(xmax+xskip+xlength/2) && x>(xmax+xskip) && y<(ymax+yskip+ylength) && y>(ymax+yskip);
	inBox[6] = x<(xmin-xskip) && x>(xmin-xskip-xlength/2) &&  y>(ymin-yskip-ylength) && y<(ymin-yskip);
	inBox[8] = x<(xmax+xskip+xlength/2) && x>(xmax+xskip) && y>(ymin-yskip-ylength) && y<(ymin-yskip);
	inBox[9] = !inBox[0] * !inBox[1] * !inBox[2] * !inBox[3] * !inBox[4] * !inBox[5] * !inBox[6] * !inBox[7] * !inBox[8];	
	inBox[10] = inBox[0] || inBox[2] || inBox[6] || inBox[8];
	inBox[11] = inBox[1] || inBox[7];
       	inBox[12] = inBox[3] || inBox[5];

	for (int i=0; i<13; ++i){
		inBox_noOtherCuts[i] = inBox[i]; // save the inBox bools before modifiying them by the cut. We will use these for defining the weightBS
		inBox[i]*=additionalCut;
	}
}

struct histDef_1D{
    TH1F* hist;
    string name;
    bool* cut;
    std::vector< double* > values;
    double* weights; 
    void clear(){
        values.clear();
    }
};

struct histDef_2D{
    TH2F* hist;
    string name;
    bool* cut;
    std::vector< double* > valuesX;
    std::vector< double* > valuesY;
    double* weights; 
    void clear(){
        valuesX.clear();
        valuesY.clear();
    }
};

class trackingGroup{
    private:
	//set< pair<Int_t,Int_t> > usedIds;
        set< pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > > usedPairMapIds;
        set< map<Particle_t, set<Int_t> > >  usedMapIds;

    public:
        TCanvas *anyCanvas;
        trackingGroup(){ cout << "Calling constructor" << endl; }; 
        ~trackingGroup(){ cout << "Calling destructor" << endl; };

        void setNameMakeCanvas(string groupName){
	    anyCanvas = new TCanvas((groupName).c_str(),"",1440, 900);
        }

        std::vector<histDef_1D> allHists_1D;
        std::vector<histDef_2D> allHists_2D;

        std::vector< set< pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > > > allUsedPairMapIds_2D;
        std::vector< set< pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > > > allUsedPairMapIds_1D;
        std::vector< set< map<Particle_t, set<Int_t> > > > allUsedMapIds_1D;
        std::vector< set< map<Particle_t, set<Int_t> > > > allUsedMapIds_2D;

        // every time we add another histogram to track we will push back another tracking set to track it.
        void insert(histDef_1D hist){
            allHists_1D.push_back(hist);
            allUsedMapIds_1D.push_back(usedMapIds);
            allUsedPairMapIds_1D.push_back(usedPairMapIds); // we only need this for 2D histograms
        }
        void insert_2D(histDef_2D hist){
            // ***** NOTE THAT IT IS BEST IF WE DO NOT MIX USING BOTH OF TYPES OF TRACKING TYPES ***** NOT WORKED ON YET
            allHists_2D.push_back(hist);
            allUsedMapIds_2D.push_back(usedMapIds);
            allUsedPairMapIds_2D.push_back(usedPairMapIds); // we only need this for 2D histograms
        }

        // for filling distance between pairs of particles like dij3FCAL
        //void fillHistograms_vectorMap( std::vector< map<Particle_t, set<Int_t> > > beingUsedPairIds){
        //    for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
        //        for (UInt_t iVal=0; iVal<allHists_1D[iHist].values.size(); ++iVal){
	//            if ( *(allHists_1D[iHist].cut) ) {
        //               allHists_1D[iHist].hist->Fill( *(allHists_1D[iHist].values[iVal]), *(allHists_1D[iHist].weights) ); 
        //            }
        //        }
        //    }
        void fillHistograms_vectorMap( std::vector< map<Particle_t, set<Int_t> > > beingUsedPairIds){
            // *********************** THIS COMMENTED CODE IS TO CHECK THE OUTPUT OF THE DIJ3 FCAL FILLING TO MAKE SURE WE ARE DOING THINGS CORRECLTY ***********************************
            //  cout << "size of {tracking,value}={"<<beingUsedPairIds.size()<<","<<allHists_1D[0].values.size()<<"}"<<endl;
            // for (UInt_t iValue=0; iValue<beingUsedPairIds.size(); ++iValue){
            //    for ( auto elem : beingUsedPairIds[iValue] ){
            //        cout << elem.first << "| ";
            //        for (auto it=elem.second.begin(); it != elem.second.end(); ++it){
            //            cout << *it << " ";
            //        }
            //    }
            //bool uniqueBool= allUsedMapIds_1D[0].find(beingUsedPairIds[iValue])==allUsedMapIds_1D[0].end();
            //bool cutBool =  *(allHists_1D[0].cut);
            //cout << " -- " << *(allHists_1D[0].values[iValue]) << " Unique? " << uniqueBool << " passCut? " << cutBool << endl; 
            //}
                
             for (UInt_t iValue=0; iValue<beingUsedPairIds.size(); ++iValue){
                 for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
	              if (allUsedMapIds_1D[iHist].find(beingUsedPairIds[iValue])==allUsedMapIds_1D[iHist].end() && *(allHists_1D[iHist].cut)  ){
	                 allUsedMapIds_1D[iHist].insert(beingUsedPairIds[iValue]); //we get a iterator which references the element of the set so we need to dereference.              
                         allHists_1D[iHist].hist->Fill( *(allHists_1D[iHist].values[iValue]), *(allHists_1D[iHist].weights) );
                      }
                 }
             }

             for (UInt_t iValue=0; iValue<beingUsedPairIds.size(); ++iValue){
                 for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
	              if (allUsedMapIds_2D[iHist].find(beingUsedPairIds[iValue])==allUsedMapIds_2D[iHist].end() && *(allHists_2D[iHist].cut)  ){
	                 allUsedMapIds_2D[iHist].insert(beingUsedPairIds[iValue]); //we get a iterator which references the element of the set so we need to dereference.              
                         allHists_2D[iHist].hist->Fill( *(allHists_2D[iHist].valuesX[iValue]),*(allHists_2D[iHist].valuesY[iValue]), *(allHists_2D[iHist].weights) );
                      }
                 }
             }
            clear_values();
        }

        // Similar to above but just for 1D distributions where we don't need to track a pair of maps for correlations
        void fillHistograms_Map( map<Particle_t, set<Int_t> >  beingUsedMap ){
             for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
	        if (allUsedMapIds_1D[iHist].find(beingUsedMap)==allUsedMapIds_1D[iHist].end() && *(allHists_1D[iHist].cut)  ){
	            allUsedMapIds_1D[iHist].insert(beingUsedMap); //we get a iterator which references the element of the set so we need to dereference.              
                    allHists_1D[iHist].hist->Fill( *(allHists_1D[iHist].values[0]), *(allHists_1D[iHist].weights) );
                }
             }
             for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
	         if (allUsedMapIds_2D[iHist].find(beingUsedMap)==allUsedMapIds_2D[iHist].end() && *(allHists_2D[iHist].cut)  ){
	            allUsedMapIds_2D[iHist].insert(beingUsedMap); //we get a iterator which references the element of the set so we need to dereference.              
                    allHists_2D[iHist].hist->Fill( *(allHists_2D[iHist].valuesX[0]),*(allHists_2D[iHist].valuesY[0]), *(allHists_2D[iHist].weights) );
                }
             }
            clear_values();
        }

        // This will be used for more complicated sets. Most generic is a pair of maps. i.e. for 2D plot of M(pi0eta) and cosThetaGJ(eta)
        void fillHistograms_pairMap( pair< map<Particle_t, set<Int_t> >, map<Particle_t, set<Int_t> > > beingUsedMapPair ){
             for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
	        if (allUsedPairMapIds_1D[iHist].find(beingUsedMapPair)==allUsedPairMapIds_1D[iHist].end() && *(allHists_1D[iHist].cut)){
	            allUsedPairMapIds_1D[iHist].insert(beingUsedMapPair); //we get a iterator which references the element of the set so we need to dereference.              
                    allHists_1D[iHist].hist->Fill( *(allHists_1D[iHist].values[0]), *(allHists_1D[iHist].weights) );
                }
             }
            for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
	        if (allUsedPairMapIds_2D[iHist].find(beingUsedMapPair)==allUsedPairMapIds_2D[iHist].end() && *(allHists_2D[iHist].cut)){
	            allUsedPairMapIds_2D[iHist].insert(beingUsedMapPair); //we get a iterator which references the element of the set so we need to dereference.              
                    allHists_2D[iHist].hist->Fill( *(allHists_2D[iHist].valuesX[0]),*(allHists_2D[iHist].valuesY[0]), *(allHists_2D[iHist].weights) );
                }
            }
            clear_values();
        }

        void drawHistograms(){
            for (UInt_t iHist=0; iHist < allHists_1D.size(); ++iHist){
                anyCanvas->Clear();
                allHists_1D[iHist].hist->Draw("HIST");
                anyCanvas->SaveAs(("graphs/"+allHists_1D[iHist].name+".png").c_str());
            }
            for (UInt_t iHist=0; iHist < allHists_2D.size(); ++iHist){
                anyCanvas->Clear();
                allHists_2D[iHist].hist->Draw();
                anyCanvas->SaveAs(("graphs/"+allHists_2D[iHist].name+".png").c_str());
            }
        }

        // use when we are starting a new event
        void clear_tracking(){
            for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
                allUsedMapIds_1D[iHist].clear();
                allUsedPairMapIds_1D[iHist].clear();
            }
            for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
                allUsedMapIds_2D[iHist].clear();
                allUsedPairMapIds_2D[iHist].clear();
            }
        }
        
        // have to use this everytime we fillHistogram or changing combos
        void clear_values(){
            for (UInt_t iHist=0; iHist<allHists_1D.size(); ++iHist){
                allHists_1D[iHist].clear();
            }
            for (UInt_t iHist=0; iHist<allHists_2D.size(); ++iHist){
                allHists_2D[iHist].clear();
            }
        }
};

    
        

class DSelector_ver20 : public DSelector
{
	public:
		DSelector_ver20(TTree* locTree = NULL) : DSelector(locTree){
                }
		virtual ~DSelector_ver20(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:
                trackingGroup group_PhNB;
                trackingGroup group_PB;
                trackingGroup group_12B_1234B;
                trackingGroup group_34B_1234B;
                trackingGroup group_12B;
                trackingGroup group_34B;
                trackingGroup group_12B_34B;
                trackingGroup group_1234B;
                trackingGroup group_1234BP;
                trackingGroup group_pairFCAL;
                trackingGroup group_pairBCAL;
                trackingGroup group_12PB;       
                trackingGroup group_34PB;       
                trackingGroup group_B;
                trackingGroup group_1234B_12PB;
                trackingGroup group_1234B_34PB;

		void Get_ComboWrappers(void);
		void Finalize(void);

                // BEAM POLARIZATION INFORMATION
                UInt_t dPreviousRunNumber;
                bool dIsPolarizedFlag; //else is AMO
                bool dIsPARAFlag; //else is PERP or AMO
                bool hasPolarizationAngle; // true if is polarized but false if through deg5ous radiator or no data. Under what cicumstances is the second one true.
                int locPolarizationAngle; // actual polarization angle
                set<UInt_t> usedRuns; //think we are over counting our filling of beam angles since the below condition only checks against previous run, doesn't work for recurrences.


		// These variables are for the teta vs Meta plots
		int num_tBins=14;
		double tMin=0;
		double tMax=2.8;
		static const int num_massBins=12;
		const int numHists = num_tBins*num_massBins;
		double mMin=1.7;
		double mMax=2.9;
		double tStep=(tMax-tMin)/num_tBins;
		double mStep=(mMax-mMin)/num_massBins;
		double teta_recCounts=0;
		double tpi0_recCounts=0;
		int idx_t_eta;
		int idx_t_pi0;
		int idx_m;
		bool passMassBin_tetaIntegrated[num_massBins];
		bool passMassBin_tpi0Integrated[num_massBins];



                Int_t eventIdx=0;
                bool isNotRepeated_eta=true;
                bool isNotRepeated_pi0=true;
                bool isNotRepeated_pi0eta=true;
                bool isNotRepeated_pi0_pi0eta=true;
                bool isNotRepeated_eta_pi0eta=true;
                Int_t uniqueSpectroscopicPi0EtaID=0;

		Int_t uniqueComboID=0;

                //string degAngle = "deg0";

                // ANALYZE CUT ACTIONS
                // // Automatically makes mass histograms where one cut is missing
                DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		TH1F* dHistThrownTopologies;
		map<TString, TH1I*> dHistInvariantMass_ThrownTopology;

                //CREATE REACTION-SPECIFIC PARTICLE ARRAYS

                //Step 0
                DParticleComboStep* dStep0Wrapper;
                DBeamParticle* dComboBeamWrapper;
                DChargedTrackHypothesis* dProtonWrapper;

                //Step 1
                DParticleComboStep* dStep1Wrapper;
                DKinematicData* dDecayingPi0Wrapper;
                DNeutralParticleHypothesis* dPhoton1Wrapper;
                DNeutralParticleHypothesis* dPhoton2Wrapper;

                //Step 2
                DParticleComboStep* dStep2Wrapper;
                DKinematicData* dDecayingEtaWrapper;
                DNeutralParticleHypothesis* dPhoton3Wrapper;
                DNeutralParticleHypothesis* dPhoton4Wrapper;

                // DEFINE YOUR HISTOGRAMS HERE
                // EXAMPLES:
		TVector3 targetCenter;
	
		// This is for calculating the masses with different starting points
		double locEtaMass_charged=1;
		double locPi0Mass_charged=1;
		double locEtaMass_target=1;
		double locPi0Mass_target=1;


                TH1F* dHist_BeamAngle;
		TH1F* dHist_Cuts;
		TH2F* dHist_checkEllipseBS[3];
		const char *cutNames[15] = {"pShowerQuality","pBeamE8GeVPlus","pUnusedEnergy","pChiSq" ,"pDeltaTRF","pdij3pass","pPhotonE","pPhotonTheta","pMagP3Proton","pzCutmin","pRProton","pMissingMassSquared","pdEdxCDCProton","pinsideEllipse", "allGeneralCutsPassed"};


// **************************************** START INITIALZING VARIABLES TO USE WITH HISTO BUILDING ********************************************//
// **************************************** START INITIALZING VARIABLES TO USE WITH HISTO BUILDING ********************************************//

// Defining the maps to hold the histograms parameters 
		std::map<int, std::vector<std::string>> histList;
		std::map<int, std::vector<Double_t>> histVals;
		std::map<int, std::vector<std::vector<Double_t> > > histVecVals;
		std::map<int, std::vector<bool>> histVecCuts;
		// we cannot make a std::vector<bool> since it gives an error about "narrowing" when we do something like a product of booleans in vector initialization. i.e. {cutsToApply*pReg1}
		std::map<int, bool> histCuts;

	        // cutString will be used by alot of histograms when defining names 
	        std::string cutString;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////// ********* DEFINING VALUES TO FILL HISTS WITH ****************/////////////////////////////////////////////

		// will be used to bin the kinematic phi and t distributions
		double lowMass;
		double upMass;
		// we need to use static const since we will use this variable to make the array. This region of the code is probably the "file scope" and declares things at compile time instead of run time. 
		// The array is in this region so it must take in a static object where the preprocessor will just replace the variables. 
		static const int numBinsMass = 13;
		double binScale;

		double lowMass_t = 0.8;
		double upMass_t = 1.8;
		static const int numBinsMass_t = 10;
		double binScale_t = (upMass_t-lowMass_t)/numBinsMass_t;

		double radToDeg = 180./TMath::Pi();
		double locWherePhoton=1;

		// Accidental subtraction variables. applyAccSub will either be = to weight or noAccSub=1.
		double weightAS=1;
		double noWeight=1;
                double weightAS_BS=1;
                double weightAS_B=1;
                double appliedWeight=1;

		//************* Other variables
		double locMissingMassSquared=1;
		double locBeamE=1;
		double locCLKinFit=1;
		double locUnusedEnergy=1;
		double locNumExtraNeutralShowers=1;
		double locChiSqKinFit=1;
		double locDOFKinFit=1;

		//************ Particle related variables  
		double locEtaE_Kin=1;
		double locPi0E_Kin=1;
		double locEtaMass_Kin=1;
		double locPi0Mass_Kin=1;
		double locEtaMass_Kin_mismatch=1;
		double locPi0Mass_Kin_mismatch=1;
		double locEtaProton_Kin=1;
		double locPi0Proton_Kin=1;
		double locPi0Eta_Kin=1;
		double locPi0Eta=1;

		double locEtaMass=1;
		double locPi0Mass=1;
		//************ Shower shape variables
		std::vector<double> E1E9_FCAL={1,1,1,1};
		std::vector<double> E9E25_FCAL={1,1,1,1};
		std::vector<double> SumU_FCAL={1,1,1,1};
		std::vector<double> SumV_FCAL={1,1,1,1};
		std::vector<double> Energy_BCALPreshower={1,1,1,1};
		std::vector<double> Energy_BCAL={1,1,1,1};
		std::vector<double> SigLong_BCAL={1,1,1,1};
		std::vector<double> SigTrans_BCAL={1,1,1,1};
		std::vector<double> SigTheta_BCAL={1,1,1,1};
		std::vector<double> DeltaPhi_BCAL={1,1,1,1};
		std::vector<double> DeltaZ_BCAL={1,1,1,1};
		std::vector<double> showerQuality_FCAL={1,1,1,1};
		std::vector<double> DOCA_FCAL={1,1,1,1};
		double locSigTheta_BCAL_proton;
		double locSigTrans_BCAL_proton;
		double locSigLong_BCAL_proton;
		double locE1E9_FCAL_proton;
		double locE9E25_FCAL_proton;
		double locSumU_FCAL_proton;
		double locSumV_FCAL_proton;
		double locEnergy_BCALPreshower_proton;
		double locEnergy_BCAL_proton;

		double pi0DetectedIn;
		double etaDetectedIn;

		//************ Charged Track
		double locPtProton=1;
		double locPzProton=1;
		double locPolarAngleProton=1;
		double locXProton=1;
		double locYProton=1;
		double locRProton=1;
		double locdzProton=1;
		double locdEdxCDCProton=1;
		double locdEdxFDCProton=1;
		double locMagP3Proton=1;

		//************ Neutral Track
		std::vector<double> photonThetas={1,1,1,1};
		std::vector<double> photonEnergies={1,1,1,1};
		std::vector<double> photonPhis={1,1,1,1};
		std::vector<double> photonXs_Kin={1,1,1,1};
		std::vector<double> photonYs_Kin={1,1,1,1};
		std::vector<double> photonZs_Kin={1,1,1,1};
		std::vector<double> photonTs_Kin={1,1,1,1};
		std::vector<double> photonXs_Shower={1,1,1,1};
		std::vector<double> photonYs_Shower={1,1,1,1};
		std::vector<double> photonZs_Shower={1,1,1,1};
		std::vector<double> photonTs_Shower={1,1,1,1};
		std::vector<double> photonDeltaTs={1,1,1,1};
		std::vector<double> photonDetectedSyss={1,1,1,1};


		// distance, angle, Z distance, phi between photon pairs. They must be vectors since we are not sure about the size of them yet so this will allow us to push_back
		//std::vector<double> dij3Vec;
		//std::vector<double> angle_ijVec; 
		//std::vector<double> deltaZ_ijVec; 
		//std::vector<double> deltaPhi_ijVec;
		// single version of the above
		double locPhotonDijFCAL=1;
		double locPhotonDijBCAL=1;
		double locPhotonAij=1;
		double locPhotonZij=1;
		double locPhotonPij=1;
		
		//*********** Kinematic variables
		double locDecayPlaneTheta=1;
		double locPhi=1;
		// Calculating kinematic variables like t and cosTheta
		double mandelstam_teta=1;
		double mandelstam_teta_Kin=1;
		double mandelstam_tpi0=1;
		double mandelstam_tpi0_Kin=1;
		double mandelstam_tp=1;
		double mandelstam_tp_pe=1;
		double mandelstam_t=1;
		double mandelstam_abst=1;
		double mandelstam_t_pe=1;
		double mandelstam_t0=1;
		// Calculate cosTheta in maybe the gottfried-jackson frame.
		double theta_pi0_lab=1;
		double theta_eta_lab=1;
		double cosTheta_decayPlane_hel=1;
		double phi_decayPlane_hel=1;
		double cosTheta_decayPlane_GJ=1;
		double cosTheta_largestEinEta_GJ=1;
		double cosTheta_largestEinPi0_CM=1;
		double phi_decayPlane_GJ=1;
		double cosTheta_pi0_hel=1;
		double theta_pi0_hel=1;
		double theta_eta_hel=1;
		double phi_pi0_hel=1;
		double theta_pi0_GJ=1;
		double phi_pi0_GJ=1;
		double cosTheta_eta_hel=1;
		double phi_eta_hel=1;
		double theta_eta_GJ=1;
		double phi_eta_GJ=1;
		double phi_eta_GJ_meas=1;
		double cosTheta_pi0eta_hel=1;
		double phi_pi0eta_hel=1;
		double cosTheta_pi0_GJ=1;
		double cosTheta_eta_GJ=1;
		double cosTheta_eta_GJ_meas=1;
		double phi_pi0eta_GJ=1;

		double vanHove_x;
		double vanHove_y;
                double q;
                double pi0_cmZ;
                double eta_cmZ;
                double recoil_cmZ;
                double omega;


		double locYDotZ_GJ=1;
		double angleBetweenPi0Eta=1;
		int delta[4] = {5,10};
		bool withinCone[4]={true,true,true,true};
		bool pi0_inCone[4]={true,true,true,true};
		bool eta_inCone[4]={true,true,true,true};
		bool largeAngle[4]={true,true,true,true};
		std::vector<double> countCone = {0,1,2,3,4};
		// calculating the cosTheta of pi0eta system, pi0, eta  in the CM framei
		double cosTheta_pi0eta_CM=1;
		double cosTheta_pi0eta_CM_meas=1;
		double cosTheta_pi0_CM=1;
		double cosTheta_eta_CM=1;
		double phi_pi0eta_CM=1;
		double phi_pi0_CM=1;
		double phi_eta_CM=1;
		double theta_pi0_CM=1;
		double theta_eta_CM=1;
		
		//*********** Timing variables, first section is about the proton
		// shifted relative to the beam
		double RFtime=1;
		double RFtimeProton=1;
		// Might not make sense to shift relative to the BeamX4.Z since we cant really track the photon in the detector. Use Protons maybe
		//double locDeltaTRF = locBeamX4.T() - (dComboWrapper->Get_RFTime() + (locBeamX4.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
		double locDeltaTRF=1;
		
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////// ********* DEFINING ADDITIONAL CUTS TO APPLY TO GRAPHS ****************/////////////////////////////////////////////
		
		
		// ******************************** DEFINING CUT THRESHOLDS ***************************************
		// Beam asymmetry stuff
		double Emin = 8.2;
		double Emax = 8.8;
		
		/////////////// Charged Track Cuts/////////
		double Rcut = 2;
		double zCutmin = 42; double zCutmax = 82;
		double dEdxCut; // not static so we need to recompute it. TMath::Power(10,-6)*(0.9+TMath::Exp(3.0-3.3*locMagP3Proton/.93827)); // The bigger then number multiplying MagP3 the sharper the cut. 
		double P3Cut = .3;
		// zooming in on 4 specific regions in the following PzPt graphs
		double Reg1Xmin = 0.6; double Reg1Xmax = 1.1; double Reg1Ymin = 0.05; double Reg1Ymax = 0.3;
		double Reg2Xmin = 0.3; double Reg2Xmax = 0.8; double Reg2Ymin = 0.35; double Reg2Ymax = 0.65;
		double Reg3Xmin = 0.1; double Reg3Xmax = 0.5; double Reg3Ymin = 0.; double Reg3Ymax = 0.35;
		double Reg4Xmin = 0.6; double Reg4Xmax = 1.5; double Reg4Ymin = 0.2; double Reg4Ymax = 0.65;
		////////////// Photons //////////////////
		double ECut = 0.1;
		double thetaCutMin = 2.5; double thetaCutMax1 = 10.3; double thetaCutMax2 = 11.5;
		double dijCut = 12.5;
		///////////// General ///////////////////
		double unusedEnergyCut = 0.010;
		double MMsqCut = 0.05;
		//double CLCut1 = 0.1;
		//double CLCut = 0.01; // CLCut is the cut we want to uses for all other cuts, these others ones are for the RFTime graph
		double ChiSqCut = 13.277;
		double chiSq100 = 100;
		//double CLCut3 = 0.001;
		//double CLCut4 = 0.0001;
		//double CLCut5 = 0.00001;
		//double CLCut6 = 0.000001;
		double RFCut = 0.5*4; // 4ns is the beam period.
		double beamECut = 6;
		double etaProtonBaryonCut;
		double pi0ProtonBaryonCut;
		// the equation of the TEllipse follows (0.135,0.47,0.0175,0.15) which should bound locPi0Mass_Kin and locEtaMass_Kin
		//double ellipseX = 0.135; double ellipseY = 0.47; double ellipseXr = 0.0175; double ellipseYr = 0.15;
		// For 25000 events root file with no vertex in kinFit we use (0.1325,0.54,0.014,0.055)
		//
		// These are the 3 sigma regions from the graph of the pi0Mass with mUE applied. 
		//  mUE
		//double ellipseX = 0.134581; double ellipseY = 0.539935; double ellipseXr = 0.024066; double ellipseYr = 0.066594;
		// mUEChiSq
		double ellipseX; double ellipseY; double ellipseXr; double ellipseYr; 
    		double ellipseXBS1; double ellipseYBS1; double ellipseXrBS1; double ellipseYrBS1;
    		double ellipseXBS2; double ellipseYBS2; double ellipseXrBS2; double ellipseYrBS2;
		double skipX, skipY;
		double ellipseXr_loose, ellipseYr_loose;
		double weightBS=1;
		double weightB=1;
		double areaRatio=1;
		
		// Beam cuts
		bool pBeamAsymE=true;
		bool pBeamE30to46=true;
		bool pBeamE46to62=true;
		bool pBeamE62to78=true;
		bool pBeamE78to94=true;
		bool pBeamE94to11=true;
		bool pBeamE8GeVPlus=true;
		bool pBeamE8to9=true;
		
		// pi0Eta specifc cuts
		bool pEtaProtonBaryonCut=true;
		bool ppi0ProtonBaryonCut=true;
		bool pBeamE=true;
		bool p_phiMassBinned[numBinsMass]; 
		double iLowMass;
		double iUpMass;
		bool p_tMassBinned[numBinsMass_t]; 
		double iLowMass_t;
		double iUpMass_t;

		bool p_massTBinned[10];

		static const int numRegions_UE=10;
		static const int numRegions_ChiSq=10;
		bool p_pi0MassEtaMassUEregion[numRegions_UE]; 
		bool p_pi0MassEtaMassChiSqregion[numRegions_ChiSq]; 
		double iUpUE;
		double iUpChiSq;
		
		// FOR THE pi0 MASS BINNED IN E(pi0) selecting out the f2(1270) region
		static const int numRegions_E=17;
		bool p_pi0MassPi0Eregion_1[numRegions_E]; 
		bool p_pi0MassPi0Eregion_2[numRegions_E]; 
		double iLowE[numRegions_E]  = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.5};
		double iUpE[numRegions_E] = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.5,9.5};
		bool pSelectf2=true;
		// General Cuts
		bool pVanHove=true;
		bool mBeamE=true;
		bool pUnusedEnergy=true;
		bool pChiSq=true;
		//bool pCLKinFit1=true;
		//bool pCLKinFit=true ;
		//bool pCLKinFit3=true;
		//bool pCLKinFit4=true;
		//bool pCLKinFit5=true;
		//bool pCLKinFit6=true;
		bool pDeltaTRF=true ;
		bool pMissingMassSquared=true;
		
		bool pUEregion0=true;
		bool pUEregion1=true;
		bool pUEregion2=true ;
		bool pUEregion3=true;
		bool pUEregion4=true;
		bool pUEregion5=true;
		bool pUEregion6=true;
		bool pUEregion7=true;
		bool pUEregion8=true;
		bool pUEregion9=true;
		// NeuUEregiontral Cuts
		bool outsideEllipse_loose=true;
		bool pinsideEllipse_loose=true;
		bool outsideEllipse=true;
		bool pinsideEllipse=true;
		bool outsideEllipseBS1=true;
		bool pinsideEllipseBS1=true;
		bool outsideEllipseBS2=true;
		bool pinsideEllipseBS2=true;
		
		bool inBox[13]={0,0,0,0,0,0,0,0,0,0,0,0,0};
		bool inBox_noOtherCuts[13]={0,0,0,0,0,0,0,0,0,0,0,0,0};

		bool pYellowBKG=true;
		bool pdij3pass=true;
		bool pPhoton1E=true;
		bool pPhoton2E=true;
		bool pPhoton3E=true;
		bool pPhoton4E=true;
		bool pPhotonE=true ;
		bool pPhoton1Theta=true;
		bool pPhoton2Theta=true;
		bool pPhoton3Theta=true;
		bool pPhoton4Theta=true;
		bool pPhotonTheta=true;
		
		// Charged Cuts
		bool pMagP3Proton=true;
		bool pzCutmin=true;
		bool pRProton=true;
		bool pdEdxCDCProton=true;
		bool pReg1=true;
		bool pReg2=true;
		bool pReg3=true;
		bool pReg4=true;
		
		//showerQuality
		bool pShowerQuality0 = true;
		bool pShowerQuality1 = true;
		bool pShowerQuality2 = true;
		bool pShowerQuality3 = true;
		bool pShowerQuality = true;


		// locWherePhoton will be set equal to photonDetectedSyss[N] 
		bool pPhotonInBCALorFCAL[4];
		bool pPhotonInFCAL[4];
		bool pPhotonInBCAL[4];
		bool pPi0InFCAL=true;
		bool pPi0InBCAL=true;
		bool pPi0InSplit=true;
		bool pEtaInFCAL=true; 
		bool pEtaInBCAL=true; 
		bool pEtaInSplit=true;
		bool pPi0InFCAL_mismatch=true;
		bool pPi0InBCAL_mismatch=true;
		bool pPi0InSplit_mismatch=true;
		bool pEtaInFCAL_mismatch=true; 
		bool pEtaInBCAL_mismatch=true; 
		bool pEtaInSplit_mismatch=true;

		bool detectorCut=true;

		bool ptpLT1=true;
		// Times it passes a cut
		int count_ShowerQuality=0;
		int count_BeamE8GeVPlus=0;
		int count_UnusedEnergy=0;
		int count_ChiSq=0;
		int count_DeltaTRF=0;
		int count_dij3pass=0;
		int count_PhotonE=0;
		int count_PhotonTheta=0;
		int count_MagP3Proton=0;
		int count_zCutmin=0;
		int count_RProton=0;
		int count_MissingMassSquared=0;
		int count_dEdxCDCProton=0;
		int count_insideEllipse=0;
		int count_allGeneralCutsPassed=0;
		int count_allGeneralCutsPassedPlusTracked=0;
		
		// location cuts
		bool inBCAL=true;
		bool inTOF=true;
		bool inFCAL=true;
		bool inSTART=true;
		bool inSYS_NULL=true;
		
		bool medBool[5]={true,true,true,true,true};

		// baryon cuts
		bool mMandelstamT_mdelta=true;
		bool mMandelstamT_mdelta_petaProton=true;
		bool mMandelstamT_mdelta_pVanHove=true;
		bool mDelta=true;

		
		// Various combinations of cuts, the majority of them will be used just for a few histograms like when showing unused energy graph we will use mUE which
		// removes the UE cut from allGeneralCutsPassed. m prefix basically stands for minus
		bool allGeneralCutsPassed=true;
		bool mMandelstamT=true;
		//bool pDiffCL=true; 
		bool pDiffUE=true; 
		bool mRProton=true;
		bool mRProtonZMin = true; 
		bool mdEdxCDC = true;
		bool mZMin = true;
		bool mMagP3 = true;
		bool mPhotonE = true;
		bool mPhotonTheta = true;
		bool mdij3 = true;
		bool mUE = true;
		bool mUEChiSq = true;
		bool mChiSq = true;
		bool mMMSq = true;
		// mEllipseRY contains both the red and yellow regions only.
		bool mEllipse = true;
		bool mEllipse_pre = true;
		bool mEllipse_pre_tAll = true;
		bool ptGT1=true;
		bool ptLT05=true;
		bool ptGT05LT1=true;
		bool mEllipse_pre_tGT1 = true;
		bool mEllipse_pre_tGT05LT1 = true;
		bool mEllipse_pre_tLT05 = true;
		bool mEllipse_pre_tAll_delta = true;
		bool mEllipseUE = true;
		bool mEllipseUE_pre = true;
		bool mEllipseUEChiSq = true;
		bool mEllipseUEChiSq_pre = true;
		bool mEllipseChiSq = true;
		bool mEllipseChiSq_pre = true;
		bool pMPi0P14=true;
		bool mMPi0P14_ellipse=true;
		bool baseCuts=true;
                bool noCut=true;

                bool dzRP=true;
                bool dzR=true;
		bool mEllipse_pre_pi0BCAL = true;
		bool mEllipse_pre_pi0FCAL = true;
		bool mEllipse_pre_pi0SPLIT = true;
		bool mEllipse_pre_etaBCAL = true;
		bool mEllipse_pre_etaFCAL = true;
		bool mEllipse_pre_etaSPLIT = true;
                bool allGen_barybkg=true;
                bool allGen_vanHove=true;
		bool cutsToApply = true; 
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////// ********* DEFINING ALL HISTOGRAMS ****************/////////////////////////////////////////////
		
		std::string pi0BinRange[3] = {"200", "0.05","0.25"} ;
		std::string etaBinRange[3] = {"300", "0.25","0.85"} ;
		int id;
       		int id_noCutDij3=0;
        	std::vector<std::vector<int>> vec_group_ids;
		std::string groupNames[16] = {"beam", "p1", "ph12b1", "ph1234", "ph12", "ph34", "ph34b1", "entire combo", "ph12p1", "ph34p1", "ph1234p1", "FCAL Pairs", "phN", "BCAL Pairs", "ph13", "ph24"};
        	std::vector<int> group_ids;
        	int groupVec_it; // this isnt used here but will be used later to denote the graphs we are filling.

		// this will index the different uniqueness tracking groups. In this updated uniqueness tracking method, we remove the need to have multiple uniquness tracking for each cut combination. So this value never changes
		//int idxCut;


		// Now lets make the array of histograms to fill
                TH1F* dHist_all1DHists[800];
                TH2F* dHist_all2DHists[800];

	ClassDef(DSelector_ver20, 0);
};

void DSelector_ver20::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
        dDecayingPi0Wrapper = dStep1Wrapper->Get_InitialParticle();
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

	//Step 2
	dStep2Wrapper = dComboWrapper->Get_ParticleComboStep(2);
        dDecayingEtaWrapper = dStep2Wrapper->Get_InitialParticle();
	dPhoton3Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(0));
	dPhoton4Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep2Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_ver20_h


