#ifndef DSelector_thrown_7_17_14_h
#define DSelector_thrown_7_17_14_h

#include <iostream>

#include "DSelector/DSelector.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_thrown_7_17_14 : public DSelector
{
	public:

		DSelector_thrown_7_17_14(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_thrown_7_17_14(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO
                int equal=0;
                int notEqual=0;
                int notPass=0;
                int pass=0;
		int beamEinf=0;
		int beamEfin=0;
                vector<Int_t> vecEqual;
                vector<Int_t> vecNotEqual;
                vector<Int_t> vecNotPassCuts;
                vector<Int_t> vecHasSecondaryDaughters;
                vector<Int_t> vecNumThrownNot7;
                Int_t eventIdx=0;
                TH1I* dHist_PID;
                TH1I* dHist_NumThrown;
		TH1F *mandelstam_tpAll; 
		TH1F *mandelstam_tpLT1 ;
		TH1F *mandelstam_tpLT06;
		TH1F *mandelstam_tpGT05LT1;
		TH2F *dHist_cosThetaVsMass_tpAll;
		TH2F *dHist_cosThetaVsMass_tpLT1;
		TH2F *dHist_cosThetaVsMass_tpLT06;
		TH2F *dHist_cosThetaVsMass_tpGT05LT1;
		TH2F *dHist_phiVsMass;
		TH1F *dHist_phi;
		TH1F *dHist_cosTheta;
		TH1F *dHist_beamE;

		bool pBeamE[12];
		bool pBeamE8to9GeV;
		TH1F *dHist_pi0eta1DBeam[12];
		TH1F *dHist_pi0eta1D;
		TH1F *dHist_phi8GeVPlus;
		TH1F *dHist_cosTheta8GeVPlus;


		int maxevent;
		int ievent;
		set<Int_t> showOutput =  {6, 82,188};

	ClassDef(DSelector_thrown_7_17_14, 0);
};

#endif // DSelector_thrown_7_17_14_h
