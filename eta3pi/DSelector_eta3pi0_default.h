#ifndef DSelector_eta3pi_h
#define DSelector_eta3pi_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"

class DSelector_eta3pi : public DSelector
{
	public:

		DSelector_eta3pi(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_eta3pi(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		// ANALYZE CUT ACTIONS
		// // Automatically makes mass histograms where one cut is missing
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		//Step 1
		DParticleComboStep* dStep1Wrapper;
		DNeutralParticleHypothesis* dPhoton1Wrapper;
		DNeutralParticleHypothesis* dPhoton2Wrapper;

		//Step 2
		DParticleComboStep* dStep2Wrapper;

		//Step 3
		DParticleComboStep* dStep3Wrapper;
		DNeutralParticleHypothesis* dPhoton3Wrapper;
		DNeutralParticleHypothesis* dPhoton4Wrapper;

		//Step 4
		DParticleComboStep* dStep4Wrapper;
		DNeutralParticleHypothesis* dPhoton5Wrapper;
		DNeutralParticleHypothesis* dPhoton6Wrapper;

		//Step 5
		DParticleComboStep* dStep5Wrapper;
		DNeutralParticleHypothesis* dPhoton7Wrapper;
		DNeutralParticleHypothesis* dPhoton8Wrapper;

		// DEFINE YOUR HISTOGRAMS HERE
		// EXAMPLES:
		TH1I* dHist_MissingMassSquared;
		TH1I* dHist_BeamEnergy;

	ClassDef(DSelector_eta3pi, 0);
};

void DSelector_eta3pi::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));

	//Step 1
	dStep1Wrapper = dComboWrapper->Get_ParticleComboStep(1);
	dPhoton1Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));

	//Step 2
	dStep2Wrapper = dComboWrapper->Get_ParticleComboStep(2);

	//Step 3
	dStep3Wrapper = dComboWrapper->Get_ParticleComboStep(3);
	dPhoton3Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep3Wrapper->Get_FinalParticle(0));
	dPhoton4Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep3Wrapper->Get_FinalParticle(1));

	//Step 4
	dStep4Wrapper = dComboWrapper->Get_ParticleComboStep(4);
	dPhoton5Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep4Wrapper->Get_FinalParticle(0));
	dPhoton6Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep4Wrapper->Get_FinalParticle(1));

	//Step 5
	dStep5Wrapper = dComboWrapper->Get_ParticleComboStep(5);
	dPhoton7Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep5Wrapper->Get_FinalParticle(0));
	dPhoton8Wrapper = static_cast<DNeutralParticleHypothesis*>(dStep5Wrapper->Get_FinalParticle(1));
}

#endif // DSelector_eta3pi_h
