#define split_gen_cxx
#include "split_gen.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void split_gen::Loop()
{
//   In a ROOT session, you can do:
//      root> .L split_gen.C
//      root> split_gen t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   outFile = new TFile("thrown_a0a2_pi0eta_nameAffix.root", "RECREATE");
   m_OutTree = new TTree("Thrown_Tree", "kin2");

   m_OutTree->Branch("Weight", new Float_t, "Weight/F");
   m_OutTree->Branch("E_Beam", new Float_t, "E_Beam/F");
   m_OutTree->Branch("Px_Beam", new Float_t, "Px_Beam/F");
   m_OutTree->Branch("Py_Beam", new Float_t, "Py_Beam/F");
   m_OutTree->Branch("Pz_Beam", new Float_t, "Pz_Beam/F");
   m_OutTree->Branch("Target_Mass", new Float_t, "Target_Mass/F");
  //NumFinalState should be int and not uint for split_mass
   m_OutTree->Branch("NumFinalState", new Int_t, "NumFinalState/I");
   m_OutTree->Branch("PID_FinalState", new Int_t[3], "PID_FinalState[NumFinalState]/I");
   m_OutTree->Branch("E_FinalState", new Float_t[3], "E_FinalState[NumFinalState]/F");
   m_OutTree->Branch("Px_FinalState", new Float_t[3], "Px_FinalState[NumFinalState]/F");
   m_OutTree->Branch("Py_FinalState", new Float_t[3], "Py_FinalState[NumFinalState]/F");
   m_OutTree->Branch("Pz_FinalState", new Float_t[3], "Pz_FinalState[NumFinalState]/F");

   // weight will always be set to 1 for thrown data!
   m_OutTree->SetBranchAddress("Weight", &m_weight);
   m_OutTree->SetBranchAddress("E_Beam", &m_eBeam);
   m_OutTree->SetBranchAddress("Px_Beam", &m_pxBeam);
   m_OutTree->SetBranchAddress("Py_Beam", &m_pyBeam);
   m_OutTree->SetBranchAddress("Pz_Beam", &m_pzBeam);
   m_OutTree->SetBranchAddress("Target_Mass", &m_TargetMass);
   m_OutTree->SetBranchAddress("NumFinalState", &NumThrown);
   m_OutTree->SetBranchAddress("PID_FinalState", Thrown__PID);
   m_OutTree->SetBranchAddress("E_FinalState", m_e);
   m_OutTree->SetBranchAddress("Px_FinalState", m_px);
   m_OutTree->SetBranchAddress("Py_FinalState", m_py);
   m_OutTree->SetBranchAddress("Pz_FinalState", m_pz);

   m_OutTree->SetBranchAddress("NumFinalState", &m_nPart);
   m_nPart = 3;

   m_OutTree->SetBranchAddress("PID_FinalState", m_PID);
   m_PID[0] = 14; m_PID[1] = 7; m_PID[2] = 17;

   m_TargetMass=0.931494;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   double countProton = 0;
   double countPi0 = 0;
   double countEta = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) { cout << "breaking since ientry<0" << endl; break; } 
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
	//cout << "NumThrown: " << NumThrown << endl;
	countProton=0;	
	countPi0=0;	
	countEta=0;	
	cout << jentry << endl;
        for (int particleN=0; particleN<NumThrown; ++particleN){
                TLorentzVector *single_P4 = (TLorentzVector *)Thrown__P4->At(particleN);
		if (Thrown__PID[particleN]==221 && Thrown__ParentIndex[particleN]==-1) {
			++countEta;
                	m_e[2] = single_P4->E();
                	m_px[2] = single_P4->Px();
                	m_py[2] = single_P4->Py();
                	m_pz[2] = single_P4->Pz();
		}
		else if (Thrown__PID[particleN]==111 && Thrown__ParentIndex[particleN]==-1) {
			++countPi0;
                	m_e[1] = single_P4->E();
                	m_px[1] = single_P4->Px();
                	m_py[1] = single_P4->Py();
                	m_pz[1] = single_P4->Pz();
		}
		else if (Thrown__PID[particleN]==2212 && Thrown__ParentIndex[particleN]==-1) {
			++countProton;
                	m_e[0] = single_P4->E();
                	m_px[0] = single_P4->Px();
                	m_py[0] = single_P4->Py();
                	m_pz[0] = single_P4->Pz();
		}
                //cout << "particleN, PID, Mass: "<< particleN << ", " << Thrown__PID[particleN]<<", "<<single_P4->M()<<endl;
		//cout << "E, px, py, pz: " << m_e[particleN] << ", " << m_px[particleN] << ", " << m_py[particleN] << ", " << m_pz[particleN] << endl;
        }
	if (countProton!=1 && countPi0!=1 && countEta !=1){
		cout << "THERE MUST BE 1 PROTON, ETA, PI0!!!!" << endl;
		exit (EXIT_FAILURE);
	}
        m_eBeam = ThrownBeam__P4->E();
        m_pxBeam = ThrownBeam__P4->Px();
        m_pyBeam = ThrownBeam__P4->Py();
        m_pzBeam = ThrownBeam__P4->Pz();
        m_weight = 1;
        //cout << "beam E, px, py, pz: " << m_eBeam << ", "<< m_pxBeam << ", "<< m_pyBeam << ", "<< m_pzBeam << endl;
	//cout << "Weight, target_mass: " << m_weight << ", " << m_TargetMass << endl;
       //cout << "Filled output tree" << endl;
       m_OutTree->Fill();


   }
   cout << "Completed loop: nbytes =" << nbytes << " nentries=" << nentries << endl;
   m_OutTree->Write();
   outFile->Close();
   cout << " nentries=" << nentries << " nb=" << nb << " nbytes=" << nbytes << endl;
}
