//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 18 12:06:37 2019 by ROOT version 6.08/06
// from TTree Thrown_Tree/Thrown_Tree
// found on file: thrownNotAmptoolsReady_pi0eta_nameAffix.root
//////////////////////////////////////////////////////////

#ifndef split_gen_h
#define split_gen_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"
#include "TClonesArray.h"

class split_gen {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          RunNumber;
   ULong64_t       EventNumber;
   Int_t           ThrownBeam__PID;
   TLorentzVector  *ThrownBeam__X4;
   TLorentzVector  *ThrownBeam__P4;
   Float_t         ThrownBeam__GeneratedEnergy;
   ULong64_t       NumPIDThrown_FinalState;
   ULong64_t       PIDThrown_Decaying;
   Float_t         MCWeight;
  //NumFinalState should be int and not uint for split_mass
   Int_t          NumThrown;
   Int_t           Thrown__ParentIndex[16];   //[NumThrown]
   Int_t           Thrown__PID[16];   //[NumThrown]
   TClonesArray    *Thrown__X4;
   TClonesArray    *Thrown__P4;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_ThrownBeam__PID;   //!
   TBranch        *b_ThrownBeam__X4;   //!
   TBranch        *b_ThrownBeam__P4;   //!
   TBranch        *b_ThrownBeam__GeneratedEnergy;   //!
   TBranch        *b_NumPIDThrown_FinalState;   //!
   TBranch        *b_PIDThrown_Decaying;   //!
   TBranch        *b_MCWeight;   //!
   TBranch        *b_NumThrown;   //!
   TBranch        *b_Thrown__ParentIndex;   //!
   TBranch        *b_Thrown__PID;   //!
   TBranch        *b_Thrown__X4;   //!
   TBranch        *b_Thrown__P4;   //!

   split_gen(TTree *tree=0);
   virtual ~split_gen();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


   Float_t m_e[3];
   Float_t m_px[3];
   Float_t m_py[3];
   Float_t m_pz[3];
   Int_t m_PID[3];
   Float_t m_eBeam;
   Float_t m_pxBeam;
   Float_t m_pyBeam;
   Float_t m_pzBeam;
   Float_t m_weight;
   Float_t m_TargetMass;
   Int_t m_nPart;

   int counter;

   TTree *m_OutTree;
   TFile *outFile;


};

#endif

#ifdef split_gen_cxx
split_gen::split_gen(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("amptools4Rebecca/deg000_gen_2017_treeFlat_DSelector.root");
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("thrownNotAmptoolsReady_a0a2_nameAffix.root");
      cout << "File input: amptools4Rebecca/deg000_gen_2017_treeFlat_DSelector.root" << endl;
      if (!f || !f->IsOpen()) {
         f = new TFile("amptools4Rebecca/deg000_gen_2017_treeFlat_DSelector.root");
      }
      f->GetObject("Thrown_Tree",tree);

   }
   Init(tree);
}

split_gen::~split_gen()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t split_gen::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t split_gen::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void split_gen::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ThrownBeam__X4 = 0;
   ThrownBeam__P4 = 0;
   Thrown__X4 = 0;
   Thrown__P4 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("ThrownBeam__PID", &ThrownBeam__PID, &b_ThrownBeam__PID);
   fChain->SetBranchAddress("ThrownBeam__X4", &ThrownBeam__X4, &b_ThrownBeam__X4);
   fChain->SetBranchAddress("ThrownBeam__P4", &ThrownBeam__P4, &b_ThrownBeam__P4);
   fChain->SetBranchAddress("ThrownBeam__GeneratedEnergy", &ThrownBeam__GeneratedEnergy, &b_ThrownBeam__GeneratedEnergy);
   fChain->SetBranchAddress("NumPIDThrown_FinalState", &NumPIDThrown_FinalState, &b_NumPIDThrown_FinalState);
   fChain->SetBranchAddress("PIDThrown_Decaying", &PIDThrown_Decaying, &b_PIDThrown_Decaying);
   fChain->SetBranchAddress("MCWeight", &MCWeight, &b_MCWeight);
   fChain->SetBranchAddress("NumThrown", &NumThrown, &b_NumThrown);
   fChain->SetBranchAddress("Thrown__ParentIndex", Thrown__ParentIndex, &b_Thrown__ParentIndex);
   fChain->SetBranchAddress("Thrown__PID", Thrown__PID, &b_Thrown__PID);
   fChain->SetBranchAddress("Thrown__X4", &Thrown__X4, &b_Thrown__X4);
   fChain->SetBranchAddress("Thrown__P4", &Thrown__P4, &b_Thrown__P4);
   Notify();
}

Bool_t split_gen::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void split_gen::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t split_gen::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef split_gen_cxx
