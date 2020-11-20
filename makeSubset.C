#include "/d/grid13/gluex/gluex_top/gluex_style.C"

void makeSubset() {
   //Get old file, old tree and set top branch address
   TFile *oldfile = new TFile("/d/grid15/ln16/pi0eta/092419/degALL_a0a2Test_treeFlat_DSelector.root");
   TTree *oldtree = (TTree*)oldfile->Get("degALL_a0a2Test_tree_flat");
   Long64_t nentries = oldtree->GetEntries();
   double Mpi0;
   double Meta;
   double Mpi0eta;
   double chiSq;
   double unusedEnergy;
   oldtree->SetBranchAddress("Mpi0",&Mpi0);
   oldtree->SetBranchAddress("Meta",&Meta);
   oldtree->SetBranchAddress("Mpi0eta",&Mpi0eta);
   oldtree->SetBranchAddress("chiSq",&chiSq);
   oldtree->SetBranchAddress("unusedEnergy",&unusedEnergy);

   //Create a new file + a clone of old tree in new file
   TFile *newfile = new TFile("/d/grid15/ln16/pi0eta/q-values/subset_a0a2.root","recreate");
   TTree *newtree = oldtree->CloneTree(0);

   for (Long64_t i=0;i<nentries; i++) {
        if (rand()%100<5){
            oldtree->GetEntry(i);
            newtree->Fill();
        }
   }
   newtree->Print();
   newtree->AutoSave();
   delete oldfile;
   delete newfile;
}
