#include "split_gen.C"
void call_split_thrown ()
{
// issue the tree->Loop() from the command line.
//
  // gROOT->ProcessLine(".L MakeAmpToolsFlat_mcthrown.C");
  //gROOT->LoadMacro("MakeAmpToolsFlat_mcthrown.C");
 split_gen t;
 t.Loop();
}
