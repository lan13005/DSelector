#include "split_selected_in_t.C"
void call_split_selected()
{
// issue the tree->Loop() from the command line.
//
  // gROOT->ProcessLine(".L MakeAmpToolsFlat.C");
  //gROOT->LoadMacro("MakeAmpToolsFlat.C");
 split_selected_in_t t;
 t.Loop();
}
