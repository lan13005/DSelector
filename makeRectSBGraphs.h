#ifndef MAKERECTSBGRAPHS_H
#include <TFile.h>

TH1F combineSBRegions(TFile* dataFile, string reg1Name, string reg2Name, string reg3Name, string reg4Name, string histName, string fileSaveName);

#define MAKERECTSBGRAPHS_H
#endif

