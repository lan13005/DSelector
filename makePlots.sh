fileLoc='"degALL_data_2017_allCuts_hists_DSelector.root"'

rm -r newGraphs
mkdir -p newGraphs

#root -l -b -q "makeGraphs.C($fileLoc)" 
#root -l -b -q "makeBaryonPlots.C($fileLoc)" 
#root -l -b -q "makeMpi0etaBinnedT.C($fileLoc)" 
root -l -b -q "makeRectSBGraphs.C($fileLoc)" 
#root -l -b -q "makeDrawDeckLines.C($fileLoc)" 
#root -l -b -q "makeSpecialtyPlots.C($fileLoc)" 
#root -l -b -q "makeMassVaryChiSq.C($fileLoc)"
