rm -rf EtaPi_fit
./divideData.pl
./fullFit_spawn.sh
./bootFit_spawn.sh
python runPlotEtaPiDeltas.py 2 1
python overlayBins.py 0 "S0+;S0+_D0+;S0+_D0+_D1+;S0+_D0+_D1+_D2+;P1+;"
python overlayBins.py 1 "S0+;S0+_D0+;S0+_D0+_D1+;S0+_D0+_D1+_D2+;P1+;"
./runBSanalysis.sh

