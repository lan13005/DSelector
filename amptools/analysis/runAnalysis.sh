rm -rf EtaPi_fit
./divideData.pl
./fullFit_spawn.sh
./bootFit_spawn.sh
python runPlotEtaPiDeltas.py 2 1
python overlayBins.py 0 "S0+;S0+_D0+;S0+_D0+_D1+;S0+_D0+_D1+_D2+;P1+;"
python overlayBins.py 1 "S0+;S0+_D0+;S0+_D0+_D1+;S0+_D0+_D1+_D2+;P1+;"
python diagnoseBS.py H0_10 plot_etapi_delta_results project_moments_polarized uncert " " True
python diagnoseBS.py S0+ plot_etapi_delta_results plot_etapi_delta err "\t" True
python plotIntensities.py True
python plotIntensities.py False
