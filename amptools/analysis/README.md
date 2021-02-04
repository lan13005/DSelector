# Current Propsed Method:
Detailed flow diagram:  
https://docs.google.com/presentation/d/1-b0ZIQYIp3QvBW7MuMyHOe1E7ZB5i-zVTKkM1gTNrMM/edit#slide=id.gb8f5457b7d_0_7

## General Order
1. tree_to_amptools to make post-DSelector trees amptools-ready
2. divideData.pl -- make sure the lowmass,upmass,nbins are the same here and in fit.cfg

Before continuing try and run "fit -c cfgFileName" in one of the bins first. It will give errors about any configuration error that the following
processes will ignore!

3. fullFit_spawn.sh 
4. bootFit_spawn.sh
5. python overlayBins.py -- to run this you need to put the contents of halld_sim folder of the this repo into halld_sim and compile it
6. python runPlotEtaPiDeltas.py -- need to update plote_etapi_delta and project_moments to your own waveset and configuration
7. python diagnoseBS.py 
8. python plotIntensities.py

## Example full analysis shell script:
runAnalysis.sh  

