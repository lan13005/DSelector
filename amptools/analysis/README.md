# Current Propsed Method:
Detailed flow diagram in https://docs.google.com/presentation/d/1-b0ZIQYIp3QvBW7MuMyHOe1E7ZB5i-zVTKkM1gTNrMM/edit#slide=id.gb8f5457b7d_0_7

## Genearal Order
tree_to_amptools to make post-DSelector trees amptools-ready
fullFit_spawn.sh 
bootFit_spawn.sh
python overlayBins.py 
python runPlotEtaPiDeltas.py
python diagnoseBS.py
python plotIntensities.py

## Example full analysis shell script:
runAnalysis.sh
Note:
overlayBins.C uses etapi_plotter program which is not in halld_sim
