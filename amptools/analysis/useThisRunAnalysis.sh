# runs the fit on the data
./divideData.pl
./fullFit_spawn.sh

# makes the individual diagnostic plots in each bin and overlays them
# python plotDiagnostics.py
# root -l -b -q plotDiagnostics.C

# outputs the amplitudes with the correct noramlization I think into a file. plot.C will plot them
plot_etapi_delta -o amplitudes.txt
root -l -b -q plot.C

