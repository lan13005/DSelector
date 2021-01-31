#!/bin/bash

rm -rf fitLogs/full
mkdir -p fitLogs/full
python fullFit.py  > fitLogs/full/driveFit.log 

#####################
# this older section runs the fit in various bins with multiple processes
# This turns out to not be so helpful since there are bin to bin correlations that we can use.
# i.e. using the converged amplitudes from the previous bin as the initialization for the current bin
#####################
#     numBins=$(grep "NUMBER_BINS" fit.cfg | cut -d "=" -f2)
#     processSpawned=$(grep "PROCESSES_TO_SPAWN" fit.cfg | cut -d "=" -f2) # so that each one takes care of 5 bins
#     binsPerProc=$(($numBins/$processSpawned))
#     
#     
#     echo "numBins, processSpawned, binsPerProc: $numBins, $processSpawned, $binsPerProc"
#     echo "---------------------------------------------------------"
#     rm -rf fitLogs/full
#     mkdir -p fitLogs/full
#     
#     # the final element of the set should be numBins-1
#     # not that for i in {0,$iters} would not work since variables are not substituted in curly braces
#     for ((iProcess=0; iProcess < $processSpawned; iProcess++));
#     do
#     	echo "iProcess: $iProcess"
#     	# drive the fit with these bin ranges and dump the output to a corresponding log file.
#     	# NOTE: (sleep 2; sleep 3) & 
#     	# vs    sleep 2 & sleep 3 &
#     	# where the frist is sequential processing of sleep in the background (using &)
#     	# and the second runs both sleep operations in parallel in the background
#     	python3 fullFit.py $iProcess $binsPerProc $processSpawned > fitLogs/full/driveFit_$iProcess.log &
#     done
#     wait
#     #echo "**The bin ranges, i.e. 0 5 will run on bin_0 to bin_4**"
#     #echo "** Below are the bins that did not converge! ** "
#     #grep 'maxIter times' fitLogs/full/driveFit_*
#     #echo "----------------------------"
#     #echo "** Below are files that had to try more than once to converge and the number of tries attemped ** "
#     #grep -c "num_lines: 0" fitLogs/full/driveFit* | awk -F: '$NF+0 > 0'
#     #echo "----------------------------"
#     #echo "** Below are files that had converged on first attempt** "
#     #grep -c "num_lines: 0" fitLogs/full/driveFit* | awk -F: '$NF+0 == 0'
#     #echo "----------------------------"

echo "DONE!"
