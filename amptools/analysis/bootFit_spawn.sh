#!/bin/bash


numBins=$(grep "NUMBER_BINS" fit.cfg | cut -d "=" -f2 | tr -d -c 0-9)
processSpawned=$(grep "PROCESSES_TO_SPAWN" fit.cfg | cut -d "=" -f2 | tr -d -c 0-9) # so that each one takes care of 5 bins
binsPerProc=$(python -c "from math import ceil; print int(ceil(1.0*$numBins/$processSpawned))")


echo "numBins, processSpawned, binsPerProc: $numBins, $processSpawned, $binsPerProc"
rm -rf fitLogs/boot/
mkdir -p fitLogs/boot/
#python removeSeedFiles.py

# the final element of the set should be numBins-1
# not that for i in {0,$iters} would not work since variables are not substituted in curly braces
#for ((i==0; i < $processSpawned; i++));
for ((iProcess=0; iProcess < $processSpawned; iProcess++));
do
        echo "iProcess: $iProcess"
	# drive the fit with these bin ranges and dump the output to a corresponding log file.
	# NOTE: (sleep 2; sleep 3) & 
	# vs    sleep 2 & sleep 3 &
	# where the frist is sequential processing of sleep in the background (using &)
	# and the second runs both sleep operations in parallel in the background
	echo "python bootFit.py $iProcess $binsPerProc $processSpawned > fitLogs/boot/driveFit_$iProcess.log &"
        python3 bootFit.py $iProcess $binsPerProc $processSpawned > fitLogs/boot/driveFit_$iProcess.log &
	echo "------------------------------"
done
wait
#echo "**The bin ranges, i.e. 0 5 will run on bin_0 to bin_4**"
#echo "** Below are the bins that did not converge! ** "
#grep 'not\ converged' fitLogs/boot/driveFit_*
#grep 'maxIter times' fitLogs/boot/driveFit_*
#echo "----------------------------"
#echo "** Below are files that had to try more than once to converge and the total number of tries attemped ** "
#grep -c "num_lines: 0" fitLogs/boot/driveFit* | awk -F: '$NF+0 > 0'
#echo "----------------------------"

echo "DONE!"
