#!/bin/bash

numBins=$(grep "NUMBER_BINS" fit.cfg | cut -d "=" -f2)
#since bash is inclusive in the loop 
numBins=$(($numBins-1))
cwd=$(pwd)
fitDir=$(grep "FIT_DIR" fit.cfg | cut -d "=" -f2)
# need to remove the single quotes that come with just using grep
fitDir=${fitDir//\'/}
# remove the extra white space
fitDir=${fitDir//\ /}

echo "current working directory: "$cwd
echo "fitDir :"$fitDir

#for i in {0..$numBins};
for ((i=0; i<=$numBins; i++));
do
	echo "making directory: $cwd/$fitDir/bin_$i/seed-analysis/"
	mkdir -p $cwd/$fitDir/bin_$i/seed-analysis/
	mv $cwd/$fitDir/bin_$i/*seed*.* $cwd/$fitDir/bin_$i/seed-analysis/
	# moving the fit and log files of the full fit back to their original place since we have moved all the seed files already
	mv $cwd/$fitDir/bin_$i/bin_${i}-full.fit $cwd/$fitDir/bin_$i/bin_${i}.fit
	mv $cwd/$fitDir/bin_$i/bin_${i}-full.log $cwd/$fitDir/bin_$i/bin_${i}.log

done
