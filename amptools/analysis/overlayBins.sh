baseDir=$(pwd)
for ((iBin=0; iBin < 59; iBin++));
do
	#echo $(pwd)
	fitDir=${baseDir}/divideRoot
	#echo $fitDir
	binDir=${fitDir}/bin_$iBin
	cd $binDir	
	#echo $(pwd)
        #python bootFit.py $iProcess $binsPerProc $processSpawned > fitLogs/boot/driveFit_$iProcess.log &
	twopi_plotter_mom bin_${iBin}.fit
	cd $baseDir
done

root -l -b -q overlayBins.C
