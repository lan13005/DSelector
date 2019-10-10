# For some reason some of the spawned proccesses of bootfit_spawn does not finish and ends prematurely. We can use this code to check what the last seed was attempted before it ended

numSeeds=$(grep "NUMBER_SEEDS" fit.cfg | cut -d "=" -f2)
numSeeds=$((numSeeds-1))
for ((i=0; i < 60; i++));
do
	counts=`ls divideRoot/bin_$i/seed-analysis/ | grep param | awk -F"-" '{print $1}' | awk -F"seed" '{print $2}' | sort -n | tail -n1`
	if [ -z "$counts" ] 
	then
		echo "FOR BIN $i: I THNINK YOU FORGOT TO MOVE ALL THE SEED FILES INTO THE SEED-ANALYSIS FOLDER!"
	fi
	if [ "$counts" -lt "$numSeeds" ]; then
		#echo "For fitLogs/boot/driveFit_$i\.log"
		#grep "fit -c" fitLogs/boot/driveFit_$i\.log | tail -1 | awk -F"seed" '{print $2}' | awk -F"." '{print $1}'
		echo "For bin_$i: $counts"
	fi
done
