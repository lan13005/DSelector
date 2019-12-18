echo "*****************************"
echo "SAVING OUTPUT TO: logTest.txt"
echo "*****************************"

#is_pi0eta=true
#
#if $is_pi0eta 
#then 
#	sed -i "s@pi0eta=false@pi0eta=true@g" DSelector_eta3pi0.h
#fi


sed -i "s@showOutput\ =\ false@showOutput\ =\ true@g" DSelector_eta3pi0.C
sed -i "s@//if(itersToRun@if(itersToRun@g" DSelector_eta3pi0.C
sed -i "s@//}//closes@}//closes@g" DSelector_eta3pi0.C
sed -i "s@useproof\ =\ 1@useproof\ =\ 0@g" runDSelector_7_17_14.C

rm -f logTest.txt
root -l -b -q runDSelector_7_17_14.C > logTest.txt

sed -i "s@showOutput\ =\ true@showOutput\ =\ false@g" DSelector_eta3pi0.C
sed -i "s@if(itersToRun@//if(itersToRun@g" DSelector_eta3pi0.C
sed -i "s@}//closes@//}//closes@g" DSelector_eta3pi0.C
sed -i "s@useproof\ =\ 0@useproof\ =\ 1@g" runDSelector_7_17_14.C


#if $is_pi0eta
#then 
#	sed -i "s@pi0eta=true@pi0eta=false@g" DSelector_eta3pi0.h
#fi


totalHists=$(grep '^Total\ Number' logTest.txt | cut -d " " -f 5)
availHists=$(grep 'dHist_all1DHists' DSelector_eta3pi0.h | cut -d "[" -f 2 | cut -d "]" -f 1)

echo "total number of histograms: $totalHists"
echo "total available histograms: $availHists"
echo "*********************************"
if [ "$totalHists" -ge "$availHists" ] 
then
	echo "NOT GOOD! MAKE AVAILHISTS LARGER!"
else
	echo "SIZE OF HIST VECTOR IS FINE!"
fi
