#!/bin/bash

############################### README ###################################
# We have to set up the DSelector and runDSelector properly. 
# 1. There is some lines that which should be uncommented when running over just the 
# polarization specific data sets in the runDSelector file. 
# 2. In the DSelector degAngle is set to some deg000,.., degAll. degAll actually doesnt do anythong
# but rather just something different to set the polarization of the beam to 
# be nonsensical so the histograms that have values to fill that depend on the
# polarization is not filled properly. 
# 3. remember to comment out the current line in the runDSelector so you don't load extra files
# 4. the root files locations split among the various polarizations should be in the
# separateRuns folder with the names padded by 0's to have 3 digits in total!

#########################################################################


# sed/// as / as a delimeter but can take anything incase our replacement has that delimeter in it.
echo ""
echo ""
echo "Make sure runDSelector_7_17_14.C is setup properly with deg000 to start"
echo "with the degAngle related code. If ready ... type anything"
read noUse


sed -i 's@degAngle="deg*"@degAngle="deg000"@g' runDSelector_7_17_14.C
sed -i 's@polarization="deg*"@polarization="deg000"@g' DSelector_ver20.C
sed -i 's@degAngle="deg*"@degAngle="deg000"@g' DSelector_ver20.C
root -l -q -b runDSelector_7_17_14.C 

sed -i 's@degAngle="deg000"@degAngle="deg045"@g' runDSelector_7_17_14.C
sed -i 's@polarization="deg000"@polarization="deg045"@g' DSelector_ver20.C
sed -i 's@degAngle="deg000"@degAngle="deg045"@g' DSelector_ver20.C
root -l -q -b runDSelector_7_17_14.C 

sed -i 's@degAngle="deg045"@degAngle="deg090"@g' runDSelector_7_17_14.C
sed -i 's@polarization="deg045"@polarization="deg090"@g' DSelector_ver20.C
sed -i 's@degAngle="deg045"@degAngle="deg090"@g' DSelector_ver20.C
root -l -q -b runDSelector_7_17_14.C 

sed -i 's@degAngle="deg090"@degAngle="deg135"@g' runDSelector_7_17_14.C
sed -i 's@polarization="deg090"@polarization="deg135"@g' DSelector_ver20.C
sed -i 's@degAngle="deg090"@degAngle="deg135"@g' DSelector_ver20.C
root -l -q -b runDSelector_7_17_14.C 

sed -i 's@degAngle="deg135"@degAngle="degAMO"@g' runDSelector_7_17_14.C
sed -i 's@polarization="deg135"@polarization="degAMO"@g' DSelector_ver20.C
sed -i 's@degAngle="deg135"@degAngle="degAMO"@g' DSelector_ver20.C
root -l -q -b runDSelector_7_17_14.C 

sed -i 's@degAngle="degAMO"@degAngle="degALL"@g' runDSelector_7_17_14.C
sed -i 's@polarization="degAMO"@polarization="degALL"@g' DSelector_ver20.C
sed -i 's@degAngle="degAMO"@degAngle="degALL"@g' DSelector_ver20.C
root -l -q -b runDSelector_7_17_14.C 

sed -i 's@degAngle="degALL"@degAngle="deg000"@g' runDSelector_7_17_14.C
sed -i 's@polarization="degALL"@polarization="deg000"@g' DSelector_ver20.C
sed -i 's@degAngle="degALL"@degAngle="deg000"@g' DSelector_ver20.C

