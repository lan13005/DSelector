# This file grabs all the amplitudes from sensiN-param_init.cfg or seedN-param_init.cfg files inside all bins and changes the format abit so that it is easier to read
# plot-sensi.py (has to be done locally since it requires  numpy and matplotlib)

from __future__ import division
import glob
import sys
import os
from loadCfg import loadCfg 


cfg=loadCfg()
nBins=cfg["NUMBER_BINS"]
lowMass=cfg["LOW_MASS"]
upMass=cfg["UPP_MASS"]
nSeeds=cfg["NUMBER_SEEDS"]
saveAmps=cfg["SAVE_AMPS"]
#seedShift=cfg["SEED_SHIFT"]

print("nBins: "+str(nBins))
print("lowMass: "+str(lowMass))
print("upMass: "+str(upMass))
print("nSeeds: "+str(nSeeds))
print("saveAmps: "+str(saveAmps))
#print("seedShift: "+str(seedShift))

os.chdir(cfg["FIT_DIR"])

# sys.argv takes the first arg automatically as the file name
#if len(sys.argv) !=2 :
#	print("There should only be one argument excluding the defaulted one which is the filename")
#	print("1. sensi - will grab all param_init.cfg files with sensi prefixed")
#	print("2. boot - will grab all param_init.cfg files with seed prefixed")
#        raise ValueError('Incorrect number of arguments')
#if sys.argv[1] == 'sensi':
#        isSensi = 1
#elif sys.argv[1] == 'boot':
#        isSensi = 0
#else:
#        raise ValueError('argument should either be sensi or boot')
isSensi=0


# need this later to generate a mass list
def linspace(start, stop, n):
    if n == 1:
        yield stop
        return
    h = (stop - start) / (n - 1)
    for i in range(n):
        yield start + h * i

# unsorted bins when we use this method, so we have to order 
myList = [bins for bins in glob.glob("bin*")]

# need to use this function to sort the list
def getNum(x):
	return int(x.split('_')[1])

myList = sorted( myList, key = getNum )
print myList
masses = list(linspace(lowMass,upMass,nBins))

if isSensi:
	nSensi = 35
	nIters = nIters
	print("plot-sensi.py, grabData.py, ampSpaceScan.pl should share the same nSensi") 
else:
	#nSeeds = 200
	nIters = nSeeds
	print("grabData.py, driveFit.pl (bootstrap) should share the same nSeeds") 


analysis_folder = "/seed-analysis/"
binN = myList[0]
print(binN+analysis_folder+"*param_init.cfg")
allParamInit = glob.glob(binN+analysis_folder+"*param_init.cfg")
allSeedFiles =[seeds.split("/")[2] for seeds in allParamInit]
seeds =[seeds.split("/")[2].split("seed")[1].split("-")[0] for seeds in allParamInit]
print(allParamInit)
print(allSeedFiles)

for i in range(len(allSeedFiles)):
	#i = i+seedShift
	counter = 0; #will iterate the bins we get in myList as we iterate the masses list
	#if isSensi:
	#	tag = "sensi"+str(i)
	#else:
	#	tag = "seed"+str(i)
	print("Writing to "+saveAmps+"-"+seeds[i]+".txt")
	# w+ will write and create the file if needed 
	with open(saveAmps+"-seed"+seeds[i]+".txt", "w+") as writeFile:
		for mass in masses: 
			#if counter > 10:
			#	break
			binN = myList[counter]
			# these sensi param_init.cfg files should be directly in the bins and not in a subfolder
			print("Reading from: "+binN+analysis_folder+allSeedFiles[i])
			with open(binN+analysis_folder+allSeedFiles[i]) as inFile:
				writeFile.write(str(mass))
				for line in inFile:
					numSplit = len(line.split(' '))-1
					# we need these two cases since there are sometimes a real tag to fix the amplitude to be real
					if numSplit == 4:
						initialize, amplitude, cartesian, real, im = ( item.strip() for item in line.split(' ', numSplit))
					if numSplit == 5:
						initialize, amplitude, cartesian, real, im, realTag = ( item.strip() for item in line.split(' ', numSplit))
					writeFile.write(" "+real+" "+im)
				writeFile.write('\n')
			counter += 1
		inFile.close()
	writeFile.close()
