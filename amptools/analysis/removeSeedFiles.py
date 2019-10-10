import glob
import os
from loadCfg import loadCfg
import sys



cfg = loadCfg()
nSeeds=cfg['NUMBER_SEEDS']
nBins=cfg['NUMBER_BINS']
shiftSeeds=cfg['SEED_SHIFT']
singleBin=cfg["SINGLE_BIN"]

print("nSeeds: "+str(nSeeds))
print("shiftSeeds: "+str(shiftSeeds))

print(range(nSeeds))
print(shiftSeeds)


if singleBin in range(nBins): 
	for seedNum in range(shiftSeeds,shiftSeeds+nSeeds):
		#print("divideRoot/bin_"+str(singleBin)+"/*/*seed"+str(seedNum)+"*.*")
		fileSeeds = glob.glob("divideRoot/bin_"+str(singleBin)+"/*/*seed"+str(seedNum)+"*.*")
		#print("removing all files in all the bins that has to do with seeds!")
		for iFile in fileSeeds:
		        print("removing "+iFile)
		        os.remove(iFile)
else:
	for seedNum in range(shiftSeeds,shiftSeeds+nSeeds):
		fileSeeds = glob.glob("divideRoot/*/*seed"+str(seedNum)+"*.*")
		print("removing all files in all the bins that has to do with seeds!")
		for iFile in fileSeeds:
		        print("removing "+iFile)
		        os.remove(iFile)
	

