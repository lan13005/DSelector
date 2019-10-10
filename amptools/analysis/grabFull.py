from __future__ import division
import glob
import os
import sys
from loadCfg import loadCfg

def linspace(start, stop, n):
    h = (stop - start) / (n - 1)
    for i in range(n):
        yield start + h * i

def getNum(x):
        return int(x.split('_')[1])


cfg=loadCfg()
nBins=cfg["NUMBER_BINS"]
lowMass=cfg["LOW_MASS"]
upMass=cfg["UPP_MASS"]
fitDir=cfg["FIT_DIR"]
seedFile=cfg["SEED_FILE"]
saveAmps=cfg["SAVE_AMPS"]

os.chdir(fitDir)
myList = [bins for bins in glob.glob("bin*")]
myList = sorted( myList, key = getNum )
print myList
masses=list(linspace(lowMass,upMass,nBins))

counter = 0;

print "Writing to "+saveAmps+".txt"
with open(saveAmps+".txt", "w+") as writeFile:
	for mass in masses:
		binN = myList[counter]
		print "opening "+binN+"/"+seedFile
		with open(binN+"/"+seedFile) as inFile:
			writeFile.write(str(mass))
			for line in inFile:
				numSplit = len(line.split(' '))-1
				if numSplit == 4:
					initialize, amplitude, cartesian, real, im = ( item.strip() for item in line.split(' ', numSplit))
				if numSplit == 5:
					initialize, amplitude, cartesian, real, im, realTag = ( item.strip() for item in line.split(' ', numSplit))
				writeFile.write(" "+real+" "+im)
			writeFile.write('\n')
			counter+=1
	inFile.close()
writeFile.close()
