from loadCfg import loadCfg
import glob 
import math 
from array import array
from ROOT import TCanvas, TGraphErrors
from ROOT import gROOT, TFile
import os
import sys
import shutil

saveFolder="amplitudes"

print(os.path.isdir(saveFolder))
if os.path.isdir(saveFolder):
	shutil.rmtree(saveFolder)
	print("removing %s folder!" % saveFolder)
	# os module only removes empty directories
	print("making %s folder!" % saveFolder)
	os.mkdir(saveFolder)
else:
	print("Error: %s folder not found!" % saveFolder)
	print("making %s folder!" % saveFolder)
	os.mkdir(saveFolder)


cfg = loadCfg()
ampNames = cfg["AMP_NAMES"]
nBins=cfg["NUMBER_BINS"]
lowMass=cfg["LOW_MASS"]
upMass=cfg["UPP_MASS"]
fitDir=cfg["FIT_DIR"]
saveAmps=cfg["SAVE_AMPS"]

# Not sure how to draw root objects with +/- prefix .... so we can just conver them to p/m
ampNamesOnlyLetters = []
for ampName in ampNames:
	ampNamesOnlyLetters.append(ampName.replace("+","p").replace("-","m"))

verbose=False

# need to use this function to sort the list
def getNum(x):
        return int(x.split('seed')[1].split('.')[0])



allAmpFiles = glob.glob(fitDir+"/"+saveAmps+"-seed*")
allAmpFiles = sorted( allAmpFiles, key = getNum )
if verbose:
	print(allAmpFiles)
nSeeds = len(allAmpFiles)
print("nSeeds: "+ str(nSeeds))

massStep = (upMass-lowMass)/nBins

# allAmps_MassSeed would contain all the cooridinates ( amps, mass, seed )
# this is so that we can loop over all the amps and mass combinations like allAmps_MassSeed[N][M] and calculate the avg and std of that 1D array easily
numAmps = int(len(ampNames))
allAmps_MassSeed = [[[0 for col in range(nSeeds)] for row in range(nBins)] for depth in range(numAmps)]

numEntries = 2*numAmps+1
print("Number of entries in each line of converged amps file: "+str(numEntries))
print("Corresponds to 1 mass value and 2 components for each ampltiude (real/imag)")

print("Getting amplitude squared values from full data...")
allAmps_Mass = [[0 for col in range(nBins)] for row in range(numAmps)]

with open(fitDir+"/"+saveAmps+".txt", "r") as ampFile:
	for massNum,massLine in enumerate(ampFile):
        	numSplit = len(massLine.split(' '))-1
		#print("numSplit: "+str(numSplit))
		# there should be numEntries-1 splits
		if (numSplit == numEntries-1):
        	#	mass, S0p_r, S0p_i, P0p_r, P0p_i, P1p_r, P1p_i, D0p_r, D0p_i, D1p_r, D1p_i, P1m_r, P1m_i, D1m_r, D1m_i = ( float(item.strip()) for item in massLine.split(' ', numSplit))
			# there will be numSplit/2 ampltiudes
			values = [float(item.strip()) for item in massLine.split(' ', numSplit)]
			if verbose:
				print(values)
			for ampNum in range(numAmps):
				amps = values[2*ampNum+1]**2 + values[2*ampNum+2]**2
				allAmps_Mass[ampNum][massNum] = float(amps)
				if verbose: 
					print("real part: "+str(values[2*ampNum+1]))
					print("imag part: "+str(values[2*ampNum+2]))
					print("amp sq: "+str(amps))
					print("------------------------------")
		else:
			raise ValueError("NOT THE CORRECT NUMBER OF ENTRIES!")
		if verbose: 
			print("getting next mass line!")

# This code is basically a loop of the above one so all debugging should be done with that one. 
print("bootstrapped amplitudes squared values are being calculated....")
print("since the seed files might come in at a random order we will create a dictionary to match the seedNum with the index in allAmps_MassSeed")
seedToIdx = {}
for idxSeed,ampFile in enumerate(allAmpFiles):
	seedNum = int(ampFile.split("seed")[1].split(".")[0])
	seedToIdx[seedNum] = idxSeed	
	sys.stdout.write("#")
	#print(ampFile+" is being loaded...")
	with open(ampFile, "r") as ampFile:
		for massNum,massLine in enumerate(ampFile):
	        	numSplit = len(massLine.split(' '))-1
			if (numSplit == numEntries-1):
				values = [float(item.strip()) for item in massLine.split(' ', numSplit)]
				for ampNum in range(numAmps):
					amps = values[2*ampNum+1]**2 + values[2*ampNum+2]**2
					print("ampNum, massNum, seedNum: {0}, {1}, {2}".format(ampNum,massNum,idxSeed))
					allAmps_MassSeed[ampNum][massNum][idxSeed] = amps
			else:
				raise ValueError("NOT THE CORRECT NUMBER OF ENTRIES!")
sys.stdout.write("\n")
sys.stdout.flush()

print("calculating the mean and population std (instead of the sample std)  of the amplitude squared values over all seeds")
allAmps_MassSeed_mean=[[0 for col in range(nBins)] for row in range(numAmps)]
allAmps_MassSeed_std=[[0 for col in range(nBins)] for row in range(numAmps)]
for ampNum in range(numAmps):
	for massNum in range(nBins):
		for seedNum in range(nSeeds):
			allAmps_MassSeed_mean[ampNum][massNum]+=allAmps_MassSeed[ampNum][massNum][seedNum]
		allAmps_MassSeed_mean[ampNum][massNum]/=nSeeds
		for seedNum in range(nSeeds):
			allAmps_MassSeed_std[ampNum][massNum]+=(allAmps_MassSeed[ampNum][massNum][seedNum]-allAmps_MassSeed_mean[ampNum][massNum])**2
		allAmps_MassSeed_std[ampNum][massNum]/=(nSeeds)
		allAmps_MassSeed_std[ampNum][massNum]=float(math.sqrt(allAmps_MassSeed_std[ampNum][massNum]))



# this is kind of the weird step since we do all this importing into python and calculating then we have to rebulild the arrays to be usuable with pyrootampFile=ROOT.TFile.Open("amplitudes.root")
ampFile=TFile.Open(saveFolder+"/amplitudes.root","RECREATE")
ampCanvas = TCanvas("ampCanvas", "amplitudes", 1200,600)

for ampNum in range(numAmps):
	amplitude = array( 'd' )
	masses = array( 'd' )
	amplitude_std = array( 'd' )
	mass_std = array( 'd' )
	for massNum in range(nBins):
		amplitude.append( allAmps_Mass[ampNum][massNum] )
		amplitude_std.append( allAmps_MassSeed_std[ampNum][massNum] )
		# masses start at lowMass up to upMass with step size massStep. We shift up to lowMass+massStep first which is the upper edge of the first bin. 
		masses.append( massStep*massNum+lowMass+massStep )	
		mass_std.append( 0.0 )
	
	gr = TGraphErrors( nBins, masses, amplitude, mass_std, amplitude_std )
	gr.SetMarkerColor( 9 )
	gr.SetMarkerStyle( 21 )
	gr.SetTitle( ampNames[ampNum] )
	gr.GetXaxis().SetTitle( 'Masses (GeV)' )
	gr.GetYaxis().SetTitle( 'Amplitude' )
	gr.Draw( 'AP' )
	ampCanvas.SaveAs(saveFolder+"/"+ampNamesOnlyLetters[ampNum]+".pdf")
	gr.Write(ampNamesOnlyLetters[ampNum])








