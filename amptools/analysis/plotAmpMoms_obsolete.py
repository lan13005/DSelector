from loadCfg import loadCfg
import glob 
import math 
from math import sqrt
from array import array
from ROOT import TCanvas, TGraphErrors, TGraph, TH1D, TLegend
from ROOT import gROOT, TFile
import os
import sys
import shutil

saveFolder="amplitudes"

makeMoments=True

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

verbose=True

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

momNames=["H00","H10","H11","H20","H21","H22","H30","H31","H32","H40","H41","H42"]
momAsym0=["H00","H20","H22","H40","H42","H44"]
momAsym1=["H00","H20","H22","H40","H42","H44"]
numMomsAsym0 = len(momAsym0)
numMoms=len(momNames)
# complexAmps = complex amplitudes matrix
complexAmps = [[0 for col in range(nBins)] for row in range(numAmps)]
with open(fitDir+"/"+saveAmps+".txt", "r") as ampFile:
	for massNum,massLine in enumerate(ampFile):
        	numSplit = len(massLine.split(' '))-1
		#print("numSplit: "+str(numSplit))
		# there should be numEntries-1 splits
		if (numSplit == numEntries-1):
		# 	The general form would take the following form 
        	#	mass, S0_r, S0_i, P0_r, P0_i, P1p_r, P1p_i, D0_r, D0_i, D1p_r, D1p_i, P1m_r, P1m_i, D1m_r, D1m_i = ( float(item.strip()) for item in massLine.split(' ', numSplit))
			# there will be numSplit/2 ampltiudes
			values = [float(item.strip()) for item in massLine.split(' ', numSplit)]
			if verbose:
				print(values)
			for ampNum in range(numAmps):
				amps = values[2*ampNum+1]**2 + values[2*ampNum+2]**2
				allAmps_Mass[ampNum][massNum] = float(amps)
				# each complex amplitude like S0+ would have a real and imaginary part corresponding to some elements of the values array
				complexAmps[ampNum][massNum] = values[2*ampNum+1]+values[2*ampNum+2]*1j
				if verbose: 
					print("Amplitue: "+ampNamesOnlyLetters[ampNum])
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
					#print("ampNum, massNum, seedNum: {0}, {1}, {2}".format(ampNum,massNum,idxSeed))
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
print("Rebuilding the python arrays to use the array module which can then be used by pyroot")
amplitudesFile=TFile.Open(saveFolder+"/amplitudes.root","RECREATE")
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

print("Calculating the summed amplitudes with errors combined in quadruture")
ampTot = array ( 'd' )
ampTot_std = array ( 'd' ) 
for massNum in range(nBins):
	ampTot.append(0)
	ampTot_std.append(0)
	for ampNum in range(numAmps):
		ampTot[massNum] += allAmps_Mass[ampNum][massNum]
		ampTot_std[massNum] += allAmps_MassSeed_std[ampNum][massNum]**2
	ampTot_std[massNum] = sqrt(ampTot_std[massNum])
gr = TGraphErrors( nBins, masses, ampTot, mass_std, ampTot_std )
gr.SetMarkerColor( 9 )
gr.SetMarkerStyle( 21 )
gr.SetTitle( "Summed Amps" )
gr.GetXaxis().SetTitle( 'Masses (GeV)' )
gr.GetYaxis().SetTitle( 'Amplitudes Summed' )
gr.Draw( 'AP' )
ampCanvas.SaveAs(saveFolder+"/summedAmps.pdf")
gr.Write("summedAmps")


print("Begin calculating the moments!")
# This step trades readibility for generality but it is OK if we keep the order of amps the same
vS0m,vP0m,vP1m,vD0m,vD1m,vP1p,vD1p=complexAmps
allMoms_Mass = [[0 for col in range(nBins)] for row in range(numMoms)]
allMoms_Mass_std = [[0 for col in range(nBins)] for row in range(numMoms)]
if makeMoments:
	for massNum in range(nBins):
		# just for the H00 moment where we can loop through all the complex amps
		for ampNum in range(numAmps):
			allMoms_Mass[0][massNum] += complexAmps[ampNum][massNum]*complexAmps[ampNum][massNum].conjugate()
			allMoms_Mass[0][massNum] = allMoms_Mass[0][massNum].real
		# ********* We will neglect which waveset the M=0 waves belongs to since it seems from the eq that the moments dont care *************
		S0 = vS0m[massNum]
		c_S0p = S0.conjugate()
		P0 = vP0m[massNum]
		c_P0 = P0.conjugate()
		P1m = vP1m[massNum]
		c_P1m = P1m.conjugate()
		D0 = vD0m[massNum]
		c_D0 = D0.conjugate()
		D1m = vD1m[massNum]
		c_D1m = D1m.conjugate()
		P1p = vP1p[massNum]
		c_P1p = P1p.conjugate()
		D1p = vD1p[massNum]
		c_D1p = D1p.conjugate()
		# calculating the moments 
		#momNames=["H00","H10","H11","H20","H21","H22","H30","H31","H32","H40","H41","H42"]
        	allMoms_Mass[1][massNum] = 2./ sqrt(3.)*(S0*c_P0).real+4./sqrt(15.)*(P0*c_D0).real+2./sqrt(5.)*(P1m*c_D1m).real+2./sqrt(5.)*(P1p*c_D1p).real
        	allMoms_Mass[2][massNum] = 2./ sqrt(6.)*(S0*c_P1m).real+2./sqrt(10.)*(P0*c_D1m).real-2./sqrt(30.)*(P1m*c_D0).real
        	allMoms_Mass[3][massNum] = 2./ sqrt(5.)*(S0*c_D0).real+2./5.*(P0*c_P0).real-1./5.*(P1m*c_P1m).real-1./5.*(P1p*c_P1p).real + 2./7.*(D0*c_D0).real+1./7.*(D1m*c_D1m).real+1./7.*(D1p*c_D1p).real
        	allMoms_Mass[4][massNum] = 2./sqrt(10.)*(S0*c_D1m).real+2./5.*sqrt(3./2.)*(P0*c_P1m).real+2./7./sqrt(2.)*(D0*c_D1m).real
        	allMoms_Mass[5][massNum] = 1./5.*sqrt(3./2.)*(P1m*c_P1m).real-1./5.*sqrt(3./2.)*(P1p*c_P1p).real+1./7.*sqrt(3./2.)*(D1m*c_D1m).real-1./7.*sqrt(3./2.)*(D1p*c_D1p).real
        	allMoms_Mass[6][massNum] = 6./7.*sqrt(3./5.)*(P0*c_D0).real-6./7./sqrt(5.)*(P1m*c_D1m).real-6./7./sqrt(5.)*(P1p*c_D1p).real
        	allMoms_Mass[7][massNum] = 4./7.*sqrt(3./5.)*(P0*c_D1m).real+6./7./sqrt(5.)*(P1m*c_D0).real
        	allMoms_Mass[8][massNum] = 2./7.*sqrt(3./2.)*(P1m*c_D1m).real-2./7.*sqrt(3./2.)*(P1p*c_D1p).real
		# even though a term is like |D0|**2 we still have to take the real part since python will interpret it as a complex with 0j 
        	allMoms_Mass[9][massNum] = 2./7.*(D0*c_D0).real-4./21.*(D1m*c_D1m).real-4./21.*(D1p*c_D1p).real
        	allMoms_Mass[10][massNum] = 2./7.*sqrt(5./3.)*(D0*c_D1m).real
        	#allMoms_Mass[11][massNum] = sqrt(10./21.)*(D1m*c_D1m).real-sqrt(10./21.)*(D1p*c_D1p).real
        	allMoms_Mass[11][massNum] = 1./6.44*(D1m*c_D1m).real-1/6.44*(D1p*c_D1p).real


ampCanvas.Clear()
moments	= []
proj_moments = []

# make sure momNames and momList follow same ordering
momentList = ["hMoment00", "hMoment10","hMoment11","hMoment20","hMoment21","hMoment22","hMoment30","hMoment31","hMoment32","hMoment40","hMoment41", "hMoment42"]
dHist_moment = TH1D("anyMoment","anyMoment",60,0,60)
momentFile = TFile.Open("moments.root")
# we want to still wrtie stuff to amplitudesFile and just keep momentFile to read. We have to change the working directory back
amplitudesFile.cd()	

for momNum,momName in enumerate(momNames):
	print(momName)
	# we will also make a TH1D to hold the moments data so we can use them to calculate beam asymmetries
	moments.append( array( 'd' ) )
	proj_moments.append( array( 'd' ))
	momentFile.GetObject(momentList[momNum],dHist_moment)
	maxY=0
	minY=0
	for massNum in range(nBins):
		moments[momNum].append(allMoms_Mass[momNum][massNum])
		if (maxY<allMoms_Mass[momNum][massNum]):
			maxY = allMoms_Mass[momNum][massNum]
		if (minY>allMoms_Mass[momNum][massNum]):
			minY = allMoms_Mass[momNum][massNum]
		# recall that the 0th bin from root is the underflow!
        	proj_moments[momNum].append(dHist_moment.GetBinContent(massNum+1))
		print(allMoms_Mass[momNum][massNum])
		print(dHist_moment.GetBinContent(massNum+1))
		print("----------------")
	gr = TGraph ( nBins, masses, proj_moments[momNum] )
	gr.SetName("project_moments")
	gr.SetMarkerColor( 2 )
	gr.SetMarkerStyle( 21 )
	gr.SetTitle( momName )
	gr.GetXaxis().SetTitle( 'Masses (GeV)' )
	gr.GetYaxis().SetTitle( momName+' Moment amplitude')
	#gr.GetYaxis().SetRangeUser(minY*1.25,maxY*1.25)
	gr.Draw( 'AP' )
	gr2 = TGraph ( nBins, masses, moments[momNum] )
	gr2.SetName("moments")
	gr2.SetMarkerColor( 9 )
	gr2.SetMarkerStyle( 21 )
	gr2.Draw( 'P' )
	legend = TLegend(0.7,0.7,0.9,0.9)
	legend.AddEntry(gr, "Moments from project_moments")
	legend.AddEntry(gr2, "Moments using Shuang's ppt")
	legend.Draw()
	ampCanvas.SaveAs(saveFolder+"/"+momName+".pdf")
	gr2.Write(momName+"_shuang")
	#gr2.Write((momName_"vincent".c_str())
