from loadCfg import loadCfg
import glob 
import math 
from math import sqrt
from array import array
import csv
from ROOT import TCanvas, TGraphErrors, TGraph, TH1D, TLegend
from ROOT import gROOT, TFile
import os
import sys
import shutil
from calcMoments import calcMoments as calcMoments

cfg=loadCfg()
rndSamp_flag=cfg["RND_SAMP"]
ampNames = cfg["AMP_NAMES"]
nBins=cfg["NUMBER_BINS"]
lowMass=cfg["LOW_MASS"]
upMass=cfg["UPP_MASS"]
fitDir=cfg["FIT_DIR"]
saveAmps=cfg["SAVE_AMPS"]
momNames=cfg["MOM_NAMES"]
projMoms=True
if projMoms:
	momentList=cfg["PROJ_MOM_NAMES"]

saveFolder="amplitudes"

makeMoments=True
# ---- Should we overlay bootstrapped errors of the first N seeds to see how inclusion of more seeds changes things ----
overlaySomeSeeds=True
firstNSeeds=50
# ---- exports the intensity and std and std_firstN for each amplitude into the csv below
if os.path.isdir("csvData"):
	shutil.rmtree("csvData")
	os.mkdir("csvData")
else:
	os.mkdir("csvData")
if rndSamp_flag:
	csvName_converged="csvData/csvData-rnd.txt"
else:
	csvName_converged="csvData/csvData-converged.txt"
# ---- 
stdVsNSeeds=True
binsToSave=[10,30,45]
csvNames_binVsStd=[]
csvNames_binVsSeededAmps=[]
# ----- We save the amplitudes 
for massIdx, massNum in enumerate(binsToSave):
	csvNames_binVsStd.append("csvData/csvData-binVsStd_bin"+str(massNum)+".txt")
	csvNames_binVsSeededAmps.append("csvData/csvData-binVsSeededAmps_bin"+str(massNum)+".txt")
	



print(os.path.isdir(saveFolder))
if os.path.isdir(saveFolder):
	shutil.rmtree(saveFolder)
	print("removing %s folder!" % saveFolder)
	# os module only removes empty directories
	print("making %s folder!" % saveFolder)
	os.mkdir(saveFolder)
	os.mkdir(saveFolder+"/moments")
else:
	print("Error: %s folder not found!" % saveFolder)
	print("making %s folder!" % saveFolder)
	os.mkdir(saveFolder)
	os.mkdir(saveFolder+"/moments")


#momNames_ph=cfg["MOM_NAMES_PH"]
print("Using amplitudes: ")
print(ampNames)
print("Using moments:")
print(momNames)
#print(momNames_ph)

ampCanvas = TCanvas("ampCanvas", "amplitudes", 1200,600)

# Not sure how to draw root objects with +/- prefix .... so we can just conver them to p/m
ampNamesOnlyLetters = []
for ampName in ampNames:
	ampNamesOnlyLetters.append(ampName.replace("+","p").replace("-","m"))

verbose=True

# need to use this function to sort the list of ampFiles
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
numMoms = int(len(momNames))
allAmps_MassSeed = [[[0 for col in range(nSeeds)] for row in range(nBins)] for depth in range(numAmps)]

numEntries = 2*numAmps+1
print("Number of entries in each line of converged amps file: "+str(numEntries))
print("Corresponds to 1 mass value and 2 components for each ampltiude (real/imag)")

allAmps_Mass = [[0 for col in range(nBins)] for row in range(numAmps)]

params = [nBins, numMoms, numAmps]

print("------------------- FINISHED WITH PRELIMINARIES\n\n")
# complexAmps = complex amplitudes matrix
complexAmps = [[0 for col in range(nBins)] for row in range(numAmps)]
with open(fitDir+"/"+saveAmps+".txt", "r") as ampFile:
	print("Opened "+fitDir+"/"+saveAmps+".txt")
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
print("-------------- GOT THE FULL DATA COMPLEX AMPS AND AMP SQ\n\n")

# This code is basically a loop of the above to grab all the boostrapped datasets, so all debugging should be done with that one. 
print("bootstrapped amplitudes squared values are being calculated....")
seedToIdx = {}
complexAmps_Seed = [[0 for col in range(nBins)] for row in range(numAmps)]
complexAmps_Seeds = [] 
for idxSeed,ampFile in enumerate(allAmpFiles):
	seedNum = int(ampFile.split("seed")[1].split(".")[0])
	seedToIdx[seedNum] = idxSeed	
	assert seedNum==idxSeed, "The values should match! seedNum, idxSeed: {0},{1}".format(seedNum,idxSeed)
	sys.stdout.write("#")
	#print(ampFile+" is being loaded...")
	with open(ampFile, "r") as ampFile:
		for massNum,massLine in enumerate(ampFile):
	        	numSplit = len(massLine.split(' '))-1
			if (numSplit == numEntries-1):
				values = [float(item.strip()) for item in massLine.split(' ', numSplit)]
				for ampNum in range(numAmps):
					amps = values[2*ampNum+1]**2 + values[2*ampNum+2]**2
					complexAmps_Seed[ampNum][massNum] = values[2*ampNum+1]+values[2*ampNum+2]*1j
					#print("ampNum, massNum, seedNum: {0}, {1}, {2}".format(ampNum,massNum,idxSeed))
					allAmps_MassSeed[ampNum][massNum][idxSeed] = amps
			else:
				raise ValueError("NOT THE CORRECT NUMBER OF ENTRIES!")
	complexAmps_Seeds.append(complexAmps_Seed)
sys.stdout.write("\n")
sys.stdout.flush()
print("-------------- GOT ALL BOOTSTRAPPED COMPLEX AMPS AND AMP SQ\n\n")



# this output helps in determining how we should sum over the seeds to get the std and mean
print("SHAPE OF allAmps_MassSeed just to check the ordering of the data")
print("Shape: {0}, nBins: {1}, nSeeds: {2}, numMoms: {3}".format([len(allAmps_MassSeed),len(allAmps_MassSeed[0]),len(allAmps_MassSeed[0][0])],nBins,nSeeds,numAmps))

print("calculating the mean and population std (instead of the sample std)  of the amplitude squared values over all seeds")
allAmps_MassSeed_mean_firstN=[[0 for col in range(nBins)] for row in range(numAmps)]
allAmps_MassSeed_std_firstN=[[0 for col in range(nBins)] for row in range(numAmps)]
allAmps_MassSeed_mean=[[0 for col in range(nBins)] for row in range(numAmps)]
allAmps_MassSeed_std=[[0 for col in range(nBins)] for row in range(numAmps)]
for ampNum in range(numAmps):
	for massNum in range(nBins):
		for seedNum in range(nSeeds):
			allAmps_MassSeed_mean[ampNum][massNum]+=allAmps_MassSeed[ampNum][massNum][seedNum]
			if seedNum<firstNSeeds:
				allAmps_MassSeed_mean_firstN[ampNum][massNum]+=allAmps_MassSeed[ampNum][massNum][seedNum]
		allAmps_MassSeed_mean[ampNum][massNum]/=nSeeds
		allAmps_MassSeed_mean_firstN[ampNum][massNum]/=firstNSeeds
		for seedNum in range(nSeeds):
			allAmps_MassSeed_std[ampNum][massNum]+=(allAmps_MassSeed[ampNum][massNum][seedNum]-allAmps_MassSeed_mean[ampNum][massNum])**2
			if seedNum<firstNSeeds:
				allAmps_MassSeed_std_firstN[ampNum][massNum]+=(allAmps_MassSeed[ampNum][massNum][seedNum]-allAmps_MassSeed_mean_firstN[ampNum][massNum])**2
		allAmps_MassSeed_std[ampNum][massNum]/=(nSeeds-1)
		allAmps_MassSeed_std[ampNum][massNum]=float(math.sqrt(allAmps_MassSeed_std[ampNum][massNum]))
		allAmps_MassSeed_std_firstN[ampNum][massNum]/=(firstNSeeds-1)
		allAmps_MassSeed_std_firstN[ampNum][massNum]=float(math.sqrt(allAmps_MassSeed_std_firstN[ampNum][massNum]))



if stdVsNSeeds:
	# This one will be to store successive stds as we include more seeds
	allAmps_MassSeed_std_binVsStd=[[[0 for col in range(nSeeds)] for row in range(numAmps)] for width in range(len(binsToSave))]
	# This one will store the raw amplitude values for each seed. We will save this to a file which we can do a personal calculation on to double check
	allAmps_MassSeed_std_binVsSeededAmps=[[[0 for col in range(nSeeds)] for row in range(numAmps)] for width in range(len(binsToSave))]
	for ampNum in range(numAmps):
		for massIdx,massNum in enumerate(binsToSave):
			# we will save all the std values as we include more seeds to the calculation. 
			# recall this maxSeed is the maximum seed number to include in a calcuation for a standard devaition. We will constantly recalculate std as we include more seeds
			checkMaxSeed=6 # we can use this to manually output the values to see if they compare to the data in the outputted text file
			for maxSeed in range(nSeeds):
				currentAvg=0
				for seedNum in range(maxSeed+1): #since maxSeed starts at 0 even though the length is 1
					currentAvg+=allAmps_MassSeed[ampNum][massNum][seedNum]
				currentAvg/=(maxSeed+1) 
				if maxSeed==checkMaxSeed and ampNum==0 and massIdx==0:
					print("currentAvg: {0}".format(currentAvg))
				for seedNum in range(maxSeed+1):								
					# at the first iteration with seedNum the value at allAmps_MassSeed_std_binVsStd[massIdx][ampNum][maxSeed] should be 0 since we initialized it like that
					# subsequent iterations would have a non zero value since we added to it BUT seedNum will no longer be 0. So hopefully this condition will just make sure
					# things are initialized properly
					if seedNum==0 and allAmps_MassSeed_std_binVsStd[massIdx][ampNum][maxSeed] != 0:
						print("The shape of allAmps_MassSeed_std_binVsStd is {0},{1},{2}".format(nSeeds,numAmps,len(binsToSave)))
						print("The error fails at massIdx,ampNum,maxSeed: {0},{1},{2}".format(massIdx,ampNum,maxSeed))
						raise ValueError("allAmps_MassSeed_std_binVsStd not set properly!")
					else:
						allAmps_MassSeed_std_binVsStd[massIdx][ampNum][maxSeed]+=(allAmps_MassSeed[ampNum][massNum][seedNum]-currentAvg)**2
				allAmps_MassSeed_std_binVsStd[massIdx][ampNum][maxSeed]/=(maxSeed+1)
				allAmps_MassSeed_std_binVsStd[massIdx][ampNum][maxSeed]=float(math.sqrt(allAmps_MassSeed_std_binVsStd[massIdx][ampNum][maxSeed]))
				if maxSeed==checkMaxSeed and ampNum==0 and massIdx==0:
					print("currentStd: {0}".format(allAmps_MassSeed_std_binVsStd[massIdx][ampNum][maxSeed]))
				# This should essentially just be a rotation of the data, the information should still be the same. Just swapping two indicies.
				allAmps_MassSeed_std_binVsSeededAmps[massIdx][ampNum][maxSeed] = allAmps_MassSeed[ampNum][massNum][seedNum]
			

	for massIdx,csvName in enumerate(csvNames_binVsStd):
		with open(csvName,"w+") as csv_file:
			csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
			for ampNum in range(numAmps):
				csv_writer.writerow(allAmps_MassSeed_std_binVsStd[massIdx][ampNum])
	for massIdx,csvName in enumerate(csvNames_binVsSeededAmps):
		with open(csvName,"w+") as csv_file:
			csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
			for ampNum in range(numAmps):
				csv_writer.writerow(allAmps_MassSeed_std_binVsSeededAmps[massIdx][ampNum])
		
		
with open(csvName_converged,"w+") as csv_file:
	csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
	# this is kind of the weird step since we do all this importing into python and calculating 
	# then we have to rebulild the arrays to be usuable with pyrootampFile=ROOT.TFile.Open("amplitudes.root")
	print("Rebuilding the python arrays to use the array module which can then be used by pyroot")
	amplitudesFile=TFile.Open(saveFolder+"/amplitudes.root","RECREATE")
	for ampNum in range(numAmps):
		amplitude = array( 'd' )
		masses = array( 'd' )
		amplitude_std = array( 'd' )
		mass_std = array( 'd' )
		amplitude_std_firstN = array( 'd' )
		mass_std_firstN = array( 'd' )
		for massNum in range(nBins):
			amplitude.append( allAmps_Mass[ampNum][massNum] )
			amplitude_std.append( allAmps_MassSeed_std[ampNum][massNum] )
			amplitude_std_firstN.append( allAmps_MassSeed_std_firstN[ampNum][massNum] )
			# masses start at lowMass up to upMass with step size massStep. We shift up to lowMass+massStep first which is the upper edge of the first bin. 
			masses.append( massStep*massNum+lowMass+massStep )	
			mass_std.append( 0.0 )
		gr = TGraphErrors( nBins, masses, amplitude, mass_std, amplitude_std )
		gr.GetYaxis().SetRangeUser(0,100000)
		gr.SetMarkerColorAlpha( 9, 0.75 )
		gr.SetMarkerStyle( 21 )
		gr.SetTitle( ampNames[ampNum] )
		gr.GetXaxis().SetTitle( 'Masses (GeV)' )
		gr.GetYaxis().SetTitle( 'Amplitude' )
		gr.Draw( 'AP' )
		if overlaySomeSeeds:
			gr_firstN = TGraphErrors( nBins, masses, amplitude, mass_std, amplitude_std_firstN )
			gr_firstN.SetMarkerColorAlpha( 4, 0)
			gr_firstN.SetLineColor( 2 )
			gr_firstN.Draw( 'P' )
		ampCanvas.SaveAs(saveFolder+"/"+ampNamesOnlyLetters[ampNum]+".png")
		gr.Write(ampNamesOnlyLetters[ampNum])
		csv_writer.writerow(amplitude)
		csv_writer.writerow(amplitude_std)
		csv_writer.writerow(amplitude_std_firstN)
	print("------------- FINISHED MAKING AMPLITUDE GRAPHS\n\n")

		
	# sum the errors in quadruture
	print("Calculating the summed amplitudes with errors combined in quadruture")
	ampTot = array ( 'd' )
	ampTot_std = array ( 'd' ) 
	ampTot_std_firstN = array ( 'd' ) 
	for massNum in range(nBins):
		ampTot.append(0) # initial value we will be adding to
		ampTot_std.append(0)
		ampTot_std_firstN.append(0)
		for ampNum in range(numAmps):
			ampTot[massNum] += allAmps_Mass[ampNum][massNum]
			ampTot_std[massNum] += allAmps_MassSeed_std[ampNum][massNum]**2
			ampTot_std_firstN[massNum] += allAmps_MassSeed_std_firstN[ampNum][massNum]**2
		ampTot_std[massNum] = sqrt(ampTot_std[massNum])
		ampTot_std_firstN[massNum] = sqrt(ampTot_std_firstN[massNum])
	gr = TGraphErrors( nBins, masses, ampTot, mass_std, ampTot_std )
	gr.SetMarkerColorAlpha( 9, 0.75 )
	gr.SetMarkerStyle( 21 )
	gr.SetTitle( "Summed Amps" )
	gr.GetXaxis().SetTitle( 'Masses (GeV)' )
	gr.GetYaxis().SetTitle( 'Amplitudes Summed' )
	gr.Draw( 'AP' )
	if overlaySomeSeeds:
		gr_firstN = TGraphErrors( nBins, masses, ampTot, mass_std, ampTot_std_firstN )
		gr_firstN.SetMarkerColorAlpha( 4, 0)
		gr_firstN.SetLineColor( 2 )
		gr_firstN.Draw( 'P' )
	ampCanvas.SaveAs(saveFolder+"/summedAmps.png")
	gr.Write("summedAmps")
	csv_writer.writerow(ampTot)
	csv_writer.writerow(ampTot_std)
	csv_writer.writerow(ampTot_std_firstN)

print("Calculating the moments from full data and avg,std from the bootstrapped data")
allMoms_MassSeed=[]
allMoms_MassSeed_ph=[]
if makeMoments:
	vS0m,vP0m,vP1m,vD0m,vD1m,vP1p,vD1p=complexAmps
	#vS0,vP0,vP1,vD0,vD1,vD2=complexAmps
	allMoms_Mass = calcMoments(params,vS0m,vP0m,vP1m,vD0m,vD1m,vP1p,vD1p)
	#allMoms_Mass = calcMoments(params,vS0,vP0,vP1,vD0,vD1,vD2)
	for i in range(nSeeds):
		vS0m,vP0m,vP1m,vD0m,vD1m,vP1p,vD1p=complexAmps_Seeds[i]
		#vS0,vP0,vP1,vD0,vD1,vD2=complexAmps_Seeds[i]
		allMoms_MassSeed.append(calcMoments(params,vS0m,vP0m,vP1m,vD0m,vD1m,vP1p,vD1p))
		#allMoms_MassSeed.append(calcMoments(params,vS0,vP0,vP1,vD0,vD1,vD2))
	# this output helps in determining how we should sum over the seeds to get the std and mean
	print("SHAPE OF allMoms_MassSeed just to check the ordering of the data")
	print("Shape: {0}, nBins: {1}, nSeeds: {2}, numMoms: {3}".format([len(allMoms_MassSeed),len(allMoms_MassSeed[0]),len(allMoms_MassSeed[0][0])],nBins,nSeeds,numMoms))
	
	print("calculating the mean and population std (instead of the sample std)  of the moments over all seeds")
	allMoms_MassSeed_mean=[[0 for col in range(nBins)] for row in range(numMoms)]
	allMoms_MassSeed_std=[[0 for col in range(nBins)] for row in range(numMoms)]
	for momNum in range(numMoms):
		for massNum in range(nBins):
			for seedNum in range(nSeeds):
				allMoms_MassSeed_mean[momNum][massNum]+=allMoms_MassSeed[seedNum][momNum][massNum]
			allMoms_MassSeed_mean[momNum][massNum]/=nSeeds
			for seedNum in range(nSeeds):
				allMoms_MassSeed_std[momNum][massNum]+=(allMoms_MassSeed[seedNum][momNum][massNum]-allMoms_MassSeed_mean[momNum][massNum])**2
			allMoms_MassSeed_std[momNum][massNum]/=(nSeeds)
			allMoms_MassSeed_std[momNum][massNum]=float(math.sqrt(allMoms_MassSeed_std[momNum][massNum]))

	ampCanvas.Clear()
	moments	= []
	moments_std = []
	
	# make sure momNames and momList follow same ordering to follow same convention as project_moments
	if projMoms:
		dHist_moment = TH1D("anyMoment","anyMoment",60,0,60)
		proj_moments = []
		momentFile = TFile.Open("moments-uniform.root")
	# we want to still wrtie stuff to amplitudesFile and just keep momentFile to read. We have to change the working directory back
	amplitudesFile.cd()	
	
	for momNum,momName in enumerate(momNames):
		# we will also make a TH1D to hold the moments data so we can use them to calculate beam asymmetries
		moments.append( array( 'd' ) )
		moments_std.append( array( 'd' ))
		if projMoms:
			proj_moments.append( array( 'd' ))
			momentFile.GetObject(momentList[momNum],dHist_moment)
		maxY=0
		minY=0
		for massNum in range(nBins):
			moments[momNum].append(allMoms_Mass[momNum][massNum])
			moments_std[momNum].append( allMoms_MassSeed_std[momNum][massNum] )
			# checking min and max values so we can set the axes manually. Sometimes needed
			if (maxY<allMoms_Mass[momNum][massNum]):
				maxY = allMoms_Mass[momNum][massNum]
			if (minY>allMoms_Mass[momNum][massNum]):
				minY = allMoms_Mass[momNum][massNum]
			# recall that the 0th bin from root is the underflow!
			if projMoms:
	        		proj_moments[momNum].append(dHist_moment.GetBinContent(massNum+1))
			if verbose:
				print(momName+" moment")
				print(allMoms_Mass[momNum][massNum])
				if projMoms:
					print(dHist_moment.GetBinContent(massNum+1))
				print("----------------")
		if projMoms:
			gr = TGraph ( nBins, masses, proj_moments[momNum] )
			gr.SetName("project_moments")
			gr.SetMarkerColor( 2 )
			gr.SetMarkerStyle( 21 )
			gr.SetTitle( momName )
			gr.GetXaxis().SetTitle( 'Masses (GeV)' )
			gr.GetYaxis().SetTitle( momName+' Moment amplitude')
			gr.GetYaxis().SetRangeUser(minY*1.25,maxY*1.25)
			gr.Draw( 'AP' )
			gr2_drawOption = "P"
		else:
			gr2_drawOption = "AP"
		gr2 = TGraphErrors ( nBins, masses, moments[momNum], mass_std, moments_std[momNum] )
		gr2.SetName(momName)
		gr2.SetTitle(momName)
		gr2.SetMarkerColor( 9 )
		gr2.SetMarkerStyle( 21 )
		gr2.Draw( gr2_drawOption )
		legend = TLegend(0.7,0.7,0.9,0.9)
		if projMoms:
			legend.AddEntry(gr, "Moments from project_moments")
		legend.AddEntry(gr2, "Shuang's Moments")
		legend.Draw()
		ampCanvas.SaveAs(saveFolder+"/moments/"+momName+".png")
		gr2.Write(momName+"_comparison")
	print("------------- FINISHED MAKING MOMENTS GRAPHS DERIVED FROM PION BEAM\n\n")


