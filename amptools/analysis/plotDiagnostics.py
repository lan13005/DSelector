import os
import subprocess
import ROOT

fitDir="EtaPi_fit"
bins=65 #in case we dont want to run over all the bins yet
fullNumBins=65
minMass=0.7
maxMass=2.0
binWidth=(maxMass-minMass)/fullNumBins


baseDir=os.getcwd()
for iBin in range(bins):
	binDir=baseDir+"/"+fitDir+"/bin_"+str(iBin)
	fitFile="bin_"+str(iBin)+".fit"
	print "moving into: "+binDir
	os.chdir(binDir)
	print "  using twopi_plotter_amp on "+fitFile
    	subprocess.Popen("twopi_plotter_amp "+fitFile+" -o diagnostic.root",shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).wait()

























	
	
