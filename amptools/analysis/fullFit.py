import os
import glob
import random
import sys
import subprocess
import time
import fileinput
from loadCfg import loadCfg
from fitHelper import updateCfg, rndSampInits

cfg=loadCfg()
realSpace_low=cfg["REAL_LOW"]
realSpace_up=cfg["REAL_UPP"]
imSpace_low=cfg["IMAG_LOW"]
imSpace_up=cfg["IMAG_UPP"]
fitName=cfg["FIT_DIR"]
seedFile=cfg["SEED_FILE"]
maxIter=cfg["MAX_ITER"]
ampNames=cfg["AMP_NAMES"]
seedAmpInit=cfg["SEED_AMP_INIT"]
nBins=cfg["NUMBER_BINS"]
lowMass=cfg["LOW_MASS"]
uppMass=cfg["UPP_MASS"]
percentDeviation=cfg["PERCENT_DEVIATION"]
baseDeviation1=cfg["BASE_DEVIATION1"]
baseDeviation2=cfg["BASE_DEVIATION2"]
rampIter=cfg["RAMP_DEVIATION_ITER"]
rndSamp_flag=cfg["RND_SAMP"]

print("realSpace_low: "+str(realSpace_low))
print("realSpace_up: "+str(realSpace_up))
print("imSpace_low: "+str(imSpace_low))
print("imSpace_up: "+str(imSpace_up))
print("fitName: "+str(fitName))
print("seedFile: "+str(seedFile))
print("maxIter: "+str(maxIter))
print("ampNames: "+str(ampNames))
print("rndSamp_flag: "+str(rndSamp_flag))

start = time.time()

workingDir = os.getcwd()
print("\n")
print("current working directory: %s" % (workingDir))
fitDir = workingDir+"/"+fitName
print("fit directory: %s" % (fitDir))
random.seed(seedAmpInit)

# WHAT SCHEME SHOULD WE CHOOSE TO INITIALIZE AMPLITUDES? 
# --------------------------------------------------------------------------------------------
# SCHEME 1: Initial bin_0. Fit. For bin i take converged amplitudes from bin i-1 for all i > 0
# --------------------------------------------------------------------------------------------
# binsArray=range(nBins)
# prevBinsArray=range(-1,nBins)
# for ix, ibin in enumerate(binsArray):
#     if prevBinsArray[ix]>0:
#         prevBinDir=fitDir+"/bin_"+str(prevBinsArray[ix]-1)
#         updateCfg(prevBinDir+"/"+seedFile,"bin_"+str(binNum)+"-full.cfg")
# --------------------------------------------------------------------------------------------
# SCHEME 2: Use a0/a2 as our anchors. Initialize and fit those bins at a0/a2 peak. Do similar thing as SCHEME 1 
# --------------------------------------------------------------------------------------------
a0mass=0.98;
binWidth=(uppMass-lowMass)/nBins
a0Bin=int((a0mass-lowMass)/binWidth)
binsArrayR=range(a0Bin+1,nBins)
binsArrayL=range(a0Bin-1,-1,-1) 
binsArray=binsArrayR+binsArrayL
prevBinsArrayR=[i-1 for i in binsArrayR]
prevBinsArrayL=[i+1 for i in binsArrayL]
prevBinsArray=prevBinsArrayR+prevBinsArrayL;
curToPrev={binsArray[i] : prevBinsArray[i] for i in range(len(binsArray))}


for binNum in [a0Bin]+binsArray:
    os.chdir(fitDir)
    os.chdir("bin_"+str(binNum))
    print("\n================   bin_%i    =================" % (binNum))
    print("current working directory: %s" % (os.getcwd()))
    rndSamp=rndSamp_flag
    binCfg = "bin_"+str(binNum)+"-full.cfg"

    if binNum in curToPrev[binNum].keys():
        prevBinDir=fitDir+"/bin_"+str(curToPrev[binNum])
        updateCfg(prevBinDir+"/"+seedFile,"bin_"+str(binNum)+"-full.cfg")
    for j in range(maxIter):
        if (j<rampIter):
            baseDeviation=baseDeviation1
        else:
            baseDeviation=baseDeviation2
        if rndSamp:
            rndSampInits("bin_"+str(binNum)+"-full.cfg",percentDeviation,baseDeviation)
        print("First 5 current initial values for amplitudes:")
        print("----------------------------------")
        print("\n".join(subprocess.check_output(['grep','^initialize',binCfg]).split("\n")[:5]))
        print("----------------------------------")
        rndSamp=True # next iteration will start randomly sampling 

        print("----------------------------------------------------------")
        callFit = "fit -c "+binCfg+" -s "+seedFile
        print("Trying with seed="+str(seedAmpInit)+" at iteration="+str(j))
        print(callFit)
        # creating a file so that I could dump the results of the fit into it
        logFile = open("bin_"+str(binNum)+"-full.log","w+")
        subprocess.Popen(callFit.split(" "), stdout=logFile).wait()
        num_lines=0
        with open("param_init.cfg","r") as param_init_cfg:
            param_lines = param_init_cfg.readlines()
            for param_line in param_lines:
                if param_line[:10]=="initialize":
                    num_lines+=1
        print("num_lines: "+str(num_lines)+" should equal the total number of waves!")
        if (num_lines>0):
            print("Bin_"+str(binNum)+" converged with this set of initialization! -- Iteration:"+str(j))
            ###############################################################################################################
            # VERY IMPORT STEP TO GET SAVE THE .FIT FILE SO THAT PROEJCT MOMENTS CAN USE THE ONE THAT IS NOT BOOTSTRAPPPED
            os.rename("bin_"+str(binNum)+".fit","bin_"+str(binNum)+"-full.fit")
            niFiles=glob.glob("*ni")
            for niFile in niFiles:
                fileName=niFile.split(".")[0]
                fileType=niFile.split(".")[1]
                with open(niFile,"r") as ni:
                    lines=ni.readlines()
                    for line in lines:
                        if "nan" in line:
                            # Normalization integrals is related to the acceptance
                            raise ValueError(niFile+" contains nan values! gen matrix is on top and acc is on bottom")
                os.rename(niFile,fileName+"-full.ni")
            ##############################################################################################################
            break
        if (num_lines==0):
            print("Bin_"+str(binNum)+" did not converge with this set of initialization! Will try to resample if we can...")
        if (num_lines==0 and j==maxIter-1):
            print("Bin_"+str(binNum)+" did not converge maxIter times!")
            raise ValueError("Bin_"+str(binNum)+" did not converge maxIter times! Have to rethink this!")

stop = time.time()
print("Execution time in seconds: %s" % (stop-start))
