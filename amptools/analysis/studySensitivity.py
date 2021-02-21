import os
import math
import glob
import random
import sys
import subprocess
import time
import fileinput
from loadCfg import loadCfg
from fitHelper import updateCfg, rndSampInits
from multiprocessing import Pool

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

def sampleUniform(source,seed):
    '''
    This function takes in a config file and resample uniformly samples in the ranges specified by min/max array
    '''
    random.seed(seed)
    print("Randomly sampling initialization paramters with seed: "+str(seed))
    # Write the output 
    sampMin=-1000
    sampMax=1000
    with open(source,"r") as src:
        lines=src.readlines();
    with open(source,"w") as src:
        for line in lines:
            if line[:10]=="initialize":
                fields=line.split(" ")
                fields[3]=str(random.uniform(sampMin,sampMax))
                if fields[4]!="0.0":
                    fields[4]=str(random.uniform(sampMin,sampMax))
                src.write(" ".join(fields).rstrip())
                src.write("\n")
            else:
                src.write(line.rstrip())
                src.write("\n")

def getAmplitudesInBin(binNum):
    os.chdir("bin_"+str(binNum))
    print("\n================   bin_%i    =================" % (binNum))
    print("current working directory: %s" % (os.getcwd()))
    rndSamp=rndSamp_flag
    binCfg = "bin_"+str(binNum)+"-full.cfg"

    with open("amplitudes.txt","w") as outFile:
        haveHeader=False
        for j in range(numIters):
            sampleUniform("bin_"+str(binNum)+"-full.cfg",seed=seeds[j])
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
            #logFile = open("bin_"+str(binNum)+"-full.log","w+")
            subprocess.Popen(callFit.split(" ")).wait()
            num_lines=0
            with open("param_init.cfg","r") as param_init_cfg:
                param_lines = param_init_cfg.readlines()
                for param_line in param_lines:
                    if param_line[:10]=="initialize":
                        num_lines+=1
            print("num_lines: "+str(num_lines)+" should equal the total number of waves!")
            if (num_lines>0):
                print("Bin_"+str(binNum)+" converged with this set of initialization! -- Iteration:"+str(j))
                resultsFile="bin_"+str(binNum)+".fit"
                print("\n\nGetting amplitude intensities")
                getAmplitudeCmd="studySensitivity "+resultsFile+" "+str(j)
                output=os.popen(getAmplitudeCmd).read()
                output=output.split("\n")
                for out in output:
                    if len(out)>0:
                        # We only want to print the header once... so lets just search for it looking for lines that have \t in the line and does not start with a number...
                        # this probably wont work for all cases
                        if len(out.split("\t"))>1:
                            if out[0].isdigit():
                                outFile.write(out+"\n")
                            if not haveHeader:
                                outFile.write(out+"\n")
                                haveHeader=True
            if (num_lines==0):
                print("Bin_"+str(binNum)+" did not converge with this set of initialization!")
    
    os.chdir("..")

def getAmplitudesInBins(binArray):
    for binNum in binArray:
        getAmplitudesInBin(binNum)

### CHOOSE BIN NUMBER
if __name__ == '__main__':
    os.chdir(fitDir)
    numIters=500
    seeds=[random.randint(1,100000) for _ in range(numIters)]
    chunkSize=2
    processes=33
    binArrays=[]
    for iChunk in range(processes):
        binArray=range(iChunk*chunkSize,(iChunk+1)*chunkSize)
        binArray=[iBin for iBin in binArray if iBin<nBins]
        binArrays.append(binArray)
    print("BINS FOR THE PROCESSES TO RUN OVER")
    print(binArrays)
    p=Pool(processes)
    p.map(getAmplitudesInBins, binArrays)
    p.terminate()


stop = time.time()
print("Execution time in seconds: %s" % (stop-start))
