import os
import sys
import subprocess
import time
import fileinput
import random
import shutil
import glob
from loadCfg import loadCfg
from fitHelper import getAmps


## max iteration of attempts to fit a specific bootstrapped dataset with each iteration changing in initialization
#maxIter = 20
## number of bootstrap data to build and fit
#nSeeds=50
## how much to vary the inital converged amplitudes.
#percentDeviation=0.5
#baseDeviation=30
## pick some common seed that is shared among all the bins that bootFit_spawn.sh
#fitName = "divideRoot"
#seedFile = "param_init.cfg"
#ampNames=["S0+","P0+","P1+","D0+","D1+","P1-","D1-"]

cfg=loadCfg()
fitName=cfg["FIT_DIR"]
seedFile=cfg["SEED_FILE"]
maxIter=cfg["MAX_ITER"]
ampNames=cfg["AMP_NAMES"]
nSeeds=cfg["NUMBER_SEEDS"]
nBins=cfg["NUMBER_BINS"]
percentDeviation=cfg["PERCENT_DEVIATION"]
baseDeviation1=cfg["BASE_DEVIATION1"]
baseDeviation2=cfg["BASE_DEVIATION2"]
rampIter=cfg["RAMP_DEVIATION_ITER"]
mySeed=cfg["STARTING_SEED"]
seedShift=cfg["SEED_SHIFT"]
seedAmpInit=cfg["SEED_AMP_INIT"]
singleBin=cfg["SINGLE_BIN"]
# this flag will split the total processes among the various seeds in a single bin 
useMultipleProcsPerBin=cfg["USE_MULTIPLE_PROCS_PER_BIN"]
# this flag will determine if we randomly sample around the converged amplitudes or use the converged amplitudes as a starting point for the bootstrapped dataset	
rndSamp_flag=cfg["RND_SAMP"]

print("fitName: "+str(fitName))
print("nBins: "+str(nBins))
print("seedFile: "+str(seedFile))
print("maxIter: "+str(maxIter))
print("ampNames: "+str(ampNames))
print("nSeeds: "+str(nSeeds))
print("percentDeviation: "+str(percentDeviation))
print("baseDeviation1: "+str(baseDeviation1))
print("baseDeviation2: "+str(baseDeviation2))
print("rampIter: "+str(rampIter))
print("mySeed: "+str(mySeed))
print("seedShift: "+str(seedShift))

if len(sys.argv)!=4:
    raise ValueError("Must take iProcess, binsPerProc, processSpawned as an argument!")
else:
    iProcess = int(sys.argv[1])
    binsPerProc = int(sys.argv[2])
    processSpawned = int(sys.argv[3])

# pre grabbing all the seeds, for some reasion the seeds might get all mixed up if I do it later.
seeds=[]
random.seed(mySeed)
for seedNum in range(nSeeds):
    seeds.append(random.randint(0,1E9))
print("allSeeds: ")
print(seeds)
seedNums=range(nSeeds)
print(seedNums)
# if useMultipleProcsPerBin is set we will split up seeds array depending on the iProcess that is passed to this program in bootfit_spawn.sh. 
# the number of total spawned processes should be a factor of the total number of seeds, nSeeds. 
if useMultipleProcsPerBin:
    if nSeeds%processSpawned==0:
        print("These are the groupings of the seeds:")
        for i in range(processSpawned):
            if iProcess == i:
                print("--------------------------------- USING THESE SEEDS ---------------------------------")
            print("Group {2}: shown below which are from seeds[{0}:{1}]".format(i*(nSeeds/processSpawned),(i+1)*(nSeeds/processSpawned),i))
            print(seeds[i*(nSeeds/processSpawned):(i+1)*(nSeeds/processSpawned)])
            if iProcess == i:
                print("-------------------------------------------------------------------------------------")
        seeds=seeds[iProcess*(nSeeds/processSpawned):(iProcess+1)*(nSeeds/processSpawned)]
        seedNums=seedNums[iProcess*(nSeeds/processSpawned):(iProcess+1)*(nSeeds/processSpawned)]
    else:
    	raise ValueError("processSpawned is not a factor of nSeeds! This would result in an untested section of the code where a process might have to handle more seeds than the others!")
	

#getting the bins that we will run over for this process
#we can also select out a bin and just run multiple processes on it over all the seeds
def getBinsArray(i,bpt,threads,singleBin):
    if singleBin in range(nBins):
        print("singleBin,nBins = {0},{1}".format(singleBin,nBins))
        print("singleBin status: Only looking a single bin!")
        array1 = [singleBin] 
    else:
        print("singleBin,nBins = {0},{1}".format(singleBin,nBins))
        print("singleBin status: Looking at multiple bins!")
        array1 = [i]
        bptm1 = bpt-1
        for iproc in range(bptm1):
            nextBin=i+(iproc+1)*threads
            if nextBin<nBins: # make sure we dont go over the max number of bins.
                array1.append(nextBin)
    return array1
binsArray = getBinsArray(iProcess,binsPerProc,processSpawned,singleBin)
print("Running over these bins for this process:")
print(binsArray)

start = time.time()

workingDir = os.getcwd()
print("\n")
print("current working directory: %s" % (workingDir))

fitDir = workingDir+"/"+fitName


print("fit directory: %s" % (fitDir))

random.seed(seedAmpInit)
#def getAmps(binDir,percentDeviation,baseDeviation,rndSamp,verbose=True):
#    # this function works in the current directory... make sure to cd
#    if(rndSamp==False):
#        print("We will use the converged values instead of randomly sampling around them!")
#        percentDeviation=0
#        baseDeviation=0
#    else:
#        print("using percentDeviation as: "+str(percentDeviation))
#        print("using baseDeviation as: "+str(baseDeviation))
#    realAmps=[]
#    imAmps=[]
#    isreal=[]
#    print("Opening: "+binDir+"/"+seedFile)
#    with open(binDir+"/"+seedFile,"r") as param_init_cfg:
#        paramLines = param_init_cfg.readlines()
#        lineCounter=0
#        ampCounter=-1
#        for paramLine in paramLines:
#            if lineCounter%2==0:
#                ampCounter+=1
#            realAmp = float(paramLine.split(" ")[3])
#            imAmp = float(paramLine.split(" ")[4])
#            lenParams = len(paramLine.split(" "))
#            real_low = float(realAmp*(1-percentDeviation)-baseDeviation)
#            real_up = float(realAmp*(1+percentDeviation)+baseDeviation)
#            im_low = float(imAmp*(1-percentDeviation)-baseDeviation)
#            im_up = float(imAmp*(1+percentDeviation)+baseDeviation)
#            realSampled = str(random.uniform(real_low,real_up))
#            if (lenParams>5) and (paramLine.split(" ")[5].rstrip()=="real"):
#                imSampled="0.0"
#                isreal.append(True)
#            else: 
#                imSampled= str(random.uniform(im_low,im_up))
#                isreal.append(False)
#            if(verbose):
#                print("Below is paramLists, check if below values match:")
#                print(paramLine.split(" "))
#                print("Amps from "+seedFile+": "+ampNames[ampCounter])
#                print("realAmp: "+str(realAmp))
#                print("imAmp: "+str(imAmp))
#                print("Sampling from the below ranges")
#                print("real_low, real_up: "+str(real_low)+", "+str(real_up))
#                print("im_low, im_up: "+str(im_low)+", "+str(im_up))
#                print("realSampled, imSampled: {0},{1}".format(realSampled,imSampled))
#            realAmps.append(realSampled)
#            imAmps.append(imSampled)
#            lineCounter+=1
#            print("----------------------------")
#        if (len(realAmps)==len(imAmps)==len(isreal)):
#            print("length of realAmps, imAmps, isreal is good!")
#        else:
#            print("length of realAmps: "+str(len(realAmps)))
#            print("length of imAmps: "+str(len(imAmps)))
#            print("length of isreal: "+str(len(isreal)))
#            raise ValueError("length of realAmps, imAmps, isreal is not good!")
#        #print("==================================\n\n")
#    return realAmps, imAmps, isreal

def writeCfg(ampNames,bin_cfg,bin_cfg_seed,seed,binDir,percentDeviation,baseDeviation,rndSamp):
    # this function works in the current directory... make sure to cd
    ampCounter=0
    shutil.copyfile(bin_cfg, bin_cfg_seed)
    bin_cfg_base=bin_cfg.split(".")[0]
    bin_cfg_seed_base=bin_cfg_seed.split(".")[0]
    sedCall="sed -i s/"+bin_cfg_base+"/"+bin_cfg_seed_base+"/g "+bin_cfg_seed
    print(sedCall)
    subprocess.Popen(sedCall.split(" ")).wait()
    realAmps,imAmps,isreal = getAmps(ampNames,binDir,seedFile,rndSamp,percentDeviation,baseDeviation,True)
    for line in fileinput.input(bin_cfg_seed,inplace=1):
        if line.strip().startswith("data"):
            line = line.replace("ROOTDataReader","ROOTDataReaderBootstrap").replace(".root",".root "+str(seed))
        if line.strip().startswith("initialize"):
            posNegTag=line.split(" ")[1].split("::")[1]
            ampName=line.split(" ")[1].split("::")[2]
            if (isreal[ampCounter]):
                line = "initialize EtaPi::"+posNegTag+"::"+ampName+" cartesian "+realAmps[ampCounter]+" "+imAmps[ampCounter]+" real\n"
            else:
                line = "initialize EtaPi::"+posNegTag+"::"+ampName+" cartesian "+realAmps[ampCounter]+" "+imAmps[ampCounter]+"\n"
            ampCounter+=1
        sys.stdout.write(line)
    return True

for binNum in binsArray:
    binDir=fitDir+"/bin_"+str(binNum)
    if os.path.exists(binDir+"/seed-analysis") and os.path.isdir(binDir+"/seed-analysis"):
        shutil.rmtree(binDir+"/seed-analysis")
        os.mkdir(binDir+"/seed-analysis")
    else:
        os.mkdir(binDir+"/seed-analysis")

for seedIdx in range(len(seedNums)):
    seed = seeds[seedIdx]
    seedNum = seedNums[seedIdx]
    print("\n**********************************************************************\n***************************** NEXT SEED *********************************\n")
    for binNum in binsArray:
        print("------------------------ NEW BIN -------------------------\n-------------------------------------------------------")
        print("Running on seed,seedNum,binNum={0},{1},{2}".format(seed,seedNum,binNum))
        binDir=fitDir+"/bin_"+str(binNum)
        os.chdir(binDir)
        print("================   bin_%i    =================" % (binNum))
        print("current working directory: %s" % (os.getcwd()))
        seedShifted=seedNum+seedShift
        bin_cfg_seed_base="bin_"+str(binNum)+"-seed"+str(seedShifted)
        bin_cfg_seed = bin_cfg_seed_base+".cfg"
        bin_cfg_base = "bin_"+str(binNum)
        bin_cfg = bin_cfg_base+"-full.cfg"
        param_init_seed = "seed-analysis/seed"+str(seedShifted)+"-"+seedFile
        
        print("Copying "+bin_cfg_seed+" from "+bin_cfg+"\n--------------------------------------------------")
        # if rndSamp_flag=True we will start with our converged amplitudes and if that doesnt converge we will randomly sample a region around the converged values
        # if rndSamp_flag=False we will not start with the converged amplitudes
        rndSamp=rndSamp_flag
        for j in range(maxIter):
            print("======================== NEW ATTEMPT =========================\n=======================================================")
            if (j<rampIter):
                baseDeviation=baseDeviation1
            else:
                baseDeviation=baseDeviation2
            writeCfg(ampNames,bin_cfg,bin_cfg_seed,seed,binDir,percentDeviation,baseDeviation,rndSamp)
            rndSamp=True # next iteration will start randomly sampling 
            # creating a file so that I could dump the results of the fit into it
            callFit = "fit -c "+bin_cfg_seed+" -s "+param_init_seed
            print(callFit)
            #DEVNULL = open(os.devnull, 'wb')
            logFileName="seed-analysis/bin_"+str(binNum)+"-seed"+str(seedShifted)+".log"
            with open(logFileName,"w+") as logFile:
                subprocess.Popen(callFit.split(" "), stdout=logFile).wait()
            num_lines=0
            if os.path.isfile(param_init_seed):	
                with open(param_init_seed,"r") as param_init_cfg:
                    param_lines = param_init_cfg.readlines()
                    for param_line in param_lines:
                        num_lines+=1
            else:
                print(param_init_seed+" is somehow not made reason below:!")
                if "STATUS=CONVERGED" in open(logFileName).read():
                    raise ValueError("Somehow the log file converged but the param_init_seed file was not made! Fix me!")
                else:
                    print("The seed was not made because the fit did not converge!")
            print("THIS IS ITERATION {0}".format(j))
            print("num_lines: "+str(num_lines)+" should equal the total number of waves!")
            if (num_lines>0):
                print("Bin_"+str(binNum)+" converged with this set of initialization! -- Iteration:"+str(j))
                ###############################################################################################################
                # VERY IMPORT STEP TO GET SAVE THE .FIT FILES
                # print("renaming bin_"+str(binNum)+".fit to seed-analysis/bin_"+str(binNum)+"-seed"+str(seedShifted)+".fit")
                print("Moving seed"+str(seedShifted)+" related files")
                os.rename("bin_"+str(binNum)+".fit","seed-analysis/bin_"+str(binNum)+"-seed"+str(seedShifted)+".fit")
                os.rename("bin_"+str(binNum)+".ni","seed-analysis/bin_"+str(binNum)+"-seed"+str(seedShifted)+".ni")
                os.rename("bin_"+str(binNum)+"-seed"+str(seedShifted)+".cfg","seed-analysis/bin_"+str(binNum)+"-seed"+str(seedShifted)+".cfg")
                ##############################################################################################################
                break
            if (num_lines==0):
                print("Bin_"+str(binNum)+" did not converge with this set of initialization!")
            if (num_lines==0 and j==(maxIter-1)):
                print("Bin_"+str(binNum)+" did not converge maxIter times!")
                raise ValueError("Bin_"+str(binNum)+" did not converge maxIter times! Have to rethink this!")

stop = time.time()
print("Execution time in seconds: %s" % (stop-start))
