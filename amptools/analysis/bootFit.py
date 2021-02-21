import os
import glob
import sys
import subprocess
import time
import fileinput
import random
import shutil
import glob
from loadCfg import loadCfg
from fitHelper import updateCfg, rndSampInits, setupBootstrap

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
    seeds.append(random.randint(0,1E5))
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

##############################
# Make all the subdirectoires
##############################
for binNum in binsArray:
    binDir=fitDir+"/bin_"+str(binNum)
    if os.path.exists(binDir+"/seed-analysis") and os.path.isdir(binDir+"/seed-analysis"):
        shutil.rmtree(binDir+"/seed-analysis")
        os.mkdir(binDir+"/seed-analysis")
    else:
        os.mkdir(binDir+"/seed-analysis")

##############################
# Run fit for all bins and all seeds
##############################
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

        # Copy the config file over and rename the fit name
        print("Copying "+bin_cfg_seed+" from "+bin_cfg+"\n--------------------------------------------------")
        shutil.copyfile(bin_cfg, bin_cfg_seed)
        bin_cfg_base=bin_cfg.split(".")[0]
        bin_cfg_seed_base=bin_cfg_seed.split(".")[0]
        sedCall="sed -i s/"+bin_cfg_base+"/"+bin_cfg_seed_base+"/g "+bin_cfg_seed
        print(sedCall)
        subprocess.Popen(sedCall.split(" ")).wait()
        # Setup the bootstrap datareader
        setupBootstrap(bin_cfg_seed,seed)

        rndSamp=rndSamp_flag
        for j in range(maxIter):
            print("======================== NEW ATTEMPT =========================\n=======================================================")
            if (j<rampIter):
                baseDeviation=baseDeviation1
            else:
                baseDeviation=baseDeviation2
            if rndSamp:
                rndSampInits(bin_cfg_seed,percentDeviation,baseDeviation)
            print("First 5 current initial values for amplitudes:")
            print("----------------------------------")
            firstCoupleLines=subprocess.check_output(['grep','^initialize',bin_cfg_seed]).decode("utf-8").split("\n")[:5]
            print("\n".join(firstCoupleLines))
            print("----------------------------------")
            rndSamp=True # next iteration will start randomly sampling 

            # creating a file so that I could dump the results of the fit into it
            callFit = "fit -c "+bin_cfg_seed+" -s "+param_init_seed
            print(callFit)
            logFileName="seed-analysis/bin_"+str(binNum)+"-seed"+str(seedShifted)+".log"
            with open(logFileName,"w+") as logFile:
                subprocess.Popen(callFit.split(" "), stdout=logFile).wait()
            num_lines=0
            if os.path.isfile(param_init_seed):	
                with open(param_init_seed,"r") as param_init_cfg:
                    param_lines = param_init_cfg.readlines()
                    for param_line in param_lines:
                        if param_line[:10]=="initialize":
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
                niFiles=glob.glob("*ni")
                for niFile in niFiles:
                    if "full" not in niFile: # recall that the "full" tag was used for the fit to the data, no resampling
                        fileName=niFile.split(".")[0]
                        fileType=niFile.split(".")[1]
                        with open(niFile,"r") as ni:
                            lines=ni.readlines()
                            for line in lines:
                                if "nan" in line:
                                    # Normalization integrals is related to the acceptance
                                    raise ValueError(niFile+" contains nan values! gen matrix is on top and acc is on bottom")
                        os.rename(niFile,"seed-analysis/"+niFile+"-seed"+str(seedShifted)+".ni")
                ##############################################################################################################
                break
            if (num_lines==0):
                print("Bin_"+str(binNum)+" did not converge with this set of initialization! Will try to resample if we can...")
            if (num_lines==0 and j==(maxIter-1)):
                print("Bin_"+str(binNum)+" did not converge maxIter times!")
                raise ValueError("Bin_"+str(binNum)+" did not converge maxIter times! Have to rethink this!")

stop = time.time()
print("Execution time in seconds: %s" % (stop-start))
