import os
import sys
import subprocess
import time
import fileinput
import random
from loadCfg import loadCfg
from fitHelper import getAmps

#if len(sys.argv)!=4:
#    raise ValueError("Must take iProcess, binsPerProc, processSpawned as an argument!")
#else:
#    iProcess = int(sys.argv[1])
#    binsPerProc = int(sys.argv[2])
#    processSpawned = int(sys.argv[3])

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
percentDeviation=cfg["PERCENT_DEVIATION"]
baseDeviation1=cfg["BASE_DEVIATION1"]
baseDeviation2=cfg["BASE_DEVIATION2"]
rampIter=cfg["RAMP_DEVIATION_ITER"]
# this flag will determine if we randomly sample around the converged amplitudes or use the converged amplitudes as a starting point for the bootstrapped dataset	
rndSamp_flag=cfg["RND_SAMP"]

print("realSpace_low: "+str(realSpace_low))
print("realSpace_up: "+str(realSpace_up))
print("imSpace_low: "+str(imSpace_low))
print("imSpace_up: "+str(imSpace_up))
print("fitName: "+str(fitName))
print("seedFile: "+str(seedFile))
print("maxIter: "+str(maxIter))
print("ampNames: "+str(ampNames))


# this is what bins we consider for this process. 
# def getBinsArray(i,bpt,threads):
#     #iProcess, binsPerThrea,Threads
#     array1 = [i]
#     bptm1 = bpt-1
#     for iproc in range(bptm1):
#         array1.append(i+(iproc+1)*threads)
#     return array1
# binsArray = getBinsArray(iProcess,binsPerProc,processSpawned)
binsArray = range(nBins)

start = time.time()

workingDir = os.getcwd()
print("\n")
print("current working directory: %s" % (workingDir))
fitDir = workingDir+"/"+fitName
print("fit directory: %s" % (fitDir))


random.seed(seedAmpInit)
def getAmpsRandomly():
    '''
    Randomly samples to initialze this bin. Using zlm amplitudes requires filling twice, once 
    for Im and once for Re parts. 
    '''
    realAmps=[]
    imAmps=[]
    isreal=[]
    for i in range(len(realSpace_low)):
        real_low = int(realSpace_low[i])
        real_up = int(realSpace_up[i])
        im_low = int(imSpace_low[i])
        im_up = int(imSpace_up[i])
        rndVal=random.uniform(real_low,real_up)
        realAmps.append(str(rndVal))
        realAmps.append(str(rndVal))
        if ((im_up-im_low)==0):
            imAmps.append("0.0")
            imAmps.append("0.0")
            isreal.append(True)
            isreal.append(True)
        else: 
            rndVal=random.uniform(im_low,im_up)
            imAmps.append(str(rndVal))
            imAmps.append(str(rndVal))
            isreal.append(False)
            isreal.append(False)
    return realAmps, imAmps, isreal



def writeCfg(ampNames,binNum,rndSamp,percentDeviation,baseDeviation):
    lineCounter=0
    ampCounter=-1 
    prevBinDir=fitDir+"/bin_"+str(binNum-1) # grabbing the converged amplitudes from the previous bin
    if binNum==0:
        realAmps, imAmps, isreal = getAmpsRandomly()
    else: 
        realAmps, imAmps, isreal = getAmps(ampNames,prevBinDir,seedFile,rndSamp,percentDeviation,baseDeviation,False)
    for line in fileinput.input(fitDir+"/bin_"+str(binNum)+"/bin_"+str(binNum)+"-full.cfg",inplace=1):
        # using fileinput seems a bit trickier. there can be some encoding issues and newline assumptions
        # need to use write instead of print and strip the lines
        if line.strip().startswith("initialize"):
            ampCounter+=1
            posNegTag=line.split(" ")[1].split("::")[1]
            ampName=line.split(" ")[1].split("::")[2]
            if (isreal[ampCounter]):
                line = "initialize EtaPi::"+posNegTag+"::"+ampName+" cartesian "+realAmps[ampCounter]+" "+imAmps[ampCounter]+" real\n"
                sys.stdout.write(line)
            else:
                line = "initialize EtaPi::"+posNegTag+"::"+ampName+" cartesian "+realAmps[ampCounter]+" "+imAmps[ampCounter]+"\n"
                sys.stdout.write(line)
            lineCounter+=1

        else:
            sys.stdout.write(line.lstrip().rstrip()+"\n")
    return True


for binNum in binsArray:
    # this is really the step where it requires us to use different processes instead of one process using multiple threads
    os.chdir(fitDir)
    os.chdir("bin_"+str(binNum))
    print("\n================   bin_%i    =================" % (binNum))
    print("current working directory: %s" % (os.getcwd()))
    rndSamp=rndSamp_flag
    for j in range(maxIter):
        if (j<rampIter):
            baseDeviation=baseDeviation1
        else:
            baseDeviation=baseDeviation2
        writeCfg(ampNames,binNum,rndSamp,percentDeviation,baseDeviation)
        rndSamp=True # next iteration will start randomly sampling 
        # creating a file so that I could dump the results of the fit into it
        print("----------------------------------------------------------")
        callFit = "fit -c bin_"+str(binNum)+"-full.cfg -s "+seedFile
        binName = "bin_"+str(binNum)+"-full.cfg"
        print(callFit)
        print("Trying with seed="+str(seedAmpInit)+" at iteration="+str(j))
        #DEVNULL = open(os.devnull, 'wb')
        logFile = open("bin_"+str(binNum)+"-full.log","w+")
        subprocess.Popen(callFit.split(" "), stdout=logFile).wait()
        #num_lines = sum(1 for line in open('param_init.cfg','r'))
        num_lines=0
        print("Initial values for the amplitudes:")
        print("----------------------------------")
        print(subprocess.check_output(['grep','initialize',binName]).rstrip())
        print("----------------------------------")
        with open("param_init.cfg","r") as param_init_cfg:
            param_lines = param_init_cfg.readlines()
            for param_line in param_lines:
                num_lines+=1
        print("num_lines: "+str(num_lines)+" should equal the total number of waves!")
        if (num_lines>0):
            print("Bin_"+str(binNum)+" converged with this set of initialization! -- Iteration:"+str(j))
            ###############################################################################################################
            # VERY IMPORT STEP TO GET SAVE THE .FIT FILE SO THAT PROEJCT MOMENTS CAN USE THE ONE THAT IS NOT BOOTSTRAPPPED
            os.rename("bin_"+str(binNum)+".fit","bin_"+str(binNum)+"-full.fit")
            os.rename("bin_"+str(binNum)+".ni","bin_"+str(binNum)+"-full.ni")
            ##############################################################################################################
            break
        if (num_lines==0 and j==maxIter-1):
            print("Bin_"+str(binNum)+" did not converge maxIter times!")
            raise ValueError("Bin_"+str(binNum)+" did not converge maxIter times! Have to rethink this!")

stop = time.time()
print("Execution time in seconds: %s" % (stop-start))
