import random
import numpy as np

def updateCfg(source,destination):
    '''
    This function reads in param_init.cfg located at <source> path and dumps it into the <destination> cfg
    file. It will replace the block of "initialize" lines and replace it directly with the data from param_init.cfg

    (Note)param_init.cfg has a weird parameter line at the bottom that does not have a newline separator
    '''
    ## We need to open and read the entire config file into memory. We can then close and begin writing to it
    initializing=False
    with open(destination,"r") as dest:
        lines = dest.readlines()
        nlines = len(lines)
        nInitLines=0
        for line in lines:
            if line[:10]=="initialize":
                nInitLines+=1 

    ## Grab the initialization parameters from param_init.cfg
    print("Updating config file")
    initLines=[]
    with open(source,"r") as src:
        lines = src.readlines()
        for line in lines:
            if line[:10]=="initialize":
                initLines.append(line)
    with open(destination,"w") as dest:
        for lineNum,line in enumerate(lines):
            # If we are finished initializing
            if (line[:10]!="initialize" or lineNum==nlines-1) and initializing: 
                print("  End of initialization block... dumping param_init.cfg")
                linesKept=0
                for initLine in initLines and linesKept<nInitLines:
                    dest.write(initLine)
                    linesKept+=1
                dest.write("\n")
                initializing=False
            # We dont write if we are still initializing
            if line[:10]=="initialize":
                initializing=True
            # leave the non-initialize lines alone
            else:
                dest.write(line.rstrip())
                dest.write("\n")

def rndSampInits(source,percentDeviation,baseDeviation):
    '''
    This function takes in a config file and resample the values based on percentDeviation and baseDeviation
    '''
    mapSampled={}
    print("Randomly sampling initialization paramters")
    # Read in the data and fill a map with the randomly sampled initializations
    with open(source,"r") as src:
        lines = src.readlines()
        for line in lines:
            if line[:10]=="initialize":
                fields=line.split(" ")
                real=np.float64(fields[3])
                imag=np.float64(fields[4])
                if real!=0:
                    real_low = (-1 if real<0 else 1)*real*(1-percentDeviation)-baseDeviation
                    real_up = (-1 if real<0 else 1)*real*(1+percentDeviation)+baseDeviation
                    realSampled = random.uniform(real_low,real_up)
                else:
                    realSampled = 0
                if imag!=0:
                    imag_low = (-1 if real<0 else 1)*imag*(1-percentDeviation)-baseDeviation
                    imag_up = (-1 if real<0 else 1)*imag*(1+percentDeviation)+baseDeviation
                    imagSampled= random.uniform(imag_low,imag_up)
                else: 
                    imagSampled = 0
                mapSampled[real]=realSampled
                mapSampled[imag]=imagSampled
    # Write the output 
    with open(source,"w") as src:
        for line in lines:
            if line[:10]=="initialize":
                fields=line.split(" ")
                real=np.float64(fields[3])
                imag=np.float64(fields[4])
                fields[3]=str(mapSampled[real])
                fields[4]=str(mapSampled[imag])
                src.write(" ".join(fields).rstrip())
                src.write("\n")
            else:
                src.write(line.rstrip())
                src.write("\n")


def setupBootstrap(source,seed):
    print("Setting up ROOTDataReaderBootstrap")
    with open(source,"r") as src:
        lines=src.readlines()
    with open(source,"w") as src:
        for line in lines:
            # the line contains ROOTDataReader but is not a comment
            if "ROOTDataReader" in line and line[0]!="#":
                fields=line.split(" ")
                fields[2]="ROOTDataReaderBootstrap"
                src.write(" ".join(fields).rstrip()+" "+str(seed))
                src.write("\n")
            else:
                src.write(line.rstrip())
                src.write("\n")
