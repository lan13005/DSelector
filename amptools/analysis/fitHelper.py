import random

def getAmps(ampNames,binDir,seedFile,rndSamp,percentDeviation,baseDeviation,verbose=True):
    '''
    This function gets the amplitudes for a specific bin based on some seedFile. 
    The seed file will be the converged amplitudes from the previous bin for the regular fits and
        the converged amplitudes of the regular fits for the bootstrapped fits
    Takes extra arguments to sample around these converged values (percentDeviation and baseDeviation)
    rndSamp is a flag to tell the program to randomly sample in the first place
    '''
    if(rndSamp==False):
        print("We will use the converged values instead of randomly sampling around them!")
        percentDeviation=0
        baseDeviation=0
    else:
        print("using percentDeviation as: "+str(percentDeviation))
        print("using baseDeviation as: "+str(baseDeviation))
    realAmps=[]
    imAmps=[]
    isreal=[]
    isAmpReal=False
    print("Opening: "+binDir+"/"+seedFile)
    with open(binDir+"/"+seedFile,"r") as param_init_cfg:
        paramLines = param_init_cfg.readlines()
        lineCounter=0
        ampCounter=-1
        for paramLine in paramLines:
            if lineCounter%2==0:
                ampCounter+=1
                realAmp = float(paramLine.split(" ")[3])
                imAmp = float(paramLine.split(" ")[4])
                lenParams = len(paramLine.split(" "))
                real_low = float(realAmp*(1-percentDeviation)-baseDeviation)
                real_up = float(realAmp*(1+percentDeviation)+baseDeviation)
                im_low = float(imAmp*(1-percentDeviation)-baseDeviation)
                im_up = float(imAmp*(1+percentDeviation)+baseDeviation)
                realSampled = str(random.uniform(real_low,real_up))
                if (lenParams>5) and (paramLine.split(" ")[5].rstrip()=="real"):
                    imSampled="0.0"
                    isAmpReal=True
                else: 
                    imSampled= str(random.uniform(im_low,im_up))
                    isAmpReal=False
            if(verbose):
                print("Below is paramLists, check if below values match:")
                print(paramLine.split(" "))
                print("Amps from "+seedFile+": "+ampNames[ampCounter])
                print("realAmp: "+str(realAmp))
                print("imAmp: "+str(imAmp))
                print("Sampling from the below ranges")
                print("real_low, real_up: "+str(real_low)+", "+str(real_up))
                print("im_low, im_up: "+str(im_low)+", "+str(im_up))
                print("realSampled, imSampled: {0},{1}".format(realSampled,imSampled))
                print("----------------------------")
            realAmps.append(realSampled)
            imAmps.append(imSampled)
            isreal.append(isAmpReal)
            #print("lc, ra, ia: "+str(lineCounter)+", "+str(realSampled)+", "+str(imSampled))
            lineCounter+=1
        if (len(realAmps)==len(imAmps)==len(isreal)):
            print("length of realAmps, imAmps, isreal is good!")
        else:
            print("length of realAmps: "+str(len(realAmps)))
            print("length of imAmps: "+str(len(imAmps)))
            print("length of isreal: "+str(len(isreal)))
            raise ValueError("length of realAmps, imAmps, isreal is not good!")
        #print("==================================\n\n")
    return realAmps, imAmps, isreal
