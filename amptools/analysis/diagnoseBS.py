import pandas as pd
import numpy as np
from loadCfg import loadCfg
import matplotlib
matplotlib.use('agg') # TKinter might not be installed so use agg
import matplotlib.pyplot as plt
import os
import sys
import shutil

#########
# Load the command line arguments
#########
def printHelp():
    print("Draws diagnostic plots for bootstrapped fits\n--------------------------------------------")
    print("Usage: python diagnoseBS.py <amp> <folder> <fileName> <uncertColName> <delim> <accumlateRunning>")
    print("This will look for <fileName>-seedN.out files in <folder> directory and build diagnostic plots for <amp>")
    print("   The data will be read in with a <delim> and error columns are taken to contain the word <uncertColName>")
    print("    if we dont care about calculating the running std we can speed things up by seeting false to flag <accumlateRunning>")
    print('i.e. python diagnoseBS.py S0+ plot_etapi_delta_results plot_etapi_delta err "\\t" True')
    print('i.e. python diagnoseBS.py H0_10 plot_etapi_delta_results project_moments_polarized uncert " " True') 

argc=len(sys.argv)
print("\n")
if argc == 7:
    chosenAmp=sys.argv[1]
    folderName=sys.argv[2]
    fileName=sys.argv[3]
    uncertColName=sys.argv[4]
    delim=sys.argv[5]
    accumlateRunning=sys.argv[6]
else:
    printHelp()
    exit()

    
#########
# Recreate the output folder and define some filetypes 
#########
cfg=loadCfg()
nBS=cfg["NUMBER_SEEDS"]

outputFolderName="bootstrapDiagnostics"
if not os.path.exists(outputFolderName):
    os.makedirs(outputFolderName)
#else:
#    shutil.rmtree(outputFolderName)
#    os.makedirs(outputFolderName)

# some more specifications for the types and delimiter
outFileType=".png"
fileType=".out"


#########
# Load the data. We can save computational time if we just do float32 instead of the default 64
#########
data=pd.read_csv("/".join([folderName,fileName+fileType]),delimiter=delim,engine='python')
data=data.astype(np.float32)
nonErrCols=[tmp for tmp in data.columns if tmp.find(uncertColName)==-1]
errCols=[tmp for tmp in data.columns if tmp.find(uncertColName)!=-1]
assert chosenAmp in data.columns



#########
# Define the functions to do our extraction of bootstrapped results
#########
def extractBSerrs(dataBS):
    '''
    Takes in the dataframe and drops all the columns that hold the amptools errors
    uses groupby+std to recompute the stds and renames the columns
    '''
    dataBS=dataBS.drop("seed",axis=1)
    dataBS=dataBS.groupby("M").std().reset_index()
    for ampErr in errCols:
        dataBS=dataBS.drop(ampErr,axis=1)
    dataBS.columns=["M"]+errCols
    return dataBS

def getBSerrors(folderName,nBS,accumlateRunning):
    '''
    Sets up the dataframes that can then be passed to extractBSerrs to calculate
    the std for the dataset
    dataBSerr: Final error dataframe
    dataBS_errRunning: A list that is nBS long containing dataframes that calculated a running std
    dataBS_raw: Giant dataframe containing all the raw amplitudes from all bootstrapped iterations
    '''
    dataBSs=[]
    dataBS_errRunning=[]
    for iSeed in range(nBS):
        if iSeed%5==0:
            print("Current seed {0}".format(iSeed))
        dataBS=pd.read_csv(folderName+"/"+fileName+"-seed"+str(iSeed)+fileType, delimiter=delim,engine='python')
        dataBS["seed"]=iSeed
        dataBSs.append(dataBS)
        if iSeed>0 and accumlateRunning: # cant calculate a std with 1 element when calculating sample std --- normalized by n-1
            tempDF=extractBSerrs(pd.concat(dataBSs,ignore_index=True))
            tempDF["seed"]=iSeed
            dataBS_errRunning.append(tempDF)
    dataBS_raw=pd.concat(dataBSs,ignore_index=True)  
    dataBSerr=extractBSerrs(dataBS_raw)
    return dataBSerr, dataBS_errRunning, dataBS_raw


## LOAD THE ACC DATA AND DO THE BOOTSTRAPPING TO REPLACE THE GIVEN ERRS
dataBS, dataBS_errRunning, dataBS_raw = getBSerrors(folderName,nBS,accumlateRunning)
data[errCols]=dataBS.drop("M",axis=1)


if accumlateRunning:
    #########
    # Output two diagnostic plots. One tracks the running BS error as a function nBS we include
    # The other tracks the distribution of the resample BS amplitudes
    #########
    print("\nDrawing diagnostic histograms")
    
    #################
    # Has bootstrapped errors stabilized and is std a good estimate of distributions spread?
    # So what does this plot tell us? 
    # We are looking for a plateau so 
    # we know the std has stabilized
    #################
    allrunningBS=pd.concat(dataBS_errRunning,ignore_index=True)
    nMasses=len(allrunningBS.M.unique())
    
    col=[tmp for tmp in allrunningBS.columns if tmp.find(chosenAmp)!=-1 and tmp.find(uncertColName)!=-1]
    if len(col)!=1:
        raise ValueError("could not find uncertainty column")
    else:
        col=col[0]
    fig,axes=plt.subplots(5,13,figsize=(16,8))
    fig.subplots_adjust(top=0.8)
    suptitle=fig.suptitle("Bootstrapped "+col+" vs nBS (max="+str(nBS)+")",size=20,y=1.05)
    axes=axes.flatten()
    for imass, mass in enumerate(allrunningBS.M.unique()):
        df=allrunningBS[allrunningBS.M==mass]
        axes[imass].plot(df.seed,df[col])
        axes[imass].axis('off')
        axes[imass].set_title(mass)
    plt.tight_layout()
    # When you save the fig, add the suptitle text object as an extra artist
    plt.savefig(outputFolderName+"/"+"runningBSerr_"+chosenAmp+outFileType, 
                bbox_extra_artists=(suptitle,), bbox_inches="tight")
    
#################
# What are we looking for here?
# We are looking to see if std is 
# a good measure of the spread of 
# the distribution. It is if it
# looks ~ Gaus.
#################

col=chosenAmp
fig,axes=plt.subplots(5,13,figsize=(16,8))
suptitle=fig.suptitle("Distibution of Bootstrapped "+col,size=20,y=1.05)
axes=axes.flatten()
for imass, mass in enumerate(dataBS.M):
    df=dataBS_raw[dataBS_raw.M==mass]
    axes[imass].hist(df[col],bins=20)
    axes[imass].axis('off')
    axes[imass].set_title(mass)
plt.tight_layout()
plt.savefig(outputFolderName+"/"+"distributionBS_"+chosenAmp+outFileType,
           bbox_extra_artists=(suptitle,), bbox_inches="tight")


#################
# Output the data to a csv we can use later
#################
print("Writing bootstrapped dataset to csv")
data.to_csv("bootstrapDiagnostics/"+fileName+"_BS.csv",index=False)
