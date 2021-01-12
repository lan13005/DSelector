import os
from loadCfg import loadCfg
import shutil

cwd=os.getcwd()
cfg=loadCfg()
fitDir=cfg["FIT_DIR"]
nBins=cfg["NUMBER_BINS"]
nSeeds=cfg["NUMBER_SEEDS"]

fitDir=cwd+"/"+fitDir
outDir="plot_etapi_delta_results"
try:
    os.mkdir(outDir)
    print("Made "+outDir)
except:
    shutil.rmtree(outDir)
    os.mkdir(outDir)
    print("remaking "+outDir)



# Run plot etapi delta on the unbootstrapped fits
def setFitFile(dirTag,fitTag):
    for iBin in range(nBins):
        binName="bin_"+str(iBin)
        baseName=fitDir+"/"+binName+"/"
        #print("Moving "+baseName+dirTag+binName+fitTag+".fit  to  "+baseName+binName+".fit")
        shutil.move(baseName+dirTag+binName+fitTag+".fit", baseName+binName+".fit")
def fixFitFile(dirTag,fitTag):
    for iBin in range(nBins):
        binName="bin_"+str(iBin)
        baseName=fitDir+"/"+binName+"/"
        #print("Moving "+baseName+binName+".fit  to  "+baseName+dirTag+binName+fitTag+".fit")
        shutil.move(baseName+binName+".fit", baseName+dirTag+binName+fitTag+".fit")
setFitFile("","-full")
os.system("plot_etapi_delta -o "+outDir+"/plot_etapi_delta.out") 
fixFitFile("","-full")
    

for iSeed in range(nSeeds):
    setFitFile("seed-analysis/","-seed"+str(iSeed))
    os.system("plot_etapi_delta -o "+outDir+"/plot_etapi_delta-seed"+str(iSeed)+".out") 
    fixFitFile("seed-analysis/","-seed"+str(iSeed))



    ## Rename original file and move the seed file in
    #for iBin in range(nBins):
    #    binName="bin_"+str(iBin)
    #    baseDir=fitDir+"/"+binName+"/"+binName

    #    file2=fitDir+"/"+binName+"/seed-analysis/"+binName+"-seed"+str(iSeed)+".fit"
    #    shutil.move(file2,srcFile)
    #
    ## Run plot eta pi delta again
    #os.system("plot_etapi_delta -o "+outDir+"/plot_etapi_delta-seed"+str(iSeed)+".out") 
    #
    ## move the seed file out and fix the original file's name
    #for iBin in range(nBins):
    #    binName="bin_"+str(iBin)
    #    baseDir=fitDir+"/"+binName+"/"+binName
    #    srcFile=baseDir+binName+".fit"
    #    shutil.move(srcFile,file2)


