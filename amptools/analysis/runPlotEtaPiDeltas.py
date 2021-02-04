import os
from loadCfg import loadCfg
import shutil
import sys


def printHelp():
    print("\nUsage: \n-----------------")
    print("Running plot_etapi_delta and/or project_moments for all the bins and gathering their results")
    print("YOU NEED TO MODIFY/COMPILE THE ABOVE PROGRAMS TO YOUR OWN WAVESET+CONFIG AS THIS IS JUST A DRIVER SCRIPT")
    print("python runPlotEtaPiDeltas.py x b")
    print("where x=0 -- only run plot_etapi_delta")
    print("      x=1 -- only run project_moments")
    print("      x=2 -- run both plot_etapi_delta and project_moments")
    print("where b=1 -- run programs with bootstrapped results also")
    print("      b=0 -- do not run programs with bootstrapped results")

args=sys.argv
nargs=len(args)
if nargs!=3:
    printHelp()
    exit()
else:
    option=int(args[1])
    runBS=bool(int(args[2]))

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
        try:
            shutil.move(baseName+dirTag+binName+fitTag+".fit", baseName+binName+".fit")
        except:
            pass
def fixFitFile(dirTag,fitTag):
    for iBin in range(nBins):
        binName="bin_"+str(iBin)
        baseName=fitDir+"/"+binName+"/"
        #print("Moving "+baseName+binName+".fit  to  "+baseName+dirTag+binName+fitTag+".fit")
        try:
            shutil.move(baseName+binName+".fit", baseName+dirTag+binName+fitTag+".fit")
        except:
            pass
setFitFile("","-full")
if option==0:
    os.system("plot_etapi_delta -o "+outDir+"/plot_etapi_delta.out") 
if option==1:
    os.system("project_moments_polarized -o "+outDir+"/project_moments_polarized.out") 
if option==2:
    os.system("plot_etapi_delta -o "+outDir+"/plot_etapi_delta.out") 
    os.system("project_moments_polarized -o "+outDir+"/project_moments_polarized.out") 
fixFitFile("","-full")
    

if runBS:
    for iSeed in range(nSeeds):
        setFitFile("seed-analysis/","-seed"+str(iSeed))
    
        ## Rename original file and move the seed file in
        #for iBin in range(nBins):
        #    binName="bin_"+str(iBin)
        #    baseDir=fitDir+"/"+binName+"/"
        #    srcFile=baseDir+binName+".fit"

        #    file2=fitDir+"/"+binName+"/seed-analysis/"+binName+"-seed"+str(iSeed)+".fit"
        #    shutil.move(file2,srcFile)
        
        # Run plot eta pi delta again
        if option==0:
            os.system("plot_etapi_delta -o "+outDir+"/plot_etapi_delta-seed"+str(iSeed)+".out") 
        if option==1:
            os.system("project_moments_polarized -o "+outDir+"/project_moments_polarized-seed"+str(iSeed)+".out") 
        if option==2:
            os.system("plot_etapi_delta -o "+outDir+"/plot_etapi_delta-seed"+str(iSeed)+".out") 
            os.system("project_moments_polarized -o "+outDir+"/project_moments_polarized-seed"+str(iSeed)+".out") 
        
        ## move the seed file out and fix the original file's name
        #for iBin in range(nBins):
        #    binName="bin_"+str(iBin)
        #    baseDir=fitDir+"/"+binName+"/"
        #    srcFile=baseDir+binName+".fit"

        #    file2=fitDir+"/"+binName+"/seed-analysis/"+binName+"-seed"+str(iSeed)+".fit"
        #    shutil.move(srcFile,file2)
        fixFitFile("seed-analysis/","-seed"+str(iSeed))


