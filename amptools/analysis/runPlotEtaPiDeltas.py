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
    os.rmdir(outDir)
    os.mkdir(outDir)
    print("remaking "+outDir)



for iSeed in range(nSeeds):
    os.system("plot_etapi_delta -o "+outDir+"/plot_etapi_delta.out") 
    
    # Rename original file and move the seed file in
    for iBin in range(nBins):
        binName="bin_"+str(iBin)
        orig=fitDir+"/"+binName+"/"+binName+".fit"
    
        shutil.move(orig,orig+".tmp")
        file2=fitDir+"/"+binName+"/seed-analysis/"+binName+"-seed"+str(iSeed)+".fit"
        shutil.move(file2,orig)
    
    os.system("plot_etapi_delta -o "+outDir+"/plot_etapi_delta-seed"+str(iSeed)+".out") 
    
    # move the seed file out and fix the original file's name
    for iBin in range(nBins):
        binName="bin_"+str(iBin)
        orig=fitDir+"/"+binName+"/"+binName+".fit"
    
        shutil.move(orig,file2)
        shutil.move(orig+".tmp",orig)


