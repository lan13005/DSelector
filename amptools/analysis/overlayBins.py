from loadCfg import loadCfg
import os
import sys
import subprocess
import time
import shutil

def printHelp():
    print("\nUsage: \n-----------------")
    print("Running etapi_plotter and plotting their results into a giant pdf")
    print("python overlayBins.py x y")
    print("where x=0 -- for all bins run etapi_plotter")
    print("      x=1 -- gather all results from etapi_plotter into overlayPlots folder")
    print("      x=2 -- do both")
    print('where y is a string -- Set to "" to do nothing')
    print("      y is a string containing _ separated amplitudes to plot that are again ; separated to group another set")
    print('      i.e "S0+_D0+;S0+_D0+_P1+" will run etapi_plotter twice. One that plots the S0+ and D0+ contributions only.')
    print("      The second will include P0+ contribution also")
    

args=sys.argv
nargs=len(args)
if nargs!=3:
    printHelp()
    exit()
else:
    option=int(args[1])
    ampString=str(args[2])
    groupVec=ampString.split(";")
    groups=[tmp if tmp!="" else tmp for tmp in groupVec]
    groupTags=["_"+tmp if tmp!="" else tmp for tmp in groupVec]
    groupsWithQuotes=['"_'+tmp+'"' if tmp!="" else '""' for tmp in groupVec]
    print("GROUPS: ")
    print(groups)

cfg=loadCfg()
fitName=cfg["FIT_DIR"]
lowMass=cfg["LOW_MASS"]
uppMass=cfg["UPP_MASS"]
nBins=cfg["NUMBER_BINS"]
cwd=os.getcwd()
fitLoc=cwd+"/"+fitName+"/"
print("Fit location: "+fitLoc)

verbose=False

def runEtaPiPlotterForAllBins():
    for iBin in range(nBins):
        binLoc=fitLoc+"bin_"+str(iBin)
        print("moving to "+binLoc)
        os.chdir(binLoc)
        print("   running etapi_plotter...")
        for igroup in range(len(groups)):
            cmd=["etapi_plotter","bin_"+str(iBin)+"-full.fit","-o","etapi_plot"+groupTags[igroup]+".root","-s",groups[igroup]]
            print("calling: "+" ".join(cmd))
            if verbose:
                subprocess.check_call(cmd)
            else:
                with open(os.devnull, 'wb') as devnull:
                    subprocess.check_call(cmd, stdout=devnull, stderr=subprocess.STDOUT)
    time.sleep(2)

def gatherPlotResultsIntoPDFs():
    sedArgs=["sed","-i",'s@nBins=.*;@nBins='+str(nBins)+';@g',"overlayBins.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@lowMass=.*;@lowMass='+str(lowMass)+';@g',"overlayBins.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@uppMass=.*;@uppMass='+str(uppMass)+';@g',"overlayBins.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    sedArgs=["sed","-i",'s@fitName=.*;@fitName="'+fitName+'";@g',"overlayBins.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    repStr="vector<string> groups="
    vecStr="{"
    vecStr+=",".join(groupsWithQuotes)
    vecStr+="}"
    print("replacing vector to "+vecStr)
    sedArgs=["sed","-i",'s@'+repStr+'.*;@'+repStr+vecStr+';@g',"overlayBins.C"]
    subprocess.Popen(sedArgs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).wait()
    if verbose:
        subprocess.check_call(["root","-l","-b","-q","overlayBins.C"])
    else:
        with open(os.devnull, 'wb') as devnull:
            subprocess.check_call(["root","-l","-b","-q","overlayBins.C"], stdout=devnull, stderr=subprocess.STDOUT)

if option==0 or option==2:
    print("RUNNING etapi_plotter FOR ALL BINS...")
    runEtaPiPlotterForAllBins()
if option==1 or option==2:
    print("Recreating overlayPlots folder...")
    if os.path.exists("overlayPlots") and os.path.isdir("overlayPlots"):
        shutil.rmtree("overlayPlots")
    os.mkdir("overlayPlots")
    print("GATHERING ALL THE RESULTS INTO A PDF...")
    gatherPlotResultsIntoPDFs()

print("DONE!")
