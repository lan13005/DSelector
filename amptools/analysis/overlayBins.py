from loadCfg import loadCfg
import os
import subprocess
import time

cfg=loadCfg()
fitName=cfg["FIT_DIR"]
lowMass=cfg["LOW_MASS"]
uppMass=cfg["UPP_MASS"]
nBins=cfg["NUMBER_BINS"]
cwd=os.getcwd()
fitLoc=cwd+"/"+fitName+"/"

verbose=False


def runTwoPiPlotterForAllBins():
    for iBin in range(nBins):
        binLoc=fitLoc+"bin_"+str(iBin)
        print("moving to "+binLoc)
        os.chdir(binLoc)
        print("   running twopiplotter...")
        if verbose:
            subprocess.check_call(["twopi_plotter_mom","bin_"+str(iBin)+".fit"])
        else:
            with open(os.devnull, 'wb') as devnull:
                subprocess.check_call(["twopi_plotter_mom","bin_"+str(iBin)+".fit"], stdout=devnull, stderr=subprocess.STDOUT)
    print("\n\nGathering all twopi_plotter_mom outputs...")
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
    if verbose:
        subprocess.check_call(["root","-l","-b","-q","overlayBins.C"])
    else:
        with open(os.devnull, 'wb') as devnull:
            subprocess.check_call(["root","-l","-b","-q","overlayBins.C"], stdout=devnull, stderr=subprocess.STDOUT)


#runTwoPiPlotterForAllBins()
gatherPlotResultsIntoPDFs()
