import os
import sys
import subprocess
import time
import fileinput
import random
from loadCfg import loadCfg

if len(sys.argv)!=4:
	raise ValueError("Must take iProcess, binsPerProc, processSpawned as an argument!")
else:
	iProcess = int(sys.argv[1])
	binsPerProc = int(sys.argv[2])
	processSpawned = int(sys.argv[3])
#realSpace_low = [-50 , -30 ,  50 , 200, -40, 200,  20]
#realSpace_up  = [ 50 ,  60 ,  200, 500,  20, 400,  40]
#imSpace_low   = [-110, -70 , -30 , 0  , -30, 0  , -40]
#imSpace_up    = [-50 , -10 ,  30 , 0  ,  30, 0  ,  20] 
#fitName = "divideRoot"
#seedFile = "param_init.cfg"
#maxIter=10
#ampNames=["S0+","P0+","P1+","D0+","D1+","P1-","D1-"]

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

print("realSpace_low: "+str(realSpace_low))
print("realSpace_up: "+str(realSpace_up))
print("imSpace_low: "+str(imSpace_low))
print("imSpace_up: "+str(imSpace_up))
print("fitName: "+str(fitName))
print("seedFile: "+str(seedFile))
print("maxIter: "+str(maxIter))
print("ampNames: "+str(ampNames))




# this is what bins we consider for this process. 
def getBinsArray(i,bpt,threads):
	#iProcess, binsPerThrea,Threads
        array1 = [i]
	bptm1 = bpt-1
        for iproc in range(bptm1):
                array1.append(i+(iproc+1)*threads)
        return array1
binsArray = getBinsArray(iProcess,binsPerProc,processSpawned)

start = time.time()

workingDir = os.getcwd()
print("\n")
print("current working directory: %s" % (workingDir))
fitDir = workingDir+"/"+fitName
print("fit directory: %s" % (fitDir))


random.seed(seedAmpInit)
def getAmps():
	realAmps=[]
	imAmps=[]
	isreal=[]
	for i in range(len(realSpace_low)):
		real_low = int(realSpace_low[i])
		real_up = int(realSpace_up[i])
		im_low = int(imSpace_low[i])
		im_up = int(imSpace_up[i])
		realAmps.append(str(random.uniform(real_low,real_up)))
		if ((im_up-im_low)==0):
			imAmps.append("0.0")
			isreal.append(True)
		else: 
			imAmps.append(str(random.uniform(im_low,im_up)))
			isreal.append(False)
	return realAmps, imAmps, isreal
def writeCfg(ampNames,binNum):
	ampCounter=0
	realAmps, imAmps, isreal = getAmps()
	for line in fileinput.input(fitDir+"/bin_"+str(binNum)+"/bin_"+str(binNum)+".cfg",inplace=1):
		if line.strip().startswith("initialize"):
			posNegTag=line.split(" ")[1].split("::")[1]
			if (isreal[ampCounter]):
				line = "initialize EtaPi::"+posNegTag+"::"+ampNames[ampCounter]+" cartesian "+realAmps[ampCounter]+" "+imAmps[ampCounter]+" real\n"
			else:
				line = "initialize EtaPi::"+posNegTag+"::"+ampNames[ampCounter]+" cartesian "+realAmps[ampCounter]+" "+imAmps[ampCounter]+"\n"
			
			ampCounter+=1
		sys.stdout.write(line)
	return True


for binNum in binsArray:
	# this is really the step where it requires us to use different processes instead of one process using multiple threads
	os.chdir(fitDir)
	os.chdir("bin_"+str(binNum))
	print("================   bin_%i    =================" % (binNum))
	print("current working directory: %s" % (os.getcwd()))
	for j in range(maxIter):
		writeCfg(ampNames,binNum)
		# creating a file so that I could dump the results of the fit into it
		print("----------------------------------------------------------")
		callFit = "fit -c bin_"+str(binNum)+".cfg -s param_init.cfg"
		binName = "bin_"+str(binNum)+".cfg"
		print(callFit)
		print("Trying with seed="+str(seedAmpInit)+" at iteration="+str(j))
		#DEVNULL = open(os.devnull, 'wb')
		logFile = open("bin_"+str(binNum)+"-full.log","w+")
		subprocess.Popen(callFit.split(" "), stdout=logFile).wait()
		#num_lines = sum(1 for line in open('param_init.cfg','r'))
		num_lines=0
		print(subprocess.check_output(['grep','initialize',binName]).rstrip())
		with open("param_init.cfg","r") as param_init_cfg:
			param_lines = param_init_cfg.readlines()
			for param_line in param_lines:
				num_lines+=1
		print("num_lines: "+str(num_lines)+" should equal the total number of waves!")
		if (num_lines>0):
			print("Bin_"+str(binNum)+" converged with this set of initialization!")
			###############################################################################################################
			# VERY IMPORT STEP TO GET SAVE THE .FIT FILE SO THAT PROEJCT MOMENTS CAN USE THE ONE THAT IS NOT BOOTSTRAPPPED
			print("renaming bin_"+str(binNum)+".fit to bin_"+str(binNum)+"-full.fit")
			os.rename("bin_"+str(binNum)+".fit","bin_"+str(binNum)+"-full.fit")
			###############################################################################################################
			break
		if (num_lines==0 and j==maxIter-1):
			raise ValueError("Bin_"+str(binNum)+" did not converge maxIter times! Have to rethink this!")

stop = time.time()
print("Execution time in seconds: %s" % (stop-start))
