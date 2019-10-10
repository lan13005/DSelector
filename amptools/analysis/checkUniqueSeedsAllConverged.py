import os
import glob
from loadCfg import loadCfg

# grab any bin and look at all the seedN-param_init.cfg files and the N values


#nBins=60
cfg=loadCfg()
nBins=cfg["NUMBER_BINS"]
fitDir=cfg["FIT_DIR"]
print("nBins: "+str(nBins))
print("fitDir: "+str(fitDir))


seed_folder = "/seed-analysis/"
allParamInit = glob.glob(fitDir+"/bin_0"+seed_folder+"*param_init.cfg")
print("allParamInits:")
print(allParamInit)
allSeedFiles =[seeds.split("/")[2] for seeds in allParamInit]
# these "seeds" are more like the labels and do not correspond to the actual seeds used in the bootstrapping data
seeds =[seeds.split("/")[3].split("-")[0].split("seed")[1]  for seeds in allParamInit]
print("seeds: ")
print(seeds)
allConverged=True
#create a set that we can check if there are duplicates of seeds used in the bootstrapping process
allSeeds=set()



filesWithNoSeed=set()
for seed in seeds:
	convergedSeed=[]
	fileSeeds = glob.glob(fitDir+"/bin*"+seed_folder+"*seed"+str(seed)+".cfg")
	if (len(fileSeeds)==0):
		raise ValueError("no seed cfgs")
	# seedNset will check if the corresponding seed files have the same seed by checking if this set is of size 1 after looping over all mass bins
	seedNset=set()
	for iFile in range(len(fileSeeds)):
		with open(fileSeeds[iFile],'r') as f:
			for line in f.readlines():
				if 'data  EtaPi' in line:
					try:
						seedNset.add(line.split(" ")[5])
					except:
						filesWithNoSeed.add(fileSeeds[iFile])
	if len(seedNset)!=1:
		print("Seeds in seed"+str(seed)+": ")
		print(list(seedNset))
		raise ValueError("error! seed"+str(seed)+ " does not have the same seeds!")
	else:
		print("seed"+str(seed)+" is setup correctly for all bin!")
		print("Seed = "+str(list(seedNset)[0]))
		if list(seedNset)[0] in allSeeds:
			raise ValueError("LETTING YOU KNOW THAT THERE IS A REPEAT OF SEEDS USED IN THE BOOTSTRAPPING! FIX THIS!")
		else:
			allSeeds.add(list(seedNset)[0])
print("SEEDS ARE ALL UNIQUE! with "+str(len(allSeeds))+" many seeds")

for binNum in range(nBins):	
	for seedNum in seeds:
		#print("seedNum: "+str(seedNum))
		fileLoc = fitDir+"/bin_"+str(binNum)+seed_folder+"seed"+str(seedNum)+"-param_init.cfg"
		if os.path.isfile(fileLoc):
			allConverged=True
		else:
			print(fileLoc+" is empty so it did not converge or the thread got killed so all preceeding bins did not have a chance to run!")
			allConverged=False


print("THESE FILES NEED TO BE FIXED FOR THE SEED INITIALIZATION IN THE INIT FILE:")
for singleFile in filesWithNoSeed:
	print(singleFile)
print("--------------------------------------------------------------------------")



if allConverged:
	print("---------------------------------------------")
	print("All the seeds for all the bins have converged!")	
				

