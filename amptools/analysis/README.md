# Current Propsed Method:
We can directly use plot_etapi_delta and it can output a txt file of all the amplitudes, with normalization integrals included (I think). So, we can still use the old fitting scripts but somehow merge with the plot_etapi_delta program and output the amplitudes. plot.C currently reads in the data from plot_etapi_delta and plots the amplitudes. I HAVE NOT CONTINUED DEVELOPMENT OF THIS CODE IN AWHILE AND THERE IS PROBABLY BUGS. THE BIGGEST THING TO NOTE IS THAT I SIMULATED A DATASET WITH AN a0 AND a2 BUT THESE RESONANCES BOTH SHOW UP IN THE S-WAVE.  

# Old Method
So the goal of this set of scripts is to run amptools in parallel to measure bootstrapped errors. We run driveFit in parallel with different seeds for the dataset such that we get N bootstrapped datasets with N different sets of converged amplitudes. The bootstrapped amplitudes uses the amplitudes that have been converged to from the fit to the full/original dataset. The data is then written to text files which is then read in by PyROOT and the moments are plotted. These moments are old, need to update with Vincent's new formulations. THESE DO NOT HAVE NORMALIZATION INTEGRALS APPLIED SO THIS METHOD IS MOSTLY DEAD BUT THE FIT DRIVING SCRIPTS REMAIN MOSTLY USEFUL.


## The majority of parameters that needs to be changed is in the fit.cfg file
countEventsInBins.C is used to get the events in each mass bin so we can see if there is a deficit which might hinder PWA
fit_etapi_moments.cfg is copied by divideData.pl to all the bins
fullFit.py with args: iProcess binsPerProc processSpawned will run a single process starting on bin=iProcess and then stepping by processSpawned with binsPerProc iterations. So 0 3 20 would run on bins 0, 20, 60
      The reason we do this striding of bins is because the nearby bins should be similar in size, we want to homegenize. This program takes the amplitude ranges defined in this cfg file and samples from them. 
fullFit_spawn.sh runs fullFit.py on all the bins. 
bootFit.py and bootFit_spawn.sh does basically the same thing as the fullFit versions of them but this time it selects various random numbers as a seed for the bootstrapping. The intial amplitudes are based off of the converged
	amplitudes from the fullFit program. It allows a percent deviation from these converged values and a base shift. These values should probably be chosen such that the amplitudes are allowed to vary abit such that we aren't kind of #       enforcing some kind of similarity but also not allowing them to change too much as to hit a new minimum. Not sure if this makes sense....
checkUniqueSeedsAllConverged.py is used to check if  bootFit_spawn.sh did the right thing and initialized the seeds properly for all the bins
grabFull.py and grabBoot.py will look inside fitDir and grab all the amplitudes into individual files, one for each seed. 
	must use moveSeedSensi.sh first to move all the seed files to proper location. We can actually remove the need for moveSeedSensi.sh to save time but it cleans things up sort of 
plotAmps.py can now run through all these files and plot the amplitudes using pyroot
### ------------------------- TYPICAL USUAGE ----------------------------
DATA SETUP - we must first the data. DSelector_thrown outputs the data into 3 files with varying t. We still have to pass this through split_gen.C using run_split_thrown.sh.
           - WE MUST USE THE FLAT TREES FOR THE RECO AND DATA. HAVE TO CHANGE THE NAMES OF THE FLAT ROOT FILE IN THE DSELECTOR!!! THEN WE PASS IT THROUGH run_split.sh to get it correct form
	        1. root -l -b -q runDSelectorThrown_7_17_14.C
      	2. mv thrownNotAmptoolsReady_pi0eta_* amptools/newAnalysis
      	3. CHANGE NAME OF degAngle in DSelector_ver20.C to data_
      	4. root -l -b -q runDSelector_7_17_14.C
      	5. mv data_treeFlat_DSelector_pi0eta.root amptools/newAnalysis
      	6. CHANGE NAME OF degAngle in DSelector_ver20.C to reco_
      	7. root -l -b -q runDSelector_7_17_14.C
      	8. mv reco_treeFlat_DSelector_pi0eta.root amptools/newAnalysis
      	9. cd amptools/newAnalysis
      	10. ./run_split.sh after changing appropriate names in split_selected_in_t.C/h
      	11. ./run_split_thrown.sh after changing appropriate names in split_gen.C/h
1. ./divideData.pl
2. ./fullFit_spawn.sh
3. ./bootFit_spawn.sh
4. python grabFull.py
5. ./moveSeedSensi.sh
6. python checkUniqueSeedsAllConverged.py
7. python grabBoot.py 
8. python plotAmps.py
9. project_moments divideRoot NONE -full moments.root
