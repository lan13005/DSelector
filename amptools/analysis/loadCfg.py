import os
import ast

def loadCfg():
	cfgfile="fit.cfg"
	cfgDict={}
	if os.path.isfile(cfgfile):
		with open(cfgfile,"r") as cfg:
			lines = cfg.readlines()
			for line in lines:
				if line.startswith("#") or len(line.strip())==0:
					continue
				variable=line.split("=")[0].strip()
				value=line.split("=")[1].strip()
				if len(value.split(","))>0:
					value=ast.literal_eval(value)
				#print(value)
				cfgDict[variable]=value
	else:
		raise ValueError("Config files does not exist!") 
	return cfgDict

#cfg = loadCfg()
#print(cfg["FIT_DIR"])
