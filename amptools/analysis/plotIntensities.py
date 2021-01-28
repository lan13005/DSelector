import pandas as pd
import math
import sys
import matplotlib
matplotlib.use('agg') # TKinter might not be installed so use agg
import matplotlib.pyplot as plt

label="GlueX acceptance"

def printHelp():
    print("Usage:\n--------")
    print("python plotIntensities.py <ampOrMom>")
    print("Will look for the csv output from diagnoseBS.py")
    print("<ampOrMom> is a flag to make the plots for moments with bootstrapped errors --- True")
    print("    or amplitudes with bootstrapped errors --- False")
    print("-------------------------------------")


argc=len(sys.argv)
if argc!=2:
    printHelp()
    exit()

if sys.argv[1]:
    # For moments
    massCol="M"
    uncertColName="uncert."
    fileName="project_moments_polarized"
    colToSetYlim="H0_00"
else:
    # For amplitudes
    massCol="M"
    uncertColName="err"
    fileName="plot_etapi_delta"
    colToSetYlim="all"

data=pd.read_csv("bootstrapDiagnostics/"+fileName+"_BS.csv")
print("Columns in dataframe:")
print(data.columns)

# Remove some columns that have some things we dont want
# in this case we try and remove 3/4 from the names which correspond to L=3/4 moments
smallLcols=[tmp for tmp in data.columns if tmp.find("3")==-1 and tmp.find("4")==-1]
data=data[smallLcols]

## Determine which columns are the mass, intensity, intensity errors
errCols=[tmp for tmp in data.columns if tmp.find(uncertColName)!=-1]
intCols=[tmp for tmp in data.columns if tmp.find(uncertColName)==-1 and tmp.find(massCol)==-1]
intCols=[tmp for tmp in intCols if tmp!=colToSetYlim]
intCols.insert(0,colToSetYlim)


### SETUP THINGS FOR THE FIGURES
nrows=int(math.ceil(math.sqrt(len(intCols))))
fig,axes=plt.subplots(nrows,nrows,figsize=(16,12))
axes=axes.flatten()
fig.suptitle(label, y=1.025,size=20)

color="indianred"
size=10
legend_loc=(1.6,)

for i in range(len(intCols)):
    data.plot(massCol,intCols[i],yerr=errCols[i],kind="scatter",ax=axes[i],c=color,s=size)
    if i==0:
        limits=axes[i].get_ylim()
    else:
        axes[i].set_ylim(limits)
    axes[i].set_ylabel("Intensity")
    axes[i].set_title(intCols[i])
    
plt.tight_layout()
fig.savefig("bootstrapDiagnostics/"+fileName+".png")
