#!/usr/bin/perl

use Cwd;

$lowMass = 0.7; #is a shared lower cutoff for all 3 datas.
$highMass = 2; #2 is the upper cutoff of the thrown data and 3ish is the upper cutoff for the reco/data
#$nBins = 65; # not sure why old me chose 65 bins... that is alot of bins and we might not have enough statistics
$nBins=65; # 26 because it is kind of small and (2-0.7)/26 = 0.05 which is nice and round

$fitName = "EtaPi_fit";
$nameAffix = "amptools";

# put a limit on the number of data events to process
# gen MC and acc MC smaples are not limited
$maxEvts = 1E9;

# this directory can be adjusted if you want to do the fit elsewhere
# but it needs to be an explicit path
$workinDir = getcwd();

# these files must exist in the workin directory.  If you don't know how
# to generate them or don't have them, see the documentation in gen_3pi
# the Simulation area of the repository
#$dataFile = "$workinDir/trees/tree_pimpip_coh_2017_amptools.root";
$dataFile = "$workinDir/a0a2_sim/dat_pi0eta_$nameAffix\_w1.root";
$accMCFile = "$workinDir/a0a2_sim/acc_pi0eta_$nameAffix.root";
$genMCFile = "$workinDir/a0a2_sim/gen_pi0eta_$nameAffix.root";

print $dataFile;
print "\n";
print $accMCFile;
print "\n";
print $genMCFile;
print "\n";

# this file sould be used for partially polarized or unpolarized beam fits
#$cfgTempl = "$workinDir/threepi_unpol_TEMPLATE.cfg";

$cfgTempl = "$workinDir/fit_etapi_moments.cfg";


### things below here probably don't need to be modified

# this is where the goodies for the fit will end up
$fitDir = "$workinDir/$fitName/";
mkdir $fitDir unless -d $fitDir;

chdir $fitDir;

print "Changing into $fitDir\n";

# use the split_mass command line tool to divide up the
# data into bins of resonance mass
@dataParts = split /\//, $dataFile; #splitting the path into sections so we can pop out the last important file name
$dataTag = pop @dataParts;
$dataTag =~ s/\.root//;
print "datatag: $dataTag\n";
# here we will use the thrown MC data to get the best represention of how PWA calculates things
system( "split_mass $dataFile $dataTag $lowMass $highMass $nBins $maxEvts -T kin:kin" );
print "splitted data!\n"; 

@accMCParts = split /\//, $accMCFile;
$accMCTag = pop @accMCParts;
$accMCTag =~ s/\.root//;
print "datatag: $accMCTag\n";
system( "split_mass $accMCFile $accMCTag $lowMass $highMass $nBins -T kin:kin" );
print "splitted data!\n"; 

@genMCParts = split /\//, $genMCFile;
$genMCTag = pop @genMCParts;
$genMCTag =~ s/\.root//;
print "datatag: $genMCTag\n";
system( "split_mass $genMCFile $genMCTag $lowMass $highMass $nBins -T kin:kin" );
#system( "split_mass $genMCFile $genMCTag $lowMass $highMass $nBins -T $nameAffix:kin" );
print "splitted data!\n"; 

# make directories to perform the fits in
for( $i = 0; $i < $nBins; ++$i ){

  mkdir "bin_$i" unless -d "bin_$i";

  system( "mv *\_$i.root bin_$i" );

  chdir "bin_$i";

#we are essentially copying fit_etapi_moments.cfg and substituting some variables. CFGOUT is going to be a config file in all of our bins. CFGIN is fit_etapi_moments.cfg. Note how fit_etapi_moments.cfg has these place holders defined (DATAFILE,ACCMCFILE,GENMCFILE ... ). They will get replaced here to fit the bin directory. 
  open( CFGOUT, ">bin_$i.cfg" );
  open( CFGIN, $cfgTempl ); 

  while( <CFGIN> ){

    # for some reason we have to print these following statements or else the substition in the following lines wont work....
    #print("dataTag: $dataTag\_$i.root\n");
    #print("accMCTag: $accMCTag\_$i.root\n");
    #print("genMCTag: $genMCTag\_$i.root\n");
    s/DATAFILE/$dataTag\_$i.root/;
    s/ACCMCFILE/$accMCTag\_$i.root/;
    s/GENMCFILE/$genMCTag\_$i.root/;
    s/NIFILE/bin_$i.ni/;
    s/FITNAME/bin_$i/;

    print CFGOUT $_;
  }

  close CFGOUT;
  close CFGIN;
  
  system( "touch param_init.cfg" );

  chdir $fitDir;
}

