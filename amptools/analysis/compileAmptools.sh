cwd=$PWD
cd ~/gluex_top/AmpTools-0.9.4/AmpTools
make
cd /d/home/ln16/gluex_top/halld_sim-4.5.0^rec1701v03/src
scons -j20 install
cd $cwd
