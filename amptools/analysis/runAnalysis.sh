rm -r divideRoot
./divideData.pl
./fullFit_spawn.sh
./bootFit_spawn.sh
./moveSeedSensi.sh
python grabFull.py
python grabBoot.py
python plot.py -b