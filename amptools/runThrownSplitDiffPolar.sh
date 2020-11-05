# sed/// as / as a delimeter but can take anything incase our replacement has that delimeter in it.

root -l -b -q call_split_thrown.C

sed -i 's@deg000_gen@deg045_gen@g' split_gen.C
sed -i 's@deg000_gen@deg045_gen@g' split_gen.h
root -l -q -b call_split_thrown.C 

sed -i 's@deg045_gen@deg090_gen@g' split_gen.C
sed -i 's@deg045_gen@deg090_gen@g' split_gen.h
root -l -q -b call_split_thrown.C 

sed -i 's@deg090_gen@deg135_gen@g' split_gen.C
sed -i 's@deg090_gen@deg135_gen@g' split_gen.h
root -l -q -b call_split_thrown.C 

sed -i 's@deg135_gen@degAMO_gen@g' split_gen.C
sed -i 's@deg135_gen@degAMO_gen@g' split_gen.h
root -l -q -b call_split_thrown.C 

sed -i 's@degAMO_gen@deg000_gen@g' split_gen.C
sed -i 's@degAMO_gen@deg000_gen@g' split_gen.h


