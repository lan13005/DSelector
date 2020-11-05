# sed/// as / as a delimeter but can take anything incase our replacement has that delimeter in it.

root -l -b -q call_split_selected.C

sed -i 's@deg000_data@deg045_data@g' split_selected_in_t.C
sed -i 's@deg000_data@deg045_data@g' split_selected_in_t.h
root -l -q -b call_split_selected.C 

sed -i 's@deg045_data@deg090_data@g' split_selected_in_t.C
sed -i 's@deg045_data@deg090_data@g' split_selected_in_t.h
root -l -q -b call_split_selected.C 

sed -i 's@deg090_data@deg135_data@g' split_selected_in_t.C
sed -i 's@deg090_data@deg135_data@g' split_selected_in_t.h
root -l -q -b call_split_selected.C 

sed -i 's@deg135_data@degAMO_data@g' split_selected_in_t.C
sed -i 's@deg135_data@degAMO_data@g' split_selected_in_t.h
root -l -q -b call_split_selected.C 

sed -i 's@degAMO_data@deg000_data@g' split_selected_in_t.C
sed -i 's@degAMO_data@deg000_data@g' split_selected_in_t.h


