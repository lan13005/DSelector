sed -i "s@nameAffix@tGT05LT1@g" split_gen.C
sed -i "s@nameAffix@tGT05LT1@g" split_gen.h

root -l -b -q call_split_thrown.C

sed -i "s@tGT05LT1@nameAffix@g" split_gen.C
sed -i "s@tGT05LT1@nameAffix@g" split_gen.h

sed -i "s@nameAffix@tLT1@g" split_gen.C
sed -i "s@nameAffix@tLT1@g" split_gen.h

root -l -b -q call_split_thrown.C

sed -i "s@tLT1@nameAffix@g" split_gen.C
sed -i "s@tLT1@nameAffix@g" split_gen.h

sed -i "s@nameAffix@tLT06@g" split_gen.C
sed -i "s@nameAffix@tLT06@g" split_gen.h

root -l -b -q call_split_thrown.C

sed -i "s@tLT06@nameAffix@g" split_gen.C
sed -i "s@tLT06@nameAffix@g" split_gen.h
