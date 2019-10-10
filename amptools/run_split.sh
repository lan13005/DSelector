echo "Using tag: $1"
#if [ "$1" != "reco" ] && [ "$1" != "data" ]; then
#  echo "not a valid argumenet... either reco or data!"
#  exit 0
#fi

sed -i "s@tag@$1@g" split_selected_in_t.C
sed -i "s@tag@$1@g" split_selected_in_t.h
root -l -b -q call_split_selected.C
sed -i "s@$1@tag@g" split_selected_in_t.C
sed -i "s@$1@tag@g" split_selected_in_t.h

