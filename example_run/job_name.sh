#!/bin/sh

temp=("310" "350" "400" "450" "500" "550" "600" "650" "700" "750")
#hbond=("0.05" "0.1" "0.2" "0.3" "0.4" "0.5" "0.6" "0.7" "0.8" "0.9" "1.0")
hbond=("0.00")
total=$((${#temp[@]} * ${#hbond[@]} - 1))
echo $total
#for i in $(seq $1 $1)
for i in $(seq $1 $1)
do
	hbond_index=$(($i / ${#temp[@]}))
	temp_index=$(($i % ${#temp[@]}))
	echo ${hbond[$hbond_index]} ${temp[$temp_index]}
done
