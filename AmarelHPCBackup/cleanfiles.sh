#!/bin/bash

kept=0 # initialize a counter for the number of files kept

for i in {1..100}; do
	oldfn="GammamnpL100gam08D6NC512W01maxNR1Sam${i}.jld2"
	oldfp="./${oldfn}"
	if [ -f "$oldfp" ]; then
		size=$(du -k "$oldfp" | cut -f1)
		if [ $size -ge $((4*1024*1024)) ]; then
			kept=$((kept+1))
			newfn="GammamnpL100gam08D6NC512W01maxNR1Sam${kept}.jld2"
			newfp="./${newfn}"
			mv "$oldfp" "$newfp"
		else
			rm "$oldfp"
		fi
	fi
done
