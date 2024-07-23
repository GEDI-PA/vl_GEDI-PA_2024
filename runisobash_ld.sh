#!/bin/bash

mkdir input
mkdir output

j=0

while getopts l:w:n: flag 
do 
	case "${flag}" in 
		l) clist=${1};;
		w) wk=${2};;
		n) N=${3};;
	esac
done 



for country in $clist
do 
	Rscript global_process_part1_2024_MAAP_step123.R $country $wk & 
	((++j==N)) && { j=0; wait; }	
done
