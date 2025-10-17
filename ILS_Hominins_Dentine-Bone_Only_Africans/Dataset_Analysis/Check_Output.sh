#!/bin/bash

for filename in ./Workspace/2_DATASETS/*; do

	shopt -s lastpipe
	ls $filename | wc -l | read NOF
	MIN=3


	#if (( $NOF >= $MIN )); then

	#echo "$filename"
	
	#fi
	
	if test -f ./Workspace/2_DATASETS/$filename/Alignment_Done; then
	
	echo "$filename"
  	echo "File exists."

	fi


done
