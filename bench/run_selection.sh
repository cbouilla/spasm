#!/bin/bash
BASE_PATH=/usr/share/hpac.imag.fr/Matrices
while read line
do
	printf "%s --- " $line
    zcat $BASE_PATH/$line | ./rank 2> /dev/null
done < ../result/matrices_linbox_gauss_is_faster_than_spasm.txt
