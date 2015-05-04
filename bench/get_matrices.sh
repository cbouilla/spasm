#!/bin/bash
BASE_URL=http://hpac.imag.fr/Matrices
while read line
do
    wget -x -nH $BASE_URL/$line
#    gunzip Matrices/$line.gz
done < matrix_list.txt
