#!/bin/bash
BASE_URL=http://hpac.imag.fr/gbla
while read line
do
    wget $BASE_URL/$line -q -O - | gunzip |./../test/gbla_in_new | ./demo_cheap > ../test/Matrix/gblatest.sms
#    gunzip Matrices/$line.gz
done < gbla_list.txt
