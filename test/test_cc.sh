#!/bin/sh

echo 1..2

#########################################
INPUT_MATRIX=$srcdir/Matrix/cc.sms
./cc 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/BIOMD0000000525.int.mpl.sms
./cc 2 < $INPUT_MATRIX
