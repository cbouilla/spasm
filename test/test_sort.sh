#!/bin/sh

echo 1..3

#########################################
INPUT_MATRIX=$srcdir/Matrix/medium.sms
./sort 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/trefethen_500.sms
./sort 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/trefethen_2000.sms
./sort 3 < $INPUT_MATRIX
