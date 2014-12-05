#!/bin/sh

echo 1..3

#########################################
INPUT_MATRIX=$srcdir/Matrix/medium
./sort 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/trefethen.500
./sort 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/trefethen.2000
./sort 3 < $INPUT_MATRIX
