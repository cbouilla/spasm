#!/bin/sh

echo 1..5

#########################################
INPUT_MATRIX=$srcdir/Matrix/small
./solve 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/medium
./solve 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/singular
./solve 3 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/upper_trapeze
./solve 4 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/lower_trapeze
./solve 5 < $INPUT_MATRIX
