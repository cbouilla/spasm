#!/bin/sh

echo 1..5

#########################################
INPUT_MATRIX=$srcdir/Matrix/small
./pluq 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/medium
./pluq 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/singular
./pluq 3 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/rectangular_h
./pluq 4 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/rectangular_l
./pluq 5 < $INPUT_MATRIX
