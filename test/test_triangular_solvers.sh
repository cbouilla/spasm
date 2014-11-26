#!/bin/sh

echo 1..4

#########################################
INPUT_MATRIX=$srcdir/Matrix/l1
./dense_solve 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/u1
./dense_solve 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/l1
./spsolve 3 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/u1
./spsolve 4 < $INPUT_MATRIX
