#!/bin/sh

echo 1..3

#########################################
INPUT_MATRIX=$srcdir/Matrix/l1
./dense_solve 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/u1
./dense_solve 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/l1
./sparse_solve 3 < $INPUT_MATRIX
