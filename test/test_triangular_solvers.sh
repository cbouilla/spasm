#!/bin/sh

echo 1..6

#########################################
INPUT_MATRIX=$srcdir/Matrix/l1
./dense_lsolve 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/lower_trapeze
./dense_lsolve 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/u1
./dense_usolve 3 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/upper_trapeze
./dense_usolve 4 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/u1
./sparse_solve 5 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/upper_trapeze
./sparse_solve 6 < $INPUT_MATRIX
