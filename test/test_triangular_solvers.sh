#!/bin/sh

echo 1..4

#########################################
INPUT_MATRIX=$srcdir/Matrix/l1
./lsolve 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/u1
./lsolve 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/l1
./spsolve 3 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/u1
./spsolve 4 < $INPUT_MATRIX
