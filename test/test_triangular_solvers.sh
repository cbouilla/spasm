#!/bin/sh

echo 1..6

#########################################
INPUT_MATRIX=$srcdir/Matrix/l1.sms
./dense_lsolve 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/lower_trapeze.sms
./dense_lsolve 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/u1.sms
./dense_usolve 3 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/upper_trapeze.sms
./dense_usolve 4 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/u1.sms
./sparse_usolve 5 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/upper_trapeze.sms
./sparse_usolve 6 < $INPUT_MATRIX
