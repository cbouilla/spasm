#!/bin/sh

echo 1..2

#########################################
INPUT_MATRIX=$srcdir/Matrix/small
./mat_perm 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/upper_trapeze
./mat_perm 2 < $INPUT_MATRIX
