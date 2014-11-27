#!/bin/sh

echo 1..2

#########################################
INPUT_MATRIX=$srcdir/Matrix/small
./lu 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/medium
./lu 2 < $INPUT_MATRIX
