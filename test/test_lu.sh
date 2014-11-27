#!/bin/sh

echo 1..1

#########################################
INPUT_MATRIX=$srcdir/Matrix/small

./lu 1 < $INPUT_MATRIX
