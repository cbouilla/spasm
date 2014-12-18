#!/bin/sh

echo 1..2

#########################################
INPUT_MATRIX=$srcdir/Matrix/dm.sms
./dm 1 < $INPUT_MATRIX

INPUT_MATRIX=$srcdir/Matrix/dm2.sms
./dm 2 < $INPUT_MATRIX
