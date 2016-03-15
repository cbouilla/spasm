#!/bin/sh

echo 1..3

INPUT_MATRIX=$srcdir/Matrix/m1.sms
./column_slices 1 < $INPUT_MATRIX

#####################################################"
INPUT_MATRIX=$srcdir/Matrix/medium.sms
./column_slices 2 < $INPUT_MATRIX

#####################################################"
INPUT_MATRIX=$srcdir/Matrix/rectangular_l.sms
./column_slices 3 < $INPUT_MATRIX

