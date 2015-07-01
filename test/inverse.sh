#!/bin/sh

echo 1..1

#########################################
INPUT_MATRIX_L=$srcdir/Matrix/l1.sms
INPUT_MATRIX_M=$srcdir/Matrix/m1.sms

./inverse_and_prod 1 $INPUT_MATRIX_L $INPUT_MATRIX_M

