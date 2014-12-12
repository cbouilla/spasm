#!/bin/sh

echo 1..1

#########################################
INPUT_MATRIX=$srcdir/Matrix/dm.sms
./dm 1 < $INPUT_MATRIX
