#!/bin/sh

echo 1..1

#########################################
INPUT_MATRIX=$srcdir/Matrix/cc.sms
./cc 1 < $INPUT_MATRIX
