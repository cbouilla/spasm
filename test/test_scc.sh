#!/bin/sh

echo 1..1

#########################################
INPUT_MATRIX=$srcdir/Matrix/scc.sms
./scc 1 < $INPUT_MATRIX
