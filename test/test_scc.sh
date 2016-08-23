#!/bin/sh

echo 1..3

#########################################
INPUT_MATRIX=$srcdir/Matrix/scc.sms
./scc 1 < $INPUT_MATRIX

INPUT_MATRIX=$srcdir/Matrix/scc2.sms
./scc 2 < $INPUT_MATRIX

INPUT_MATRIX=$srcdir/Matrix/mat364.sms
./scc 3 < $INPUT_MATRIX