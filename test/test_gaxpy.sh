#!/bin/sh

echo 1..1

#########################################
INPUT_MATRIX=$srcdir/Matrix/m1.sms
OUT=$srcdir/Output/gaxpy.1
EXPECTED=$srcdir/Expected/gaxpy.1
./gaxpy < $INPUT_MATRIX > $OUT

if diff -w $OUT $EXPECTED >/dev/null ; then
  echo 'ok 1 - matrix-vector product'
else
  echo 'not ok 1 - matrix-vector product'
fi
