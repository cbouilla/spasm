#!/bin/sh

echo 1..1

#########################################
INPUT_MATRIX=$srcdir/Matrix/small
OUT=$srcdir/Output/lu.1
EXPECTED=$srcdir/Expected/lu.1

./lu < $INPUT_MATRIX >

if diff -w $OUT $EXPECTED >/dev/null ; then
  echo 'ok 1 - LU (square, invertible)'
else
  echo 'not ok 1 - LU (square, invertible)'
fi
