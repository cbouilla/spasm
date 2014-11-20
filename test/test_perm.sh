#!/bin/sh

echo 1..2

#########################################
INPUT_MATRIX=$srcdir/Matrix/t1
OUT=$srcdir/Output/perm.1
EXPECTED=$srcdir/Expected/perm.1
./mat_perm 1 < $INPUT_MATRIX > $OUT

if diff -w $OUT $EXPECTED >/dev/null ; then
  echo 'ok 1 - permuting the rows of a matrix'
else
  echo 'not ok 1 - permuting the rows of a matrix'
fi

#########################################
INPUT_MATRIX=$srcdir/Matrix/t1
OUT=$srcdir/Output/perm.2
EXPECTED=$srcdir/Expected/perm.2
./mat_perm 2 < $INPUT_MATRIX > $OUT

if diff -w $OUT $EXPECTED >/dev/null ; then
  echo 'ok 2 - permuting the rows and columns of a matrix'
else
  echo 'not ok 2 - permuting the rows and columns of a matrix'
fi
