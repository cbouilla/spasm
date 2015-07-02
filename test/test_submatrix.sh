#!/bin/sh

echo 1..2


# test 1 : reading/writing a matrix in triplet form
INPUT_MATRIX=$srcdir/Matrix/singular.sms
OUT=$srcdir/Output/submatrix.1
EXPECTED=$srcdir/Expected/submatrix.1

./submatrix 1 < $INPUT_MATRIX > $OUT

if diff -w $OUT $EXPECTED >/dev/null ; then
    echo 'ok 1 - range submatrix'
else
  echo 'not ok 1 - range submatrix'
fi


# test 2 : with sorted rows
INPUT_MATRIX=$srcdir/Matrix/m1.sms
OUT=$srcdir/Output/submatrix.2
EXPECTED=$srcdir/Expected/submatrix.2

./sorted_submatrix 2 < $INPUT_MATRIX > $OUT

if diff -w $OUT $EXPECTED >/dev/null ; then
    echo 'ok 2 - range submatrix (sorted)'
else
  echo 'not ok 2 - range submatrix (sorted)'
fi
