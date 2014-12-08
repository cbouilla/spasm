#!/bin/sh

echo 1..2

INPUT_MATRIX=$srcdir/Matrix/t1.sms

# test 1 : reading/writing a matrix in triplet form
OUT=$srcdir/Output/io.1
EXPECTED=$srcdir/Expected/io.1

./io 1 < $INPUT_MATRIX > $OUT

if diff -w $OUT $EXPECTED >/dev/null ; then
  echo 'ok 1 - reading and writing a matrix'
else
  echo 'not ok 1 - reading and writing a matrix'
fi

# test 2 : reading a matrix in triplet form, writing it in CSR
OUT=$srcdir/Output/io.2
EXPECTED=$srcdir/Expected/io.2

./io 2 < $INPUT_MATRIX > $OUT

if diff -w $OUT $EXPECTED >/dev/null ; then
  echo 'ok 2 - reading a matrix, converting it to CSR and writing it'
else
  echo 'not ok 2 -  reading a matrix, converting it to CSR and writing it'
fi
