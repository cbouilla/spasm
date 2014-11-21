#!/bin/sh

echo 1..2

#########################################
INPUT_MATRIX=$srcdir/Matrix/g1
OUT=$srcdir/Output/reach.1
EXPECTED=$srcdir/Expected/reach.1
./reach 0 < $INPUT_MATRIX > $OUT

if diff -w $OUT $EXPECTED >/dev/null ; then
  echo 'ok 1 - DFS'
else
  echo 'not ok 1 - DFS'
fi

#########################################
OUT=$srcdir/Output/reach.2
EXPECTED=$srcdir/Expected/reach.2
./reach 1 < $INPUT_MATRIX > $OUT

if diff -w $OUT $EXPECTED >/dev/null ; then
  echo 'ok 2 - DFS'
else
  echo 'not ok 2 - DFS'
fi
