#!/bin/sh

echo 1..17

#########################################
INPUT_MATRIX=$srcdir/Matrix/small.sms
./lu 1 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/medium.sms
./lu 2 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/singular.sms
./lu 3 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/rectangular_h.sms
./lu 4 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/rectangular_l.sms
./lu 5 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/small.sms
./permuted_lu 6 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/medium.sms
./permuted_lu 7 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/singular.sms
./permuted_lu 8 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/rectangular_h.sms
./permuted_lu 9 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/rectangular_l.sms
./permuted_lu 10 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/small.sms
./lu_find_pivot 11 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/medium.sms
./lu_find_pivot 12 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/singular.sms
./lu_find_pivot 13 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/rectangular_h.sms
./lu_find_pivot 14 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/rectangular_l.sms
./lu_find_pivot 15 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/BIOMD0000000525.int.mpl.sms
./super_find_pivot 16 < $INPUT_MATRIX

#########################################
INPUT_MATRIX=$srcdir/Matrix/BIOMD0000000424.int.mpl.sms
./super_find_pivot 17 < $INPUT_MATRIX
