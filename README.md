SpaSM (Sparse direct Solver Modulo _p_)
=======================================

SpaSM is a software library devoted to sparse gaussian elimination modulo a small prime _p_. 
It is available under the General Public License Version 2 or later (GPLv2+).

The algorithms used in SpaSM are described in [this paper](http://cristal.univ-lille.fr/~bouillag/pub/CASC16.pdf).

Features
--------

The core of the library is an implementation of the GPLU algorithm, heavily inspired by CSparse, and 
adapted to the context of exact computation. On top of this, we designed a hybrid left-and-right looking algorithm. 
This allows several kind of useful operations on sparse matrices:
  * LU and PLUQ factorization
  * Rank computation
  * Solution of linear systems
  * Kernel basis
  * Permutation to block triangular form

Finally, the library does I/O of matrices in SMS format, which makes it somewhat compatible with LinBox.

A set of demonstration programs is provided (see the `bench` folder).

Installation
------------

In brief:
```./configure <options> && make && make check```

If you do not have the `configure` script, try:
```autoreconf -i```

SpaSM does not rely on any third-party software, but is capable of using:
  * (METIS)[http://glaros.dtc.umn.edu/gkhome/metis/metis/overview] to find row separators.
  * (FFLAS-FFPACK)[https://github.com/linbox-team/fflas-ffpack] for dense rank computation.
  * (LinBox)[https://github.com/linbox-team/linbox] for other rank algorithms.

The most commonly used option include:
- `--with-metis=<path>` : build the METIS interface
- `--with-fflas-ffpack=<path>` : enable the tools relying on dense rank computation
- `--with-linbox=<path>` : build the linbox wrappers (for comparison purpose)

Demonstration scripts
---------------------

All SpaSM demonstration scripts read a matrix in [SMS format](http://hpac.imag.fr/) on the standard input.

For instance, these commands (run inside the `bench/` folder) will compute the rank of several large matrices in a few seconds:
```
curl http://hpac.imag.fr/Matrices/Margulies/kneser_10_4_1.sms.gz | gunzip - | ./rank_hybrid
curl http://hpac.imag.fr/Matrices/Homology/mk13.b5.135135x270270.sms.gz | gunzip - | ./transpose | ./rank_hybrid 2
curl http://hpac.imag.fr/Matrices/Mgn/M0,6.data/M0,6-D9.sms.gz | gunzip - | ./rank_hybrid
```

When matrices have a special shape (very "tall"), a different algorithm can be used:
```
curl http://hpac.imag.fr/Matrices/Relat/relat8.sms.gz | gunzip - | ./ffpack_narrow 
curl http://hpac.imag.fr/Matrices/Relat/relat9.sms.gz | gunzip - | ./ffpack_narrow 
```

Citing SpaSM
------------

If by any luck your research depends on the SpaSM library, please consider citing the project as

```
@manual{spasm,
title = {{SpaSM}: a Sparse direct Solver Modulo $p$},
author = {The SpaSM group},
edition = {v1.0},
year = {2016},
note = {\url{http://github.com/cbouilla/spasm}}
}
```

Contact and discussion
----------------------

