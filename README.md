SpaSM (Sparse direct Solver Modulo _p_)
=======================================

SpaSM is a software library devoted to sparse gaussian elimination modulo a small prime _p_. 
It is available under the General Public License Version 3 or later (GPLv3+).

The algorithms used in SpaSM are described in [CASC'16](http://www-almasty.lip6.fr/~bouillaguet/pub/CASC16.pdf) and [PASCO'17](http://www-almasty.lip6.fr/~bouillaguet/pub/PASCO17.pdf).


Features
--------

The core of the library is an implementation of the GPLU algorithm, heavily inspired by 
[Tim Davis](http://faculty.cse.tamu.edu/davis/)'s [CSparse](http://faculty.cse.tamu.edu/davis/publications_files/CSparse.zip), and 
adapted to the context of exact computation. On top of this, we designed new strategies to search for structural pivots. 
This allows several kind of useful operations on sparse matrices:
  * LU and PLUQ factorization
  * Rank computation
  * Solution of linear systems
  * Kernel basis
  * Permutation to block triangular form
  * Reduced Row-Echelon Form

Finally, the library does I/O of matrices in [SMS format](http://hpac.imag.fr/), which makes it 
somewhat compatible with [LinBox](http://linalg.org/).

It is also capable of reading the [MatrixMarket format](https://math.nist.gov/MatrixMarket/), which is arguably better.

A set of demonstration programs is provided (see the `tools/` folder).


Installation
------------

In brief:
```make```

This requires [cmake](https://cmake.org). The executables can be found in `build/tools`.

SpaSM does not rely on any third-party software, but is capable of using:
  * [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) to find row separators.
  * [FFLAS-FFPACK](https://github.com/linbox-team/fflas-ffpack) for dense rank computation.
  * [LinBox](https://github.com/linbox-team/linbox) for other rank algorithms.
  * [Lemon](https://lemon.cs.elte.hu/trac/lemon) to find maximum matchings on non-bipartite graphs.

SpaSM uses OpenMP to exploit multicore machines.

The most commonly used option include:
- `--with-metis=<path>` : build the METIS interface
- `--with-fflas-ffpack=<path>` : enable the tools relying on dense rank computation
- `--with-linbox=<path>` : build the linbox wrappers (for comparison purpose)
- `--with-lemon=<path>` : build the lemon matching tool

Demonstration scripts
---------------------

All SpaSM demonstration scripts read a matrix in [SMS format](http://hpac.imag.fr/) on the standard input.

For instance, these commands (run inside the `build/tools/` folder) will compute the rank of several large matrices in a few seconds:
```
curl http://hpac.imag.fr/Matrices/Margulies/kneser_10_4_1.sms.gz | gunzip - | ./rank
curl http://hpac.imag.fr/Matrices/Homology/mk13.b5.135135x270270.sms.gz | gunzip - | ./rank
curl http://hpac.imag.fr/Matrices/G5/IG5-17.sms.gz | gunzip - | ./rank
```

It would be necessary to disable greedy pivot search for this one:
```
curl http://hpac.imag.fr/Matrices/Mgn/M0,6.data/M0,6-D9.sms.gz | gunzip - | ./rank
```

When matrices have many empty rows/columns, they can/have to be removed with the `stack` utility:
```
curl http://hpac.imag.fr/Matrices/Relat/relat8.sms.gz | gunzip - | ./stack | ./rank
curl http://hpac.imag.fr/Matrices/Relat/relat9.sms.gz | gunzip - | ./stack | ./rank
```

Finding good pivots is crucial for the performance of any kind of sparse elimination procedure. The pivot-finding code is still a bit naïve. Sometimes it will find much more pivots, much faster, if the matrices are flipped around a vertical axis with the `vertical_swap` utility:
```
curl http://hpac.imag.fr/Matrices/GL7d/GL7d14.sms.gz | gunzip - | ./vertical_swap | ./rank --sparse-threshold 0.01
...
curl http://hpac.imag.fr/Matrices/GL7d/GL7d22.sms.gz | gunzip - | ./vertical_swap | ./rank --sparse-threshold 0.01
```

Citing SpaSM
------------

If by any luck your research depends on the SpaSM library, please consider citing the project as

```
@manual{spasm,
title = {{SpaSM}: a Sparse direct Solver Modulo $p$},
author = {The SpaSM group},
edition = {v1.3},
year = {2023},
note = {\url{http://github.com/cbouilla/spasm}}
}
```

Contact and discussion
----------------------

Please email <charles.bouillaguet@lip6.fr> for any questions.