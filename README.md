SpaSM (Sparse direct Solver Modulo _p_)
=======================================

SpaSM is a software library devoted to sparse gaussian elimination modulo a small prime _p_. 
It is available under the General Public License Version 3 or later (GPLv3+).

This is "research-quality" software. While we try to make sure that it actually works, we cannot make any promises.

Features
--------

The core of the library is a multithreaded function that computes the row echelon form of a sparse matrix modulo a word-sized prime. 

This enables several of useful computations on sparse matrices modulo _p_:
  * Basis of the row space
  * Basis of the kernel
  * Row Echelon form and Reduced Row Echelon form
  * Rank

In addition, SpaSM is capable of computing a full PLUQ factorization, but this is slower than simple echelonization. This enables:
  * Solution of linear systems
  * Extraction of a square submatrix of maximal rank
  * Production of rank certificates

SpaSM works with all odd 32-bit prime moduli.  While it would be possible to work with _p = 2_, this would require non-trivial tweaks to the code. 

Finally, SpaSM contains code to compute the Dulmage-Mendelson decomposition (permutation of a matrice to block triangular form) and several other useful functions.

The following algorithms algorithms are used in SpaSM:
  * [Gilbert-Peierls sparse triangular solving with sparse right-hand side](https://doi.org/10.1137/0909058)
  * [Faugère-Lachartre pivot selection](http://www-almasty.lip6.fr/~bouillaguet/pub/CASC16.pdf)
  * [Improved greedy pivot selection](http://www-almasty.lip6.fr/~bouillaguet/pub/PASCO17.pdf)
  * [Dense linear algebra mod _p_](https://hal.science/hal-00018223/)
  * [Linear-size rank certificates](https://prism.ucalgary.ca/server/api/core/bitstreams/b00bb76d-12bf-41c2-9fed-88cd774d3b29/content)

Initial versions of SpaSM were heavily influence by
[Tim Davis](http://faculty.cse.tamu.edu/davis/)'s [CSparse](http://faculty.cse.tamu.edu/davis/publications_files/CSparse.zip). 

Spasm does I/O of matrices in either the :
  * [SMS format](http://hpac.imag.fr/) --- which makes it somewhat compatible with [LinBox](http://linalg.org/).  
  * [MatrixMarket format](https://math.nist.gov/MatrixMarket/), with kind `matrix coordinate integer general`.
(The actual input format is autodetected).

A set of demonstration programs is provided (see the `tools/` folder). They can be used to compute the rank, the RREF, a kernel basis, or a Dumlage-Mendelson decomposition of a sparse matrix.


Installation
------------

In brief:
```make```

This requires [cmake](https://cmake.org) and [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/). The executables can be found in `build/tools`.

SpaSM relies on two third-party libraries, which are required at compile-time:
  * [Givaro](https://github.com/linbox-team/givaro)
  * [FFLAS-FFPACK](https://github.com/linbox-team/fflas-ffpack)
  
Under Debian linux or ubuntu, installing the `fflas-ffpack` package is sufficient to compile.

SpaSM uses OpenMP to exploit multicore machines.

Demonstration scripts
---------------------

All SpaSM demonstration scripts read a matrix in [SMS format](http://hpac.imag.fr/) on the standard input.

For instance, these commands (run inside the `build/tools/` folder) will compute the rank of several large matrices in a few seconds:
```
curl https://hpac.imag.fr/Matrices/Margulies/kneser_10_4_1.sms.gz | gunzip - | ./rank
curl https://hpac.imag.fr/Matrices/Homology/mk13.b5.135135x270270.sms.gz | gunzip - | ./rank
curl https://hpac.imag.fr/Matrices/G5/IG5-17.sms.gz | gunzip - | ./rank
```

It would be necessary to disable greedy pivot search for this one:
```
curl https://hpac.imag.fr/Matrices/Mgn/M0,6.data/M0,6-D9.sms.gz | gunzip - | ./rank
```

When matrices have many empty rows/columns, they can/have to be removed with the `stack` utility:
```
curl https://hpac.imag.fr/Matrices/Relat/relat8.sms.gz | gunzip - | ./stack | ./rank
curl https://hpac.imag.fr/Matrices/Relat/relat9.sms.gz | gunzip - | ./stack | ./rank
```

Finding good pivots is crucial for the performance of any kind of sparse elimination procedure. The pivot-finding code is still a bit naïve. Sometimes it will find much more pivots, much faster, if the matrices are flipped around a vertical axis with the `vertical_swap` utility:
```
curl https://hpac.imag.fr/Matrices/GL7d/GL7d14.sms.gz | gunzip - | ./vertical_swap | ./rank --sparse-threshold 0.01
...
curl https://hpac.imag.fr/Matrices/GL7d/GL7d22.sms.gz | gunzip - | ./vertical_swap | ./rank --sparse-threshold 0.01
```

Dealing with large matrices
---------------------------

Sparse Gaussian elimination is more an art than a science.  In some cases, it will inevitably fail (the matrix will fill and the process will grind to a halt). 

However, with some expertise, it may be possible to deal with potentially larger problems than what the auto-pilot is capable of. Don't hesitate to get in touch.

Citing SpaSM
------------

If by any luck your research depends on the SpaSM library, please consider citing the project as

```
@manual{spasm,
title = {{SpaSM}: a Sparse direct Solver Modulo $p$},
author = {Charles Bouillaguet},
edition = {v1.3},
year = {2023},
note = {\url{http://github.com/cbouilla/spasm}}
}
```

Contact and discussion
----------------------

Please email <charles.bouillaguet@lip6.fr> for any questions.