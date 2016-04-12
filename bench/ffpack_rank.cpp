/* Copyright (c) FFLAS-FFPACK
* Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
* ========LICENCE========
* This file is part of the library FFLAS-FFPACK.
*
* FFLAS-FFPACK is free software: you can redistribute it and/or modify
* it under the terms of the  GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
* ========LICENCE========
*/
#include <fflas-ffpack/fflas-ffpack-config.h>
#include <givaro/modular-balanced.h>
#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/utils/timer.h>
#include <fflas-ffpack/utils/Matio.h>
#include <fflas-ffpack/utils/args-parser.h>
#include <fflas-ffpack/ffpack/ffpack.h>


#include <iostream>

extern "C" {
#include "spasm.h"
}

using namespace FFLAS;

int main(int argc, char** argv) {
	typedef Givaro::Modular<int> Ring;
	int prime, i, j, k, n, m, nz, rank, *Ti, *Tj, *Tx;
	
	prime = 42013;
	Ring F(prime);
	spasm_triplet *T;
	Ring::Element * A;

	T = spasm_load_sms(stdin, prime);
	A = fflas_new(F, T->n, T->m);

	m = T->m;
    n = T->n;
    Ti = T->i;
    Tj = T->j;
    Tx = T->x;
    nz = T->nz;

    for(k=0; k<nz; k++) {
    	F.init(*(A + Ti[k]*n + Tj[k]), Tx[k]);
    }

    rank = FFPACK::Rank(F, n, m, A, n);
	std::cout << rank << std::endl;

	fflas_delete(A);
    return 0;
}

