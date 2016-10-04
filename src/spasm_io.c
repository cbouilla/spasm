#include <assert.h>
#include <math.h>
#include <err.h>

#include "spasm.h"
#include "mmio.h"


/*
 * load a matrix in SMS format from f (an opened file). 
 * set prime == -1 to avoid loading values.
 */
spasm_triplet *spasm_load_sms(FILE * f, int prime) {
	int i, j;
	spasm_GFp x;
	char type;
	
	assert(f != NULL);

	double start = spasm_wtime();
	if (fscanf(f, "%d %d %c\n", &i, &j, &type) != 3)
		errx(1, "[spasm_load_sms] bad SMS file (header)\n");

	if (prime != -1 && type != 'M')
		errx(1, "[spasm_load_sms] only ``Modular'' type supported\n");

	fprintf(stderr, "[IO] loading %d x %d SMS matrix modulo %d... ", i, j, prime);
	fflush(stderr);

	/* allocate result */
	spasm_triplet *T = spasm_triplet_alloc(i, j, 1, prime, prime != -1);

	while (fscanf(f, "%d %d %d\n", &i, &j, &x) == 3) {
		if (i == 0 && j == 0 && x == 0)
			break;
		assert(i != 0);
		assert(j != 0);
		spasm_add_entry(T, i - 1, j - 1, x);
	}

	char nnz[16];
	spasm_human_format(T->nz, nnz);
	fprintf(stderr, "%s NNZ [%.1fs]\n", nnz, spasm_wtime() - start);
	return T;
}

/*
 * Load a matrix in MatrixMarket sparse format.
 * Heavily inspired by the example program:
 *     http://math.nist.gov/MatrixMarket/mmio/c/example_read.c
 */
spasm_triplet *spasm_load_mm(FILE * f, int prime) {
	MM_typecode matcode;
	int n, m, nnz;

	double start = spasm_wtime();
	if (mm_read_banner(f, &matcode) != 0) 
		errx(1, "Could not process Matrix Market banner.\n");

	if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode) || !mm_is_general(matcode))
		errx(1, "Market Market type: [%s] not supported", mm_typecode_to_str(matcode));
	
	if (mm_read_mtx_crd_size(f, &n, &m, &nnz) !=0)
		errx(1, "Cannot read matrix size");

	fprintf(stderr, "[IO] loading %d x %d MTX [%s] modulo %d...", n, m, mm_typecode_to_str(matcode), prime);
	fflush(stderr);


	if (mm_is_pattern(matcode))
		prime = -1;

	spasm_triplet *T = spasm_triplet_alloc(n, m, nnz, prime, prime != -1);

	for (int i = 0; i < nnz; i++) {
		int u, v, w;
		double x, y;

		if (mm_is_pattern(matcode)) {
			fscanf(f, "%d %d\n", &u, &v);
			spasm_add_entry(T, u - 1, v - 1, 1);
		} else if (mm_is_integer(matcode)) {
			fscanf(f, "%d %d %d\n", &u, &v, &w);
			spasm_add_entry(T, u - 1, v - 1, w);
		} else if (mm_is_real(matcode)) {
			fscanf(f, "%d %d %lg\n", &u, &v, &x);
			spasm_add_entry(T, u - 1, v - 1, (int) (100000 * x));
		} else if (mm_is_complex(matcode)) {
			fscanf(f, "%d %d %lg %lg\n", &u, &v, &y, &x);
			spasm_add_entry(T, u - 1, v - 1, (int) (1000 * (y + 100 * x)));
		} else {
			errx(1, "Don't know how to read matrix");
		}
	}

	char s_nnz[16];
	spasm_human_format(T->nz, s_nnz);
	fprintf(stderr, "%s NNZ [%.1fs]\n", s_nnz, spasm_wtime() - start);
	return T;
}

/* load a matrix in old GBLA format */
spasm *spasm_load_gbla_old(FILE * f, int with_values) {
	int m, n, p;
	long long int nnz;
	spasm *M;

	/* get rows */
	if (fread(&n, sizeof(uint32_t), 1, f) != 1)
		err(1, "error reading m");

	/* get columns */
	if (fread(&m, sizeof(uint32_t), 1, f) != 1)
		err(1, "error reading n");

	if ((fread(&p, sizeof(uint32_t), 1, f) != 1))
		err(1, "error reading p");

	/* get number of nonzero elements */
	if (fread(&nnz, sizeof(uint64_t), 1, f) != 1)
		err(1, "error reading nnz");

	if (nnz >> 31ull)
		errx(2, "more than 2^31 NNZ. You are going to hit limitations of spasm... Sorry\n");

	M = spasm_csr_alloc(n, m, nnz, p, with_values);
	int * const Mp = M->p;
	int * const Mj = M->j;
	spasm_GFp * const Mx = M->x;

	/*
	 * dirty hack : we load the 16-bit coefficient in the area dedicated
	 * to storing the column indices (unused at this point). This avoids
	 * allocating extra memory
	 */
	uint16_t *buffer = (uint16_t *) M->j;
	size_t r = nnz;
	if (fread(buffer, sizeof(uint16_t), nnz, f) != r)
		err(1, "error while reading the coefficients");
	if (with_values)
		for (int i = 0; i < nnz; i++)
			Mx[i] = buffer[i];
	/* read the column indices */
	r = nnz;
	if (fread(Mj, sizeof(uint32_t), nnz, f) != r)
		err(1, "error while reading the column indices");
	/* read the row lengths */
	r = n;
	if (fread(&Mp[1], sizeof(uint32_t), n, f) != r)
		err(1, "error while reading the row lengths (read %zu items)", r);
	/* sum-prefix */
	Mp[0] = 0;
	for (int i = 1; i <= n; i++)
		Mp[i] += Mp[i - 1];
	return M;
}


#define VERMASK (1U<<31)

/*
 * load a matrix in new GBLA format. Only the pattern is loaded, not the
 * values.
 */
spasm *spasm_load_gbla_new(FILE * f) {
	unsigned int i, j, m, n;
	uint64_t nnz;
	uint16_t p;
	spasm *M;
	uint32_t b;

	/* get header */
	if (fread(&b, sizeof(uint32_t), 1, f) != 1)
		err(1, "error reading b");
	if ((b & VERMASK) != VERMASK)
		errx(1, "wrong format version");
	b = b ^ VERMASK;
	if (((b >> 1) & 3) != 1)
		errx(1, "field elements are not on 16 bits");
	/* get rows */
	if (fread(&n, sizeof(uint32_t), 1, f) != 1)
		err(1, "error reading m");
	/* get columns */
	if (fread(&m, sizeof(uint32_t), 1, f) != 1)
		err(1, "error reading n");
	if ((fread(&p, sizeof(uint16_t), 1, f) != 1))
		err(1, "error reading p");
	/* get number of nonzero elements */
	if (fread(&nnz, sizeof(uint64_t), 1, f) != 1)
		err(1, "error reading nnz");
	if (nnz >> 31ull)
		errx(2, "more than 2^31 NNZ. You are going to hit limitations of spasm... Sorry\n");
	M = spasm_csr_alloc(n, m, nnz, (int)p, SPASM_IGNORE_VALUES);
	int * const Mp = M->p;
	int * const Mj = M->j;

	/* read the row lengths */
	size_t r = n;
	if (fread(&Mp[1], sizeof(uint32_t), n, f) != r)
		err(1, "error while reading the row lengths (read %zu items)", r);
	/* sum-prefix */
	Mp[0] = 0;
	for (unsigned int i = 1; i <= n; i++)
		Mp[i] += Mp[i - 1];

	/*
	 * read (and ignore) the polmaps. Dirty hack: we send them to the
	 * area for the column indices
	 */
	r = n;
	if (fread(Mj, sizeof(uint32_t), n, f) != r)
		err(1, "error while reading the polmaps (read %zu items)", r);

	uint64_t k;
	/* get number of compressed columns element */
	if (fread(&k, sizeof(uint64_t), 1, f) != 1)
		err(1, "error reading k");

	/* read compressed columns indices */
	uint32_t *buffer = spasm_malloc(k * sizeof(uint32_t));
	r = k;
	if (fread(buffer, sizeof(uint32_t), k, f) != r)
		err(1, "error while reading the compressed column ids (read %zu items)", r);

	/* uncompress columns indices */
	i = 0;
	j = 0;
	const uint32_t MASK = 0x80000000;
	while (i < k) {
		uint32_t col = buffer[i++];

		if (col & MASK) {	/* single column */
			Mj[j++] = col ^ MASK;
		} else {
			uint32_t x = buffer[i++];
			for (p = 0; p < x; p++) {
				Mj[j++] = col + p;
			}
		}
	}
	assert(i == k);
	assert(j == nnz);
	free(buffer);

	return M;
}


/*
 * save a matrix in SMS format. TODO : change name to spasm_csr_save
 */
void spasm_save_csr(FILE * f, const spasm * A) {
	assert(f != NULL);

	int *Aj = A->j;
	int *Ap = A->p;
	spasm_GFp *Ax = A->x;
	int n = A->n;
	int m = A->m;
	int prime = A->prime;

	fprintf(f, "%d %d M\n", n, m);
	for (int i = 0; i < n; i++)
		for (int px = Ap[i]; px < Ap[i + 1]; px++) {
			spasm_GFp x = (Ax != NULL) ? Ax[px] : 1;
			x = (x > prime / 2) ? x - prime : x;
			fprintf(f, "%d %d %d\n", i + 1, Aj[px] + 1, x);
		}
	fprintf(f, "0 0 0\n");
}

/*
 * save a matrix in SMS format. TODO : change name to spasm_triplet_save
 */
void spasm_save_triplet(FILE * f, const spasm_triplet * A) {
	assert(f != NULL);

	int *Ai = A->i;
	int *Aj = A->j;
	spasm_GFp *Ax = A->x;
	int nz = A->nz;

	fprintf(f, "%d %d M\n", A->n, A->m);
	for (int px = 0; px < nz; px++)
		fprintf(f, "%d %d %d\n", Ai[px] + 1, Aj[px] + 1, (Ax != NULL) ? Ax[px] : 1);
	fprintf(f, "0 0 0\n");
}


/* Saves a PBM (bitmap) of specified dimensions of A */
void spasm_save_pbm(FILE * f, int x, int y, const spasm * A) {
	int i, j, k, n, m, t, p;
	int *Aj, *Ap, *w;

	assert(f != NULL);
	assert(A != NULL);

	Aj = A->j;
	Ap = A->p;
	n = A->n;
	m = A->m;
	x = spasm_min(x, m);
	y = spasm_min(y, n);

	w = spasm_malloc(x * y * sizeof(int));
	for (j = 0; j < x * y; j++) {
		w[j] = 0;
	}

	fprintf(f, "P1\n");
	fprintf(f, "%d %d\n", x, y);

	for (i = 0; i < n; i++) {
		k = ((long long int)i) * ((long long int)y) / ((long long int)n);
		int *ww = w + k * x;
		for (p = Ap[i]; p < Ap[i + 1]; p++) {
			t = ((long long int)Aj[p]) * ((long long int)x) / ((long long int)m);
			ww[t] = 1;
		}
	}

	/* print row */
	t = 0;
	for (k = 0; k < y; k++) {
		for (j = 0; j < x; j++) {
			fprintf(f, "%d ", w[k * x + j]);
			t++;
			if ((t & 40) == 0) {
				fprintf(f, "\n");
			}
		}
	}

	free(w);
}

/* Saves a PGM (graymap) of specified dimensions of A */
void spasm_save_pgm(FILE * f, int x, int y, const spasm * A) {
	int i, j, k, n, m, t, p;
	int *Aj, *Ap, *w;
	double max;

	assert(f != NULL);
	assert(A != NULL);

	Aj = A->j;
	Ap = A->p;
	n = A->n;
	m = A->m;
	x = spasm_min(x, m);
	y = spasm_min(y, n);

	w = spasm_malloc(x * y * sizeof(int));
	for (j = 0; j < x * y; j++) {
		w[j] = 0;
	}

	fprintf(f, "P2\n");
	fprintf(f, "%d %d\n", x, y);
	fprintf(f, "255\n");

	for (i = 0; i < n; i++) {
		k = ((long long int)i) * ((long long int)y) / ((long long int)n);
		int *ww = w + k * x;
		for (p = Ap[i]; p < Ap[i + 1]; p++) {
			t = ((long long int)Aj[p]) * ((long long int)x) / ((long long int)m);
			ww[t]++;
		}
	}

	/* scan for max */
	max = 0;
	for (j = 0; j < x * y; j++) {
		if (w[j] > max) {
			max = w[j];
		}
	}

	/* print row */
	t = 0;
	for (k = 0; k < y; k++) {
		for (j = 0; j < x; j++) {
			double intensity = 1.0 - exp(0.1 * log(w[k * x + j] / max));
			assert(0 <= intensity && intensity <= 1.0);
			fprintf(f, "%.0f ", 255.0 * intensity);
			t++;
			if ((t & 31) == 0) {
				fprintf(f, "\n");
			}
		}
	}

	free(w);
}

#if 0
static void render_block(FILE * f, int m, int *Ap, int *Aj, spasm_partition * CC, spasm_partition ** SCC, int *rr, int *cc,
	          int ri, int rj, int ci, int cj, int *colors, int *pixel) {
	int u, i, j, k, l, t, p, CC_n, CC_m, last_CC_row;

	t = 0;
	k = 0;			/* in which CC are we ? */
	l = 0;			/* in which SCC are we ? */
	assert(ci < cj);
	last_CC_row = CC->rr[CC->nr];

	for (i = rr[ri]; i < rr[rj]; i++) {

		spasm_vector_set(pixel, 0, m, 0xFFFFFF);	/* white */

		/* jump CC / SCC */
		while (i < last_CC_row && i >= CC->rr[k + 1]) {
			k++;
			l = 0;
		}
		while (i < last_CC_row && i >= SCC[k]->rr[l + 1]) {
			l++;
		}

		/* are we in a matched row ? */
		CC_n = CC->rr[k + 1] - CC->rr[k];
		CC_m = CC->cc[k + 1] - CC->cc[k];

		if (i < CC->rr[k + 1] - CC_m) {
			/* unmatched BG */
			for (j = CC->cc[k]; j < CC->cc[k + 1]; j++) {
				pixel[j] = colors[0];
			}
		} else {
			/* diagonal block of the SCC */
			for (j = SCC[k]->cc[l]; j < SCC[k]->cc[l + 1]; j++) {
				pixel[j] = colors[1];
			}

			/*
			 * stuff on the right of the diagonal block, inside
			 * the CC
			 */
			for (j = SCC[k]->cc[l + 1]; j < CC->cc[k + 1]; j++) {
				pixel[j] = colors[2];
			}

			for (u = cj; u < 4; u++) {
				/* put the rest of the matrix */
				for (j = cc[u]; j < cc[u + 1]; j++) {
					pixel[j] = colors[3 + u - cj];
				}
			}
		}

		/* scatters row i to black pixels */
		for (p = Ap[i]; p < Ap[i + 1]; p++) {
			pixel[Aj[p]] = 0;
		}

		/* dump the "pixel" array */
		for (j = 0; j < m; j++) {
			fprintf(f, "%d %d %d ", (pixel[j] >> 16) & 0xFF, (pixel[j] >> 8) & 0xFF, pixel[j] & 0xFF);
			t++;
			if ((t & 7) == 0) {
				fprintf(f, "\n");
			}
		}
	}
	fprintf(f, "\n");
}

/*
 * Saves a PPM (color pixmap) of specified dimensions of A, with an optional
 * DM decomposition
 */
void spasm_save_ppm(FILE * f, const spasm * A, const spasm_dm * X) {
	int i, j, n, m, t;
	int *Aj, *Ap, *rr, *cc, *pixel;
	spasm_cc *H, *S, *V;

	int colors[13] = {0, 0xFF0000, 0xFF6633, 0xCC0000, 0x990000,
		0xFFFF66, 0xFFCC00, 0xCC9900,
	0x669933, 0x99FF99, 0x33CC00};

	assert(f != NULL);
	assert(A != NULL);
	assert(X != NULL);

	Aj = A->j;
	Ap = A->p;
	n = A->n;
	m = A->m;

	pixel = spasm_malloc(m * sizeof(int));

	rr = X->DM->rr;
	cc = X->DM->cc;

	fprintf(f, "P3\n");
	fprintf(f, "%d %d\n", m, n);
	fprintf(f, "255\n");

	H = X->H;
	S = X->S;
	V = X->V;

	t = 0;

	/* --- H ---- */
	if (H != NULL) {
		render_block(f, m, Ap, Aj, H->CC, H->SCC, rr, cc, 0, 1, 0, 2, colors, pixel);
	}
	/* --- S ---- */
	if (S != NULL) {
		render_block(f, m, Ap, Aj, S->CC, S->SCC, rr, cc, 1, 2, 2, 3, colors + 4, pixel);
	}
	/* --- V ---- */
	if (V != NULL) {
		render_block(f, m, Ap, Aj, V->CC, V->SCC, rr, cc, 2, 4, 3, 4, colors + 8, pixel);
	} else {
		/* if V is empty, we need to print the empty rows */
		t = 0;
		for (i = rr[2]; i < rr[4]; i++) {
			for (j = 0; j < m; j++) {
				fprintf(f, "255 255 255 ");
				t++;
				if ((t & 7) == 0) {
					fprintf(f, "\n");
				}
			}
		}
		fprintf(f, "\n");
	}


	fprintf(f, "\n");
	free(pixel);
}
#endif