#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <err.h>
#include <string.h>

#include "spasm.h"
#include "mmio.h"

static bool read_line(const char *fn, i64 line, char *buffer, int size, SHA256_CTX *ctx, FILE *f)
{
	if (fgets(buffer, size, f) == NULL) {
		if (feof(f))
			return 1;
		err(1, "[%s] impossible to read line %" PRId64, fn, line);
	}
	int l = strlen(buffer);
	if (l == 0)
		errx(1, "[%s] empty line %" PRId64, fn, line);
 	if (buffer[l - 1] != '\n' && !feof(f))
 		errx(1, "[%s] line %" PRId64 " too long (> %d)", fn, line, size);
 	if (ctx != NULL)
 		spasm_SHA256_update(ctx, buffer, l);
 	return 0;
}

/*
 * load a matrix in SMS format from f (an opened file, possibly stdin). 
 * set prime == -1 to avoid loading values.
 * if hash != NULL, then the SHA256 of the input matrix is written in hash (32 bytes)
 */
spasm_triplet *spasm_load_sms(FILE * f, i64 prime, u8 *hash)
{
	assert(f != NULL);
	double start = spasm_wtime();
	int i, j;
	char type;
	char buffer[256];
	SHA256_CTX ctx_always;
	spasm_SHA256_init(&ctx_always);
	SHA256_CTX *ctx = (hash != NULL) ? &ctx_always : NULL;

	/* Process header */
	i64 line = 0;
	bool eof = read_line("spasm_load_sms", line, buffer, 256, ctx, f);
	if (eof)
		errx(1, "[spasm_load_sms] empty file\n");
	if (sscanf(buffer, "%d %d %c\n", &i, &j, &type) != 3)
		errx(1, "[spasm_load_sms] bad SMS file (header)\n");
	if (prime != -1 && type != 'M')
		errx(1, "[spasm_load_sms] only ``Modular'' type supported\n");

	fprintf(stderr, "[IO] loading %d x %d SMS matrix modulo %" PRId64 "... ", i, j, prime);
	fflush(stderr);

	/* allocate result */
	spasm_triplet *T = spasm_triplet_alloc(i, j, 1, prime, prime != -1);

	i64 x;
	bool end = 0;
	for (;;) { 
		line += 1;
		eof = read_line("spasm_load_sms", line, buffer, 256, ctx, f);
		if (end && eof)
			break;
		if (end && !eof) {
			warn("[spasm_load_sms] garbage detected near end of file");
			continue;
		}
		if (!end && eof)
			errx(1, "[spasm_load_sms] premature end of file");

		if (sscanf(buffer, "%d %d %" SCNd64 "\n", &i, &j, &x) != 3)
			errx(1, "parse error line %" PRId64, line);
		if (i == 0 && j == 0 && x == 0)
			end = 1;
		if (!end)
			spasm_add_entry(T, i - 1, j - 1, x);
	}

	char nnz[16];
	spasm_human_format(T->nz, nnz);
	fprintf(stderr, "%s NNZ [%.1fs]\n", nnz, spasm_wtime() - start);
	if (ctx != NULL) {
		spasm_SHA256_final(hash, ctx);
		fprintf(stderr, "[spasm_load_sms] sha256(matrix) = ");
		for (int i = 0; i < 32; i++)
			fprintf(stderr, "%02x", hash[i]);
		fprintf(stderr, " / size = %" PRId64" bytes\n", (((i64) ctx->Nh) << 29) + ctx->Nl / 8);
	}
	return T;
}

/*
 * Load a matrix in MatrixMarket sparse format.
 * Heavily inspired by the example program:
 *     http://math.nist.gov/MatrixMarket/mmio/c/example_read.c
 */
spasm_triplet *spasm_load_mm(FILE *f, i64 prime)
{
	MM_typecode matcode;
	int n, m, nnz;

	double start = spasm_wtime();
	if (mm_read_banner(f, &matcode) != 0) 
		errx(1, "Could not process Matrix Market banner.\n");

	if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode))
		errx(1, "Matrix Market type: [%s] not supported", mm_typecode_to_str(matcode));
	
	int symmetric = mm_is_symmetric(matcode);
	int skew = mm_is_skew(matcode);
	if (!mm_is_general(matcode) && !symmetric && !skew)
		errx(1, "Matrix market type [%s] not supported",  mm_typecode_to_str(matcode));
	if (mm_read_mtx_crd_size(f, &n, &m, &nnz) != 0)
		errx(1, "Cannot read matrix size");
	fprintf(stderr, "[IO] loading %d x %d MTX [%s] modulo %" PRId64 ", %d nnz...", n, m, mm_typecode_to_str(matcode), prime, nnz);
	fflush(stderr);
	
	if (mm_is_pattern(matcode))
		prime = -1;

	spasm_triplet *T = spasm_triplet_alloc(n, m, nnz, prime, prime != -1);

	for (int i = 0; i < nnz; i++) {
		int u, v;
		i64 w;

		if (mm_is_pattern(matcode)) {
			if (2 != fscanf(f, "%d %d\n", &u, &v))
				errx(1, "parse error entry %d\n", i);
			spasm_add_entry(T, u - 1, v - 1, 1);
		} else if (mm_is_integer(matcode)) {
			if (3 != fscanf(f, "%d %d %" SCNd64 "\n", &u, &v, &w))
				errx(1, "parse error entry %d\n", i);
			spasm_add_entry(T, u - 1, v - 1, w);
		} else {
			errx(1, "Don't know how to read matrix");
		}
	}

	if (symmetric || skew)
		nnz *= 2;

	if (symmetric) {
		int mult = skew ? -1 : 1;
		int nz = T->nz;
		for (int px = 0; px < nz; px++)
			if (T->j[px] != T->i[px])
				spasm_add_entry(T, T->j[px], T->i[px], (T->x != NULL) ? (mult * T->x[px]) : 1);
	}

	char s_nnz[16];
	spasm_human_format(T->nz, s_nnz);
	fprintf(stderr, "%s NNZ [%.1fs]\n", s_nnz, spasm_wtime() - start);
	return T;
}

/*
 * save a matrix in SMS format. TODO : change name to spasm_csr_save
 */
void spasm_save_csr(FILE *f, const spasm *A)
{
	assert(f != NULL);
	const int *Aj = A->j;
	const i64 *Ap = A->p;
	const spasm_ZZp *Ax = A->x;
	int n = A->n;
	int m = A->m;

	fprintf(f, "%d %d M\n", n, m);
	for (int i = 0; i < n; i++)
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			i64 x = (Ax != NULL) ? Ax[px] : 1;
			fprintf(f, "%d %d %" PRId64 "\n", i + 1, Aj[px] + 1, x);
		}
	fprintf(f, "0 0 0\n");
}

/*
 * save a matrix in SMS format. TODO : change name to spasm_triplet_save
 */
void spasm_save_triplet(FILE *f, const spasm_triplet *A)
{
	assert(f != NULL);
	const int *Ai = A->i;
	const int *Aj = A->j;
	const spasm_ZZp *Ax = A->x;
	i64 nz = A->nz;
	fprintf(f, "%d %d M\n", A->n, A->m);
	for (i64 px = 0; px < nz; px++)
		fprintf(f, "%d %d %d\n", Ai[px] + 1, Aj[px] + 1, (Ax != NULL) ? Ax[px] : 1);
	fprintf(f, "0 0 0\n");
}

/* Saves a PBM (bitmap) of specified dimensions of A.
 * Mode: 1 --> create a PBM file (B/W bitmap)
 *       2 --> create a PGM file (gray levels) 
 *       3 --> create a PNM file (colors)
 */
void spasm_save_pnm(const spasm *A, FILE *f, int x, int y, int mode, spasm_dm *DM)
{
	const int *Aj = A->j;
	const i64 *Ap = A->p;
	int n = A->n;
	int m = A->m;
	x = spasm_min(x, m);
	y = spasm_min(y, n);
	int *w = spasm_malloc(x * y * sizeof(int));
	assert(f != NULL);
	assert((mode-1)*(mode-2)*(mode-3) == 0);
	assert((mode != 3) || (DM != NULL));
	for (int i = 0; i < x*y; i++)
		w[i] = 0;

	fprintf(f, "P%d\n", mode);
	fprintf(f, "%d %d\n", x, y);
	if (mode > 1)
		fprintf(f, "255\n");

	for (i64 i = 0; i < n; i++) {
		i64 k = i * y / n;
		int *ww = w + k * x;
		for (i64 px = Ap[i]; px < Ap[i + 1]; px++) {
			int t = ((i64) Aj[px]) * x / m;
			ww[t] += 1;
		}
	}

	double max = 0;
	for (int j = 0; j < x * y; j++)
		max = spasm_max(max, w[j]);

	int bgcolor[3][3] = {
		{0xFF0000, 0xCC0000, 0x990000},
		{0xFFFFFF, 0xFFCC00, 0xCC9900},    /* 0xFFFF66 inside SCC */
	    {0xFFFFFF, 0xFFFFFF, 0x33CC00}
	};

	int limits_h[2] = {};
	int limits_v[2] = {};
	int *r = NULL, *c = NULL;

	if (DM != NULL) {
		int *rr = DM->rr;
		int *cc = DM->cc;
		limits_h[0] = ((long long int) cc[2]) * ((long long int) x) / ((long long int) m);
		limits_h[1] = ((long long int) cc[3]) * ((long long int) x) / ((long long int) m);
		limits_v[0] = ((long long int) rr[1]) * ((long long int) y) / ((long long int) n);
		limits_v[1] = ((long long int) rr[2]) * ((long long int) y) / ((long long int) n);
		fprintf(stderr, "limits_v = %d, %d\n", limits_v[0], limits_v[1]);
		r = DM->r;
		c = DM->c;
	}

	int t = 0;
	double intensity;
	int scc = 0;
	int scc_left_edge = 0;
	int scc_right_edge = 0;
	int scc_bot_edge = 0;
	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			switch(mode) {
			case 1:
				fprintf(f, "%d ", (w[i * x + j] > 0) ? 1 : 0);
				break;
			case 2: 
				intensity = 1.0 - exp(0.1 * log(w[i * x + j] / max));
				/* intensity = 1.0 - ((double) w[i * x + j]) / max; */
				assert(0 <= intensity && intensity <= 1.0);
				fprintf(f, "%.0f ", 255.0 * intensity);
				break;
			case 3:
			        ;/* find out which blocks we are in */
				int block_h = 0;
				int block_v = 0;
				if (limits_v[0] <= i)
					block_v = (i < limits_v[1]) ? 1 : 2;
				if (limits_h[0] <= j)
					block_h = (j < limits_h[1]) ? 1 : 2;

				int bg = bgcolor[block_v][block_h];

				if (block_h == 1 && block_v == 1) {
					/* inside S, we need to take care of SCC's */
					while(scc_bot_edge <= i) {
						scc_left_edge = scc_right_edge;
						++scc;
						scc_right_edge = ((long long int) c[scc]) * ((long long int) x) / ((long long int) m);
						scc_bot_edge = ((long long int) r[scc]) * ((long long int) y) / ((long long int) n);
					}

					if (j < scc_left_edge)
						bg = 0xFFFFFF;
					else if (j < scc_right_edge)
						bg += 0x003366;
				}

				int pixel = (w[i * x + j] > 0) ? 0 : bg;
				fprintf(f, "%d %d %d ", (pixel >> 16) & 0xFF, (pixel >> 8) & 0xFF, pixel & 0xFF);
			}
			
			if (((t++) & 31) == 0)
				fprintf(f, "\n");
		}
	}

	free(w);
}
