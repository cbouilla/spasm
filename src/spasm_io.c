#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <err.h>
#include <string.h>
#include <ctype.h>

#include "spasm.h"

static bool read_line(const char *fn, i64 line, char *buffer, int size, spasm_sha256_ctx *ctx, FILE *f)
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

static void validate_mm_header(const char *buffer)
{
	char mtx[1024], crd[1024], data_type[1024], storage_scheme[1024];
	
	if (sscanf(buffer, "%%%%MatrixMarket %s %s %s %s", mtx, crd, data_type, storage_scheme) != 4)
        errx(1, "incomplete MatrixMarket header");

    for (char *p = mtx; *p != '\0'; *p = tolower(*p),p++);  /* convert to lower case */
    for (char *p = crd; *p != '\0'; *p = tolower(*p),p++);  
    for (char *p = data_type; *p != '\0'; *p = tolower(*p),p++);
    for (char *p = storage_scheme; *p != '\0'; *p = tolower(*p),p++);

    /* first field should be "mtx" */
    if (strcmp(mtx, "matrix") != 0)
        errx(1, "unsupported MatrixMarket object type %s (I only know about ``matrix'')", mtx);
    
    if (strcmp(crd, "coordinate") != 0)
    	errx(1, "unsupported MatrixMarket format %s (I only know about ``coordinate'')", crd);

    if (strcmp(data_type, "integer") != 0)
        errx(1, "unsupported MatrixMarket data type %s (I only know about ``integer'')", data_type);
    
    if (strcmp(storage_scheme, "general") != 0)
        errx(1, "unsupported MatrixMarket storage scheme %s (I only know about ``general'')", storage_scheme);
}

/*
 * load a matrix in SMS format from f (an opened file, possibly stdin). 
 * set prime == -1 to avoid loading values.
 * if hash != NULL, then the SHA256 of the input matrix is written in hash (32 bytes)
 */
struct spasm_triplet *spasm_triplet_load(FILE * f, i64 prime, u8 *hash)
{
	assert(f != NULL);
	double start = spasm_wtime();
	int i, j;
	i64 nnz = 1;
	char buffer[1024];
	char hnnz[16];

	spasm_sha256_ctx ctx_always;
	spasm_SHA256_init(&ctx_always);
	spasm_sha256_ctx *ctx = (hash != NULL) ? &ctx_always : NULL;

	/* Process header */
	i64 line = 0;
	bool eof = read_line("spasm_triplet_load", line, buffer, 1024, ctx, f);
	if (eof)
		errx(1, "[spasm_triplet_load] empty file\n");

	/* MatrixMarket ? */
	bool mm = 0;
	if (strncmp(buffer, "%%MatrixMarket", 14) == 0) {
		/* MatrixMarket Header */ 
		mm = 1;
		validate_mm_header(buffer);
		/* skip comments, read dimensions */
		for (;;) {
			line += 1;
			eof = read_line("spasm_triplet_load", line, buffer, 1024, ctx, f);
			if (eof)
				errx(1, "premature EOF on line %" PRId64 " (expected matrix dimensions)", line);
			if (buffer[0] != '%')
				break;
		}
		if (sscanf(buffer, "%d %d %" SCNd64 "\n", &i, &j, &nnz) != 3)
			errx(1, "[spasm_triplet_load] bad MatrixMarking dimensions (line %" PRId64 ")\n", line);
		spasm_human_format(nnz, hnnz);
		logprintf("[IO] loading %d x %d MatrixMarket matrix modulo %" PRId64 " with %s non-zero... ", 
			i, j, prime, hnnz);
		fflush(stderr);
	} else {
		/* SMS header */
		char type;
		if (sscanf(buffer, "%d %d %c\n", &i, &j, &type) != 3)
			errx(1, "[spasm_triplet_load] bad SMS file (header)\n");
		if (prime != -1 && type != 'M')
			errx(1, "[spasm_triplet_load] only ``Modular'' type supported\n");
		logprintf("[IO] loading %d x %d SMS matrix modulo %" PRId64 "... ", i, j, prime);
		fflush(stderr);
	}

	/* allocate result */
	struct spasm_triplet *T = spasm_triplet_alloc(i, j, nnz, prime, prime != -1);

	i64 x;
	bool end = 0;
	i64 entries = 0;
	for (;;) { 
		line += 1;
		eof = read_line("spasm_triplet_load", line, buffer, 1024, ctx, f);
		if (end && eof)
			break;
		if (end && !eof) {
			warn("[spasm_load] garbage detected near end of file");
			continue;
		}
		if (!end && eof)
			errx(1, "[spasm_triplet_load] premature end of file (line %" PRId64 ", read %" PRId64" nz)", line, entries);

		if (sscanf(buffer, "%d %d %" SCNd64 "\n", &i, &j, &x) != 3)
			errx(1, "parse error line %" PRId64, line);
		if (i == 0 && j == 0 && x == 0) {
			if (mm)
				errx(1, "SMS end marker in MatrixMarket file");
			end = 1;
		}
		if (!end) {
			spasm_add_entry(T, i - 1, j - 1, x);
			entries += 1;
		}
		if (mm && entries == nnz)
			end = 1;
	}

	if (!mm) {
		spasm_triplet_realloc(T, -1);
		spasm_human_format(T->nz, hnnz);
		logprintf("%s non-zero [%.1fs]\n", hnnz, spasm_wtime() - start);
	} else {
		logprintf("[%.1fs]\n", spasm_wtime() - start);
	}

	if (ctx != NULL) {
		spasm_SHA256_final(hash, ctx);
		logprintf("[spasm_triplet_load] sha256(matrix) = ");
		for (int i = 0; i < 32; i++)
			logprintf("%02x", hash[i]);
		logprintf(" / size = %" PRId64" bytes\n", (((i64) ctx->Nh) << 29) + ctx->Nl / 8);
	}
	return T;
}

/*
 * save a matrix in SMS format. TODO : change name to spasm_csr_save
 */
void spasm_csr_save(const struct spasm_csr *A, FILE *f)
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
void spasm_triplet_save(const struct spasm_triplet *A, FILE *f)
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
void spasm_save_pnm(const struct spasm_csr *A, FILE *f, int x, int y, int mode, struct spasm_dm *DM)
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
		logprintf("limits_v = %d, %d\n", limits_v[0], limits_v[1]);
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
