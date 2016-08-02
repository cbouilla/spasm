#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include "spasm.h"

/* generate .pbm or .pgm (=netlib) images. */

int main(int argc, char **argv) {
	spasm_triplet *T;
	spasm *A;
	int ch, w = -1, h = -1, gray = 0;
	double mpix = -1, n, m;

	/* options descriptor */
	struct option longopts[5] = {
		{"gray", no_argument, NULL, 'g'},
		{"width", required_argument, NULL, 'w'},
		{"height", required_argument, NULL, 'h'},
		{"mpixels", required_argument, NULL, 'm'},
		{NULL, 0, NULL, 0}
	};

	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'g':
			gray = 1;
			break;
		case 'w':
			w = atoi(optarg);
			break;
		case 'h':
			h = atoi(optarg);
			break;
		case 'm':
			mpix = atof(optarg);
			break;
		default:
			printf("Unknown option\n");
			exit(1);
		}
	}
	argc -= optind;
	argv += optind;

	T = spasm_load_sms(stdin, -1);
	A = spasm_compress(T);
	spasm_triplet_free(T);
	n = A->n;
	m = A->m;

	/* compute the output size */
	if (mpix > 0 && (w > 0 || h > 0)) {
		fprintf(stderr, "--mpixels and --width/--height are mutually exclusive\n");
		exit(1);
	}
	if (mpix < 0 && w < 0 && h < 0) {
		w = m;
		h = n;
	}
	if (mpix > 0 && w < 0 && h < 0) {
		double alpha = sqrt((mpix * 1e6) / (n * m));
		w = alpha * m;
		h = alpha * n;
		fprintf(stderr, "[bitmap] targeting %.1f Mpixels ; w=%d, h=%d\n", mpix, w, h);
	}
	if (mpix < 0 && w < 0 & h > 0)
		w = m * (h / n);

	if (mpix < 0 && w > 0 & h < 0)
		h = n * (w / m);

	w = spasm_min(w, m);
	h = spasm_min(h, n);

	if (gray)
		spasm_save_pgm(stdout, w, h, A);
	else
		spasm_save_pbm(stdout, w, h, A);

	spasm_csr_free(A);
	return 0;
}
