#include <assert.h>
#include <stdio.h>
#include <getopt.h>
#include "spasm.h"

/* computes a chain decomposition of the graph described by a symmetric matrix. */

int main() 
{
	int prime = -1;	
	spasm_triplet * T = spasm_load_mm(stdin, prime);
	assert (T->n == T->m);

	spasm * A = spasm_compress(T);
	spasm_triplet_free(T);
	int n = A->n;
	int *Ap = A->p;
	int *Aj = A->j;

	/**************** DFS ****************/
	int *stack = malloc(sizeof(int) * n);
	int *vmark = malloc(sizeof(int) * n);
	int *pstack = malloc(sizeof(int) * n);
	int *parent = malloc(sizeof(int) * n);
	for (int i = 0; i < n; i++)
		parent[i] = -1;
	for (int i = 0; i < n; i++)
		vmark[i] = 0;
	int top = n;
	int head = 0;
	stack[0] = 0;
	pstack[0] = Ap[0];
	vmark[0] = 1;
	while (head >= 0) {
		int i;
	start:
		i = stack[head];
		for (int px = pstack[head]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			if (vmark[j])
				continue;
			pstack[head] = px + 1;
			head += 1;
			stack[head] = j;
			pstack[head] = Ap[j];
			parent[j] = i;
			vmark[j] = 1;
			goto start;
		}
		top -= 1;
		stack[top] = stack[head];
		head -= 1;
	}
	if (top != 0) {
		fprintf(stderr, "non-connected graph\n");
		exit(42);
	}
	free(pstack);
	/**************** Form Oriented Graph ****************/
	int *level = malloc(sizeof(int) * n);
	for (int i = 0; i < n; i++)
		level[stack[i]] = i;
	T = spasm_triplet_alloc(n, n, A->nzmax, A->prime, 0);
	for (int i = 0; i < n; i++)
		for (int px = Ap[i]; px < Ap[i + 1]; px++) {
			int j = Aj[px];
			if ((level[j] >= level[i]  && parent[j] != i) || parent[i] == j)
				spasm_add_entry(T, i, j, 0);
		}
	spasm *G = spasm_compress(T);
	spasm_triplet_free(T);
	int *Gp = G->p;
	int *Gj = G->j;
	for (int i = 0; i < n; i++)
		for (int px = Gp[i]; px < Gp[i + 1]; px++)
			if (parent[i] == Gj[px]) {
				spasm_swap(Gj, Gp[i], px);
				break;
			}

	/**************** Chain Decomposition ****************/
	int n_cut_vertex = 0;
	int n_bridge = 0;

	int *emark = malloc(sizeof(int) * Gp[n]);
	for (int i = 0; i < n; i++)
		vmark[i] = -1;
	for (int i = 0; i < Gp[n]; i++)
		emark[i] = -1;

	int chain = 0;
	for (int i = 0; i < n; i++) {
		int u = stack[i];
		for (int px = (parent[u] >= 0) ? Gp[u] + 1 : Gp[u]; px < Gp[u + 1]; px++) {
			int v = Gj[px];
			if (vmark[u] < 0) {
				vmark[u] = chain;
				if (i > 0)
					// printf("%d is a cut-vertex (start of a cycle)\n", u + 1);
					n_cut_vertex++;
			}
			emark[px] = chain;
			while (vmark[v] < 0) {
				vmark[v] = chain;
				assert(emark[Gp[v]] < 0);
				emark[Gp[v]] = chain;
				v = Gj[Gp[v]];
			}
			chain += 1;
		}
	}
	//assert(chain == spasm_nnz(G) - n + 1);

	/*************** reap of the rewards *********************/
	for (int i = 0; i < n; i++)
		for (int px = Gp[i]; px < Gp[i + 1]; px++)
			if (emark[px] < 0)
				// printf("edge %d -- %d is a bridge (and both vertices are cut vertices)\n", i + 1, Gj[px] + 1);
				n_bridge++;

	printf("%d ; %d\n", n_cut_vertex, n_bridge);
	spasm_csr_free(G);
	spasm_csr_free(A);
	return 0;
}
