add_library(libspasm spasm_util.c mmio.c spasm_io.c spasm_triplet.c
	spasm_GFp.c spasm_transpose.c spasm_submatrix.c 
	spasm_scatter.c spasm_gaxpy.c spasm_permutation.c spasm_reach.c 
	spasm_triangular.c spasm_lu.c spasm_solutions.c spasm_pivots.c 
	spasm_matching.c spasm_dm.c spasm_cc.c spasm_scc.c spasm_kernel.c 
	spasm_concatenate.c spasm_dense_lu.c spasm_schur.c)
target_link_libraries(libspasm PUBLIC OpenMP::OpenMP_C)
target_link_libraries(libspasm PUBLIC m)

