add_library(libspasm 
	# common stuff, IO, utilities
	spasm_util.c mmio.c spasm_io.c spasm_triplet.c spasm_GFp.c spasm_scatter.c
	spasm_spmv.c spasm_transpose.c spasm_permutation.c

	# core algo: triangular solver
	spasm_reach.c spasm_triangular.c 
	
	# core algo: echelonization
	spasm_pivots.c spasm_schur.c spasm_ffpack.cpp 
	spasm_echelonize.c spasm_rref.c spasm_kernel.c
	
	# core algo: dulmage-mendelson
	spasm_submatrix.c spasm_matching.c spasm_dm.c spasm_scc.c 

	# needs work
	# spasm_solutions.c
)

target_link_libraries(libspasm PUBLIC OpenMP::OpenMP_C)
target_link_libraries(libspasm PUBLIC m)
target_link_libraries(libspasm PUBLIC PkgConfig::GIVARO)
target_link_libraries(libspasm PUBLIC PkgConfig::FFLAS_FFPACK)
