include_directories (../src)

add_custom_target(check ${CMAKE_CTEST_COMMAND} -LE LongTest)                # create "make check"
add_compile_options(-Wno-unused-parameter)

add_custom_target(longtests ${CMAKE_CTEST_COMMAND} -L LongTest)            # create "make long-tests"

############# custom functions to declare and run tests easily ################

add_library(spasmtest test_tools.c)

function(spasm_declare_test name)
    add_executable(test_${name} EXCLUDE_FROM_ALL ${name}.c)
    target_link_libraries(test_${name} PUBLIC spasm)
    target_link_libraries(test_${name} PUBLIC spasmtest)
    target_link_libraries(test_${name} PUBLIC OpenMP::OpenMP_C)
    add_dependencies(check test_${name})
endfunction()

function(spasm_run_test name modulus test_matrix)
    add_test(NAME ${name}-${modulus}-${test_matrix} 
             COMMAND sh -c "./test_${name} --modulus ${modulus} < ${CMAKE_CURRENT_SOURCE_DIR}/Matrix/${test_matrix}")
    set_tests_properties(${name}-${modulus}-${test_matrix} PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")
    set_tests_properties(${name}-${modulus}-${test_matrix} PROPERTIES FAIL_REGULAR_EXPRESSION "not ok")
    set_tests_properties(${name}-${modulus}-${test_matrix} PROPERTIES TIMEOUT 1)
endfunction()

function(spasm_test_expected_output name test_matrix expected_output)
    add_test(NAME ${name}-${test_matrix} 
             COMMAND sh -c "./test_${name} < ${CMAKE_CURRENT_SOURCE_DIR}/Matrix/${test_matrix} | diff - ${CMAKE_CURRENT_SOURCE_DIR}/Expected/${expected_output}")
endfunction()

set(GITLAB_URL "https://gitlab.lip6.fr/almasty/spasm-test-matrices/-/raw/main")

function(spasm_long_test name test_matrix opts)
    add_test(NAME ${name}-${test_matrix}-${opts}
             COMMAND sh -c "curl --silent ${GITLAB_URL}/${test_matrix} | zcat | ./test_${name} ${opts}")
    set_tests_properties(${name}-${test_matrix}-${opts} PROPERTIES LABELS "LongTest")
    set_tests_properties(${name}-${test_matrix}-${opts} PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP")
    set_tests_properties(${name}-${test_matrix}-${opts} PROPERTIES FAIL_REGULAR_EXPRESSION "not ok")
endfunction()

spasm_long_test(lu kneser_10_4_1.sms.gz "")
spasm_long_test(lu boundary_C_6_9.sms.gz "")

set(ALL_MODULI
    3
    257
    65537
    67108859
    189812507
    4294967291
)

function(spasm_run_tests_mod name test_matrices)
    foreach (test_matrix ${test_matrices})
        foreach (modulus ${ALL_MODULI})
            spasm_run_test(${name} ${modulus} ${test_matrix})
        endforeach (modulus)
    endforeach (test_matrix)
endfunction()

function(spasm_run_tests name test_matrices)
    foreach (test_matrix ${test_matrices})
            spasm_run_test(${name} 257 ${test_matrix})
    endforeach (test_matrix)
endfunction()

########################## now the actual tests ###############################

set(ALL_TEST_MATRICES
    BIOMD0000000424.int.mpl.sms
    cc.sms
    chains.sms 
    dm.sms     
    example.sms  
    G2.sms  
    lower_trapeze.sms  
    mat364.sms  
    p3.sms             
    rectangular_l.sms  
    scc3.sms  
    singular2.sms  
    singular.sms    
    t1.sms         
    # trefethen_2000.sms  
    # tridiagonal.sms  
    upper_trapeze.sms
    BIOMD0000000525.int.mpl.sms  
    chains2.sms  
    dm2.sms     
    empty.sms  
    g1.sms       
    l1.sms  
    m1.sms             
    medium.sms  
    rectangular_h.sms  
    scc2.sms           
    scc.sms   
    singular3.sms  
    small.sms       
    tree_test.sms  
    trefethen_500.sms   
    u1.sms           
    void.sms
)

######### GFp arith

spasm_declare_test(GFp)
add_test(NAME GFp COMMAND test_GFp)

########## prng / hash

spasm_declare_test(prng)
add_test(NAME prng COMMAND sh -c "./test_prng | diff - ${CMAKE_CURRENT_SOURCE_DIR}/Expected/prng")

spasm_declare_test(sha)
add_test(NAME sha COMMAND sh -c "./test_sha | diff - ${CMAKE_CURRENT_SOURCE_DIR}/Expected/hash")

########## Dulmage-Mendelson / matching / SCC / etc.

spasm_declare_test(dm)
spasm_declare_test(matching)
spasm_declare_test(scc)

set(DM_TEST_MATRICES
    scc.sms
    scc2.sms
    scc3.sms
    mat364.sms
    t1.sms
    dm.sms
    dm2.sms
    BIOMD0000000424.int.mpl.sms 
    BIOMD0000000525.int.mpl.sms
)
spasm_run_tests(matching "${DM_TEST_MATRICES}")
spasm_run_tests(scc      "${DM_TEST_MATRICES}")
spasm_run_tests(dm       "${DM_TEST_MATRICES}")

########## permutations

spasm_declare_test(vec_perm)
spasm_declare_test(mat_perm)

add_test(NAME vec_perm COMMAND test_vec_perm)
spasm_run_test(mat_perm 65537 small.sms)
spasm_run_test(mat_perm 65537 upper_trapeze.sms)

########## transpose / submatrix / etc.

spasm_declare_test(transpose)
spasm_run_tests(transpose "${ALL_TEST_MATRICES}")

spasm_declare_test(spmv)
spasm_test_expected_output(spmv m1.sms gaxpy.1)

spasm_declare_test(submatrix)
spasm_test_expected_output(submatrix singular.sms submatrix.1)

########### FFPACK

spasm_declare_test(dense_rref_ffpack)
spasm_declare_test(dense_lu_ffpack)

spasm_run_tests_mod(dense_lu_ffpack   "${ALL_TEST_MATRICES}")
spasm_run_tests_mod(dense_rref_ffpack   "${ALL_TEST_MATRICES}")

########## triangular solvers

set(LOWER_TEST_MATRICES
    l1.sms
    lower_trapeze.sms
)
set(UPPER_TEST_MATRICES
    u1.sms
    upper_trapeze.sms
)

spasm_declare_test(dense_lsolve)
spasm_declare_test(sparse_lsolve)
spasm_declare_test(dense_usolve)
spasm_declare_test(sparse_usolve)
spasm_declare_test(sparse_lu_usolve)
spasm_declare_test(sparse_utsolve)

spasm_run_tests_mod(dense_lsolve "${LOWER_TEST_MATRICES}")
spasm_run_tests_mod(sparse_lsolve "${LOWER_TEST_MATRICES}")
spasm_run_tests_mod(dense_usolve "${UPPER_TEST_MATRICES}")
spasm_run_tests_mod(sparse_usolve "${UPPER_TEST_MATRICES}")
spasm_run_tests_mod(sparse_lu_usolve "${ALL_TEST_MATRICES}")
spasm_run_tests_mod(sparse_utsolve "${ALL_TEST_MATRICES}")

########## schur complement

spasm_declare_test(schur)
spasm_declare_test(schur_dense)
spasm_run_tests_mod(schur       "${ALL_TEST_MATRICES}")
spasm_run_tests_mod(schur_dense "${ALL_TEST_MATRICES}")

########## echelonization

spasm_declare_test(echelonize)
spasm_run_tests_mod(echelonize       "${ALL_TEST_MATRICES}")

########## kernel

spasm_declare_test(kernel)
spasm_run_tests_mod(kernel "${ALL_TEST_MATRICES}")

########## lu / solve / dense LU

spasm_declare_test(lu)
spasm_declare_test(solve)
spasm_declare_test(gesv)

spasm_run_tests_mod(lu                "${ALL_TEST_MATRICES}")
spasm_run_tests_mod(solve             "${ALL_TEST_MATRICES}")
spasm_run_tests_mod(gesv              "${ALL_TEST_MATRICES}")

########## certificates

spasm_declare_test(rank_cert)
spasm_run_tests_mod(rank_cert    "${ALL_TEST_MATRICES}")
