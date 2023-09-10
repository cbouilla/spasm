include_directories (../src)

############### tools

add_executable(bitmap bitmap.c)
target_link_libraries(bitmap PUBLIC spasm)

add_executable(vertical_swap vertical_swap.c)
target_link_libraries(vertical_swap PUBLIC spasm)

add_executable(stack stack.c)
target_link_libraries(stack PUBLIC spasm)

add_executable(transpose transpose.c)
target_link_libraries(transpose PUBLIC spasm)

add_executable(mtx_reader mtx_reader.c)
target_link_libraries(mtx_reader PUBLIC spasm)


############# core algorithms

add_executable(dm dm.c)
target_link_libraries(dm PUBLIC spasm)

add_executable(rank rank.c)
target_link_libraries(rank PUBLIC spasm)

add_executable(echelonize echelonize.c)
target_link_libraries(echelonize PUBLIC spasm)

add_executable(kernel kernel.c)
target_link_libraries(kernel PUBLIC spasm)
