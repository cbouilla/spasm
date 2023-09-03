include_directories (../src)

############### tools

add_executable(bitmap bitmap.c)
target_link_libraries(bitmap PUBLIC libspasm)

add_executable(vertical_swap vertical_swap.c)
target_link_libraries(vertical_swap PUBLIC libspasm)

add_executable(stack stack.c)
target_link_libraries(stack PUBLIC libspasm)

add_executable(transpose transpose.c)
target_link_libraries(transpose PUBLIC libspasm)

add_executable(mtx_reader mtx_reader.c)
target_link_libraries(mtx_reader PUBLIC libspasm)


############# core algorithms

add_executable(dm dm.c)
target_link_libraries(dm PUBLIC libspasm)

add_executable(find_pivots find_pivots.c)
target_link_libraries(find_pivots PUBLIC libspasm)

add_executable(rank rank.c)
target_link_libraries(rank PUBLIC libspasm)

add_executable(echelonize echelonize.c)
target_link_libraries(echelonize PUBLIC libspasm)

add_executable(kernel kernel.c)
target_link_libraries(kernel PUBLIC libspasm)
