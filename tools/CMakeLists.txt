include_directories (../src)

add_executable(bitmap bitmap.c)
target_link_libraries(bitmap PUBLIC libspasm)

add_executable(dm dm.c)
target_link_libraries(dm PUBLIC libspasm)

add_executable(vertical_swap vertical_swap.c)
target_link_libraries(vertical_swap PUBLIC libspasm)

add_executable(stack stack.c)
target_link_libraries(stack PUBLIC libspasm)

add_executable(transpose transpose.c)
target_link_libraries(transpose PUBLIC libspasm)

add_executable(rank rank.c)
target_link_libraries(rank PUBLIC libspasm)

add_executable(mtx_reader mtx_reader.c)
target_link_libraries(mtx_reader PUBLIC libspasm)

add_executable(kernel kernel.c)
target_link_libraries(kernel PUBLIC libspasm)
