include_directories (../src)

if(APPLE)
    include_directories (/opt/homebrew/include)
    set(ARGP argp)
endif()

add_library(spasmtools common.c)

############### tools

add_executable(bitmap bitmap.c)
target_link_libraries(bitmap PUBLIC spasm)

add_executable(vertical_swap vertical_swap.c)
target_link_libraries(vertical_swap PUBLIC spasm)

add_executable(stack stack.c)
target_link_libraries(stack PUBLIC spasm)

add_executable(transpose transpose.c)
target_link_libraries(transpose PUBLIC spasm)

############# core algorithms

add_executable(dm dm.c)
target_link_libraries(dm PUBLIC spasm)

add_executable(rank rank.c)
target_link_libraries(rank PUBLIC spasm spasmtools ${ARGP})

add_executable(echelonize echelonize.c)
target_link_libraries(echelonize PUBLIC spasm spasmtools)

add_executable(kernel kernel.c)
target_link_libraries(kernel PUBLIC spasm spasmtools ${ARGP})

add_executable(solve solve.c)
target_link_libraries(solve PUBLIC spasm spasmtools ${ARGP})

add_executable(check_cert check_cert.c)
target_link_libraries(check_cert PUBLIC spasm spasmtools ${ARGP})
