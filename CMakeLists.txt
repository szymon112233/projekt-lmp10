cmake_minimum_required(VERSION 3.3)
project(lmp10)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    gaus/matrix.c
    gaus/matrix.h
    gaus/piv_ge_solver.h
    gaus/piv_ge_solver.c
    gaus/pivot.c
    gaus/pivot_ge_solver.c
    test/gen.c
    aproksymator_na_bazie.c
    interpolator.c
    main.c
    makespl.h
    points.c
    points.h
    prosta.c
    splines.c
    splines.h)

add_executable(lmp10 ${SOURCE_FILES})