cmake_minimum_required(VERSION 3.13)
project(DA_495 C)

set(CMAKE_C_STANDARD 99)

add_executable(DA_495 main.c ArrayMem.c ArrayMem.h Functions.c Functions.h mt19937ar.h mt19937ar.c SelectFunctions.c SelectFunctions.h DA.c DA.h Util.c Util.h)