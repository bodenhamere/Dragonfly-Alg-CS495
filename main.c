/**
 * @file main.c
 */
/**
 * @author Emily Bodenhamer
 *  CWU ID 41119306
 *  CS 471 Optimization Project 4
 *  Date 5/10/2019
 *
 *  This project implements three meta-heuristic optimization algorithms.
 *  Particle Swarm Optimization (PSO), Firefly Algorithm (FA), and Harmony Search Algorithm (HS).
 *
 */
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mt19937ar.h"
#include "Util.h"
#include "ArrayMem.h"
#include "DA.h"

int main() {
    // declare variables
    init_genrand(time(0));
    int NS = 0; //!< number of solutions
    int DIM = 0; //!< number of dimensions
    int iterations = 0; //!< number of iterations
    int fitnessCallCounter = 0; //!< function call counter

    // open our file
    FILE *DAFile;
    DAFile = fopen("DAFile.txt", "r");
    // check if our file exists
    if (DAFile == NULL) {
        printf("Error opening DAFile.txt file, "
               "make sure the file is in the cmake-build-debug folder\n");
        exit(-1);
    }

    // get the data from the files
    while (!feof(DAFile)) {
        // read values and save them to the DA struct
        fscanf(DAFile, "%d %d %d",
               &iterations, &NS,
               &DIM);
    }
    // close the file
    fclose(DAFile);

    // create structure object
    initData *myData;
    myData = (initData *) calloc(1, sizeof(initData));

    // starting DA
    readInput(myData, NS, DIM, iterations, fitnessCallCounter);

    return 0;
    // # simulations = > 30 < 50
    // 18 different *50 = 900
    // write a script to launch each simulation
    // create a python script to call it internally
}