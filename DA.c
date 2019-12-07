/**
 * @file main.c
 */
/**
 * @author Emily Bodenhamer
 *  CWU ID 41119306
 *  CS 495
 *  Date 9/15/2019
 *
 *  This project implements the Dragonfly meta-heuristic optimization algorithm.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "DA.h"
#include "Util.h"
#include "ArrayMem.h"
#include "SelectFunctions.h"
#include "mt19937ar.h"
#include <math.h>
#include <string.h>

/**
 * Function that reads from the DA input file to
 * initialize the structure variables for the DA algorithm
 * @param iterations number of iterations to run the algorithm for
 * @param fitnessCounter count how many times the fitness function was called
 * @param myData Data structure pointer for initializing the bounds and which function to run
 * @param NS number of solutions for a data type
 * @param DIM number of dimensions for a data type
 */
void readInput(initData *myData, int NS, int DIM, int iterations, int fitnessCounter) {
    // initialize our DA struct
    DA *myDA;
    myDA = (DA *) calloc(1, sizeof(DA));

    // create the file to fill the population
    FILE *funFile;
    funFile = fopen("Functions.txt", "r");
    // file to output results
    FILE *DAOut;
    DAOut = fopen("DA.csv", "w+");

    // check if our file exists
    if (funFile == NULL) {
        printf("Error opening Functions.txt file, make sure the file is in the cmake-build-debug folder\n");
        exit(1);
    }

    // get the data from the files
    while (!feof(funFile)) {
        // read values and save them to the data struct
        fscanf(funFile, "%lf %lf %d",
               &myData->min, &myData->max,
               &myData->functionNumber);

        // initializing the size of population and fitness
        myData->population = createDblArray(NS, DIM);
        myData->fitness = singleArray(NS);

        // initialize the size of the arrays
        myDA->step = createDblArray(NS, DIM);
        myDA->worstArr = singleArray(2);
        myDA->bestArr = singleArray(2);
        // Initialize the size of the neighboring arrays
        myDA->neighborsStep = createDblArray(NS, DIM);
        myDA->neighborsPop = createDblArray(NS, DIM);
        myDA->sVector = singleArray(DIM);
        myDA->aVector = singleArray(DIM);
        myDA->cVector = singleArray(DIM);
        myDA->eVector = singleArray(DIM);
        myDA->fVector = singleArray(DIM);

        // fill our population and step vectors with random numbers
        // and get the fitness of our population
        myData->population = fillIn(myData->population, NS, DIM, myData->min, myData->max);
        myDA->step = fillIn(myDA->step, NS, DIM, myData->min, myData->max);
        myData->fitness = getFun(myData->fitness, myData->population, NS, DIM, myData->functionNumber);


        // enemy and food value
        findWorst(myDA->worstArr, myData->fitness, NS);
        findBest(myDA->bestArr, myData->fitness, NS);
        myDA->enemy = myDA->worstArr[0];
        myDA->food = myDA->bestArr[0];

        // enemy and food position
        myDA->foodPos = myDA->worstArr[1];
        myDA->enemyPos = myDA->bestArr[1];
        // start the algorithm
        startDA(myDA, myData, NS, DIM, iterations, fitnessCounter, DAOut);

        // free memory
        free(myDA->worstArr);
        free(myDA->bestArr);
        freeMem(NS, myDA->step);
        free(myData->fitness);
        freeMem(NS, myData->population);
        freeMem(NS, myDA->neighborsStep);
        freeMem(NS, myDA->neighborsPop);
        free(myDA->sVector);
        free(myDA->aVector);
        free(myDA->cVector);
        free(myDA->fVector);
        free(myDA->eVector);
        free(myDA->o);
    }
    fclose(funFile);
    fclose(DAOut);
}

/**
 * Function that runs the DA algorithm
 * @param fileOut output file for the fitness of the population
 * @param NS number of solutions for a data type
 * @param DIM number of dimensions for a data type
 * @param iterations number of iterations to run the algorithm for
 * @param myData Data structure pointer for initializing the bounds and which function to run
 * @param fitnessCounter count how many times the fitness function was called
 * @param myDA the DA structure
*/
void startDA(DA *myDA, initData *myData, int NS, int DIM, int iterations, int fitnessCounter, FILE *fileOut) {
    // start the clock
    clock_t start;
    start = clock();
    // arrays to hold the best and worst solutions
    double *holdb = singleArray(iterations);
    double *holdw = singleArray(iterations);
    for (int t = 0; t < iterations; t++) {
        // update weights and radius
        updateWeights(myDA, myData, t, iterations);

        for (int i = 0; i < NS; i++) {
            myDA->o = singleArray(DIM);

            // update neighbors
            findNeighbors(myDA, myData, i, DIM, NS);

            // calculate S, A, C, F and E
            separation(myDA, myData, DIM, i);
            alignment(myDA, DIM, i);
            cohesion(myDA, myData, DIM, i);
            distraction(myDA, myData, i, DIM);
            attraction(myDA, myData, i, DIM);

            //update velocity and position vectors
            if (greaterR3(myDA, DIM)){
                if (myDA->numNeighbors > 1) {
                    updateStepPosition2(myDA, myData, i, DIM);
                } else {
                    // perform random walk
                    randomWalk(myDA, myData, i, DIM);
                }
            } else {
                updateStepPosition(myDA, myData, i, DIM);
            }


            // update the fitness of our new population
            myData->fitness[i] = getFunSingle(myData->fitness[i], myData->population[i], DIM, myData->functionNumber);
            fitnessCounter++;

            if (myData->fitness[i] < myDA->food) {
                myDA->food = myData->fitness[i];
                // enemy and food position
                myDA->foodPos = i;
            }
            if (myData->fitness[i] > myDA->enemy) {
                myDA->enemy = myData->fitness[i];
                // enemy and food position
                myDA->enemyPos = i;
            }

            for (int j = 0; j < NS; ++j) {
                memset(myDA->neighborsStep[j],0, sizeof(double));
                memset(myDA->neighborsPop[j],0, sizeof(double));
            }
            memset(myDA->sVector,0, sizeof(double));
            memset(myDA->aVector,0, sizeof(double));
            memset(myDA->cVector,0, sizeof(double));
            memset(myDA->fVector,0, sizeof(double));
            memset(myDA->eVector,0, sizeof(double));
            memset(myDA->o,0, sizeof(double));

        }
        printf("(%d)Best: ", t);
        printf("%lf ", myDA->food);
        printf("Worst: ");
        printf("%lf ", myDA->enemy);
        printf("\n");

        holdb[t] = myDA->food;
        holdw[t] = myDA->enemy;

    }
    printSingle(fileOut, holdb, iterations);
    printSingle(fileOut, holdw, iterations);
    start = (((clock() - start)));
    printf("\nExperiment for %d took, %lf, Counter, %d\n",
           myData->functionNumber, ((((double) start) / CLOCKS_PER_SEC) * 1000), fitnessCounter);

    fprintf(fileOut, "\nExperiment for %d took, %lf, Counter, %d\n",
            myData->functionNumber, ((((double) start) / CLOCKS_PER_SEC) * 1000), fitnessCounter);
    free(holdb);
    free(holdw);
}

/**
 * Function that updates the weights for the DA algorithm
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param iter number of iterations for the algorithm
 * @param maxIter max number of iterations for the algorithm
 */
void updateWeights(DA *myDA, initData *myData, int iter, int maxIter) {
    myDA->r = (myData->max - myData->min) / 4 + ((myData->max - myData->min) * (iter / maxIter) * 2);
    myDA->w = 0.9 - iter * ((0.9 - 0.4) / maxIter);
    myDA->weight = 0.1 - iter * ((0.1 - 0) / (maxIter / 2));
    if (myDA->weight < 0) {
        myDA->weight = 0;
    }
    myDA->s = 2 * genrand_real1() * myDA->weight;
    myDA->a = 2 * genrand_real1() * myDA->weight;
    myDA->c = 2 * genrand_real1() * myDA->weight;
    myDA->f = 2 * genrand_real1();
    myDA->e = myDA->weight;
}

/**
 * Function that finds the neighboring solutions
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param NS number of solutions for a data type
 * @param DIM number of dimensions for a data type
 */
void findNeighbors(DA *myDA, initData *myData, int i, int DIM, int NS) {
    int index = 0;
    myDA->numNeighbors = 0;

    for (int k = 0; k < NS; k++) {
        distance(myDA, myData, i, k, DIM);
        if (lessR(myDA, DIM)) {
            index++;
            myDA->numNeighbors++;
            for (int j = 0; j < DIM; ++j) {
                myDA->neighborsPop[index][j] = myData->population[k][j];
                myDA->neighborsStep[index][j] = myDA->step[k][j];
            }
        }
    }
}

/**
 * Function that finds the distances between the dragonflies
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param j the next dragonfly
 * @param DIM number of dimensions for a data type
 */
void distance(DA *myDA, initData *myData, int i, int j, int DIM) {
    for (int k = 0; k < DIM; k++) {
        myDA->o[k] = sqrt(pow((myData->population[i][k] - myData->population[j][k]), 2));
    }
}

/**
 * Function that calculates the separation factor
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void separation(DA *myDA, initData *myData, int DIM, int i) {
    if (myDA->numNeighbors > 1) {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            for (int k = 0; k < DIM; ++k) {
                myDA->sVector[k] += myDA->neighborsPop[j][k] - myData->population[i][k];
            }
        }
        for (int k = 0; k < DIM; ++k) {
            myDA->sVector[k] = -myDA->sVector[k];
        }
    }
}

/**
 * Function that calculates the alignment factor
 * @param myDA the DA structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void alignment(DA *myDA, int DIM, int i) {
    if (myDA->numNeighbors > 1) {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            for (int k = 0; k < DIM; k++) {
                myDA->aVector[k] += myDA->neighborsStep[j][k];
            }
        }
        for (int k = 0; k < DIM; k++) {
            myDA->aVector[k] = myDA->aVector[k] / myDA->numNeighbors;
        }
    } else {
        for (int j = 0; j < DIM; ++j) {
            myDA->aVector[j] = myDA->step[i][j];
        }
    }
}

/**
 * Function that calculates the cohesion factor
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void cohesion(DA *myDA, initData *myData, int DIM, int i) {
    if (myDA->numNeighbors > 1) {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            for (int k = 0; k < DIM; k++) {
                myDA->cVector[k] += myDA->neighborsPop[j][k];
            }
        }
        for (int k = 0; k < DIM; k++) {
            myDA->cVector[k] = myDA->cVector[k] / myDA->numNeighbors;
        }
    } else {
        for (int j = 0; j < DIM; ++j) {
            myDA->cVector[j] = myDA->neighborsPop[i][j];
        }
    }
    for (int j = 0; j < DIM; ++j) {
        myDA->cVector[j] -= myData->population[i][j];
    }
}

/**
 * Function that calculates the distraction factor
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void distraction(DA *myDA, initData *myData, int i, int DIM) {
    distance(myDA, myData, i, myDA->enemyPos, DIM);
    if (lessR2(myDA, DIM)) {
        for (int j = 0; j < DIM; ++j) {
            myDA->eVector[j] = myData->population[myDA->enemyPos][j] + myData->population[i][j];
        }
    }
}

/**
 * Function that calculates the attraction factor
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void attraction(DA *myDA, initData *myData, int i, int DIM) {
    distance(myDA, myData, i, myDA->foodPos, DIM);
    if (lessR2(myDA, DIM)) {
        for (int j = 0; j < DIM; ++j) {
            myDA->fVector[j] = myData->population[myDA->foodPos][j] - myData->population[i][j];
        }
    }
}

/**
 * Function that updates the velocity vector and population
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void updateStepPosition(DA *myDA, initData *myData, int i, int DIM) {
    for (int t = 0; t < DIM; ++t) {
        // velocity matrix
        myDA->step[i][t] = (myDA->s * myDA->sVector[t] + myDA->a * myDA->aVector[t] +
                                myDA->c * myDA->cVector[t] + myDA->f * myDA->fVector[t] +
                                myDA->e * myDA->eVector[t]) + myDA->w * myDA->step[i][t];

        // if the new position is outside the range of
        // the bounds, then make it equal to the bounds
        checkBounds(myData, myDA->step[i][t]);

        // position matrix
        myData->population[i][t] = myData->population[i][t] + myDA->step[i][t];

        // if the new population is outside the range of
        // the bounds, then make it equal to the bounds
        checkBounds(myData,myData->population[i][t]);

    }
}

/**
 * Function that updates the velocity vector and population
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void updateStepPosition2(DA *myDA, initData *myData, int i, int DIM) {
    for (int t = 0; t < DIM; ++t) {
        // velocity matrix
        myDA->step[i][t] = myDA->w * myDA->step[i][t] + genrand_real1() * myDA->sVector[t] + genrand_real1() * myDA->aVector[t] +
                            genrand_real1() * myDA->cVector[t];

        // if the new position is outside the range of
        // the bounds, then make it equal to the bounds
        checkBounds(myData, myDA->step[i][t]);

        // position matrix
        myData->population[i][t] = myData->population[i][t] + myDA->step[i][t];

        // if the new population is outside the range of
        // the bounds, then make it equal to the bounds
        checkBounds(myData,myData->population[i][t]);

    }
}

/**
 * if our vector distance is within the radius of our population
 * then the vector is within the distance of our population
 * @param myDA the DA structure pointer
 * @param DIM number of dimensions for a data type
 */
int lessR(DA *myDA, int DIM) {
    int counter = 0;
    for (int i = 0; i < DIM; ++i) {
        if (myDA->o[i] <= myDA->r && myDA->o[i] != 0) {
            counter++;
        }
    }
    if (counter == DIM-1) {
        return 1;
    }
    return 0;
}

/**
 * if our vector distance is within the radius of our population
 * then the vector is within the distance of our population
 * @param myDA the DA structure pointer
 * @param DIM number of dimensions for a data type
 */
int lessR2(DA *myDA, int DIM) {
    int counter = 0;
    for (int i = 0; i < DIM; ++i) {
        if (myDA->o[i] <= myDA->r) {
            counter++;
        }
    }
    if (counter == DIM-1) {
        return 1;
    }
    return 0;
}

/**
 * if our vector distance is within the radius of our population
 * then the vector is within the distance of our population
 * @param myDA the DA structure pointer
 * @param DIM number of dimensions for a data type
 */
int greaterR3(DA *myDA, int DIM) {
    for (int i = 0; i < DIM; ++i) {
        if (myDA->o[i] > myDA->r) {
            return 1;
        }
    }
    return 0;
}

/**
 * Function that updates the population using Random Walk
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void randomWalk(DA *myDA, initData *myData, int i, int DIM) {
    for (int t = 0; t < DIM; ++t) {
        myData->population[i][t] = myData->population[i][t] + levyFlight(DIM) * myData->population[i][t];
        myDA->step[i][t] = 0;
        // if the new position is outside the range of
        // the bounds, then make it equal to the bounds
        checkBounds(myData,myData->population[i][t] );

    }
}

/**
 * Function that calculates a part of the random walk equation
 * @param DIM number of dimensions for a data type
 */
double levyFlight(int DIM) {
    double beta = 1.5;
    double rand = genrand_real1();
    double rand2 = genrand_real1();
    double R = factorial(DIM - 1);
    double sigma = pow(
            ((R * (1 + beta) * sin((M_PI * beta) / 2)) / R * ((1 + beta) / 2) * beta * powf(2, (beta - 1) / 2)),
            1 / beta);
    double o = 0.01 * (rand * sigma) / pow(fabs(rand2), 1 / beta);
    return o;
}

/**
 * Function that calculates the factorial of a number
 * @param DIM number of dimensions for a data type
 */
int factorial(int DIM) {
    double fact = 1;
    for (int i = 1; i <= DIM; i++) {
        fact = fact * i;
    }
    return DIM;
}
