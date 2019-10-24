//
// Created by emyli on 10/16/2019.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "DA.h"
#include "Util.h"
#include "ArrayMem.h"
#include "SelectFunctions.h"
#include "mt19937ar.h"
#include <math.h>
#include <rpcndr.h>
#include <stdbool.h>

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

        // fill our population and step vectors with random numbers
        // and get the fitness of our population
        myData->population = fillIn(myData->population, DIM, NS, myData->min, myData->max);
        //myData->fitness = getFun(myData->fitness, myData->population, DIM, NS, myData->functionNumber);
        myDA->step = createDblArray(NS, DIM);
        myDA->step = fillIn(myDA->step, DIM, NS, myData->min, myData->max);

        // initialize the size of the arrays
        myDA->sVector = singleArray(NS);
        myDA->aVector = singleArray(NS);
        myDA->cVector = singleArray(NS);
        myDA->fVector = singleArray(NS);
        myDA->eVector = singleArray(NS);

        // start the algorithm
        startDA(myDA, myData, NS, DIM, iterations, fitnessCounter, DAOut);
    }
    // close the file
    fclose(funFile);
    fclose(DAOut);
}

void startDA(DA *myDA, initData *myData, int NS, int DIM, int iterations, int fitnessCounter, FILE* fileOut) {
    // start the clock
    clock_t start;
    start = clock();
    for (int t = 0; t < iterations; ++t) {
        for (int i = 0; i < NS; ++i) {
                // update weights and radius
                updateWeights(myDA, myData, t, iterations);

                // calculate objective values
                myData->fitness = getFun(myData->fitness, myData->population, DIM, NS, myData->functionNumber);
                myDA->enemy = findWorst(myData->fitness, NS);
                myDA->food = findBest(myData->fitness, NS);

                // calculate S, A, C, F and E
                calculateVectors(myDA, myData, NS, i);

                // euclidean distance - is this used to update the radius ?
                //update velocity and position vectors
                if (NS > 1) { // is this the correct comparison?
                    updateStepPositon(myDA, myData, NS, i, DIM);
                } else {
                    // perform random walk
                    randomWalk(myDA, myData, i, j, DIM);
                }
                // find the fitness of our new population
                myData->fitness[i] = getFunSingle(myData->fitness[i], myData->population[i], DIM,
                                                  myData->functionNumber);
                fitnessCounter++;

                if (myData->fitness[i] < myDA->food) {
                    myDA->food = myData->fitness[i];
                }
                if (myData->fitness[i] > myDA->enemy) {
                    myDA->enemy = myData->fitness[i];
                }

        }
        fprintf(fileOut,"%lf, ", myDA->food);
    }
    start = (((clock() - start)));
    fprintf(fileOut, "\nExperiment for %d took, %lf, Counter, %d\n",
           myData->functionNumber, ((((double) start) / CLOCKS_PER_SEC) * 1000), fitnessCounter);
}

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

void calculateVectors(DA *myDA, initData *myData, int NS, int i) {
    separation(myDA, myData, NS, i);
    alignment(myDA, myData, NS, i);
    cohesion(myDA, myData, NS, i);
    attraction(myDA, myData, i);
    distraction(myDA, myData, i);
}

void separation(DA *myDA, initData *myData, int DIM, int i) {
    double separation = 0;
    if(myDA->numNeighbors >1){
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            for (int k = 0; k < DIM; ++k) {
                separation += myDA->neighborsPop[j][k] - myData->population[i][k];
            }
        }
        myDA->sVector[i] = -separation;
    } else {
        myDA->sVector[i] = 0;
    }

}

void alignment(DA *myDA, initData *myData, int NS, int i) {
    double alignment = 0;
    for (int j = 0; j < NS; ++j) {
        alignment += myData->fitness[j]; //velocity? myData->fitness[i][j];
    }
    myDA->aVector[i] = alignment / NS;
}

void cohesion(DA *myDA, initData *myData, int NS, int i) {
    double cohesion = 0;
    for (int j = 0; j < NS; ++j) {
        cohesion += myData->fitness[j];
    }
    myDA->cVector[i] = cohesion / NS - myData->fitness[i];
}

void attraction(DA *myDA, initData *myData, int i) {
    myDA->fVector[i] = myDA->food - myData->fitness[i];
}

void distraction(DA *myDA, initData *myData, int i) {
    myDA->eVector[i] = myDA->enemy - myData->fitness[i];
}

void updateStepPositon(DA *myDA, initData *myData, int NS, int i, int DIM) {
    for (int t = 0; t < DIM; ++t) {
        // velocity matrix
        myDA->step[i + 1][t] = (myDA->s * myDA->sVector[i] + myDA->a * myDA->aVector[i] +
                                myDA->c * myDA->cVector[i] + myDA->f * myDA->fVector[i] +
                                myDA->e * myDA->eVector[i]) + myDA->w * myDA->step[i][t];

        // if the new position is outside the range of
        // the bounds, then make it equal to the bounds
        if (myDA->step[i][t] > myData->max)
            myDA->step[i][t] = myData->max;
        if (myDA->step[i][t] < myData->min)
            myDA->step[i][t] = myData->min;

        // position matrix
        myData->population[i + 1][t] = myData->population[i][t] + myDA->step[i + 1][t];
    }
}

void randomWalk(DA *myDA, initData *myData, int i, int j, int DIM) {
    for (int t = 0; t < DIM; ++t) {
        myData->population[i + 1][t] = myData->population[i][j] + levyFlight(DIM) * myData->population[i][j];
        // if the new position is outside the range of
        // the bounds, then make it equal to the bounds
        if (myData->population[i][t] > myData->max)
            myData->population[i][t] = myData->max;
        if (myData->population[i][t] < myData->min)
            myData->population[i][t] = myData->min;
    }
}

double levyFlight(int DIM) {
    double beta = 1.5;
    double rand = genrand_real1();
    double rand2 = genrand_real1();
    double R = factorial(DIM - 1);
    double sigma = pow(
            ((R * (1 + beta) * sinf((M_PI * beta) / 2)) / R * ((1 + beta) / 2) * beta * powf(2, (beta - 1) / 2)),
            1 / beta);
    double o = 0.01 * (rand * sigma) / pow(fabs(rand2), 1 / beta);
    return o;
}

int factorial(int DIM) {
    for (int i = 1; i <= DIM; i++) {
        DIM = DIM * i;
    }
    return DIM;
}

double distance(DA *myDA, initData *myData, int i, int j, int DIM){
    for (int k = 0; k < DIM; ++k) {
        myDA->o[i] = powf(sqrtf(myData->population[i][k] - myData->population[j][k]),2);
    }
}

// find the neighboring solutions
double findNeighbors(DA *myDA, initData *myData, int i, int DIM, int NS){
    int index = 0;
    myDA->numNeighbors = 0;
    for (int k = 0; k < NS; ++k) {
        distance(myDA,myData,i,k,DIM);
        if(lessR(myDA,DIM)){
            index++;
            myDA->numNeighbors++;
        for (int j = 0; j < DIM; ++j) {
                myDA->neighborsPop[index][j] = myData->population[k][j];
                myDA->neighborsStep[index][j] = myDA->step[k][j];
            }
        }
    }
}

// if our vector distance is within the radius of our population
// then the vector is within the distance of our population
boolean lessR(DA *myDA, int DIM){
    int counter = 0;
    for (int i = 0; i < DIM; ++i) {
        if (myDA->o[i] <= myDA->r && myDA->o[i] != 0){
            counter++;
        }
    }
    if (counter == DIM){
        return true;
    }
    return false;
}