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
        myData->fitness = getFun(myData->fitness, myData->population, DIM, NS, myData->functionNumber);
        myDA->step = createDblArray(NS, DIM);
        myDA->step = fillIn(myDA->step, DIM, NS, myData->min, myData->max);

        // initialize the size of the arrays
        myDA->sVector = singleArray(NS);
        myDA->aVector = singleArray(NS);
        myDA->cVector = singleArray(NS);
        myDA->fVector = singleArray(NS);
        myDA->eVector = singleArray(NS);
        myDA->o = singleArray(DIM);

        // start the algorithm
        startDA(myDA, myData, NS, DIM, iterations, fitnessCounter, DAOut);
    }
    // close the file
    fclose(funFile);
    fclose(DAOut);
}

void startDA(DA *myDA, initData *myData, int NS, int DIM, int iterations, int fitnessCounter, FILE *fileOut) {
    // start the clock
    clock_t start;
    start = clock();
    for (int t = 0; t < iterations; ++t) {
        for (int i = 0; i < NS - 1; ++i) {
            // update weights and radius
            updateWeights(myDA, myData, t, iterations);

            // enemy and food value
            myDA->enemy = findWorst(myData->fitness, NS);
            myDA->food = findBest(myData->fitness, NS);
            // enemy and food position
            myDA->foodPos = bestIndex(myData->fitness, NS);
            myDA->enemyPos = worstIndex(myData->fitness, NS);

            // update neighbors
            findNeighbors(myDA, myData, i, DIM, NS);

            // calculate S, A, C, F and E
            calculateVectors(myDA, myData, NS, DIM, i);

            //update velocity and position vectors
            if (myDA->numNeighbors > 1) {
                updateStepPosition(myDA, myData, i, DIM);
            } else {
                // perform random walk
                randomWalk(myDA, myData, i, 0, DIM);
            }

            // update the fitness of our new population
            myData->fitness[i] = getFunSingle(myData->fitness[i], myData->population[i], DIM, myData->functionNumber);
            fitnessCounter++;

            if (myData->fitness[i] < myDA->food) {
                myDA->food = myData->fitness[i];
            }
            if (myData->fitness[i] > myDA->enemy) {
                myDA->enemy = myData->fitness[i];
            }

        }
        printf("(%d)Best: ", t);
        printf("%lf ", myDA->food);
        printf("Worst: ");
        printf("%lf ", myDA->food);
        printf("\n");

        // fprintf(fileOut, "%lf, ", myDA->food);
    }
    start = (((clock() - start)));
    printf("\nExperiment for %d took, %lf, Counter, %d\n",
           myData->functionNumber, ((((double) start) / CLOCKS_PER_SEC) * 1000), fitnessCounter);

//    fprintf(fileOut, "\nExperiment for %d took, %lf, Counter, %d\n",
//            myData->functionNumber, ((((double) start) / CLOCKS_PER_SEC) * 1000), fitnessCounter);
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

void calculateVectors(DA *myDA, initData *myData, int NS, int DIM, int i) {
    separation(myDA, myData, DIM, i);
    alignment(myDA, DIM, i);
    cohesion(myDA, myData, DIM, i);
    attraction(myDA, myData, i, NS, DIM);
    distraction(myDA, myData, i, NS, DIM);
}

void separation(DA *myDA, initData *myData, int DIM, int i) {
    if (myDA->numNeighbors > 1) {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            for (int k = 0; k < DIM; ++k) {
                myDA->sVector[j] += myDA->neighborsPop[j][k] - myData->population[i][k];
            }
            myDA->sVector[j] = -myDA->sVector[j];
        }
    } else {
        for (int j = 0; j < myDA->numNeighbors; ++j) {

            myDA->sVector[j] = 0;
        }
    }

}

void alignment(DA *myDA, int DIM, int i) {
    if (myDA->numNeighbors > 1) {

        for (int j = 0; j < myDA->numNeighbors; ++j) {
            for (int k = 0; k < DIM; ++k) {
                myDA->aVector[j] += myDA->neighborsStep[j][k];
            }
            myDA->aVector[j] = myDA->aVector[j] / myDA->numNeighbors;
        }

    } else {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            myDA->aVector[j] = myDA->step[i][j];
        }
    }
}

void cohesion(DA *myDA, initData *myData, int DIM, int i) {
    if (myDA->numNeighbors > 1) {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            for (int k = 0; k < DIM; ++k) {

                myDA->cVector[j] += myDA->neighborsPop[j][k];
            }
            myDA->cVector[j] = myDA->cVector[j] / myDA->numNeighbors;
            myDA->cVector[j] -= myData->population[i][j];
        }

    } else {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            myDA->cVector[j] = 0;
        }
    }
}

void attraction(DA *myDA, initData *myData, int i, int DIM) {
    distance(myDA->o, myData->population, myData->population, i, myDA->foodPos, DIM);
    if (lessR2(myDA, DIM)) {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            myDA->fVector[j] = myData->population[myDA->foodPos][j] - myData->population[i][j];
        }
    } else {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            myDA->fVector[j] = 0;
        }
    }

}

void distraction(DA *myDA, initData *myData, int i, int DIM) {
    distance(myDA->o, myData->population, myData->population, i, myDA->enemyPos, DIM);
    if (lessR2(myDA, DIM)) {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            myDA->eVector[j] = myData->population[myDA->enemyPos][j] + myData->population[i][j];
        }
    } else {
        for (int j = 0; j < myDA->numNeighbors; ++j) {
            myDA->fVector[j] = 0;
        }
    }
}

void updateStepPosition(DA *myDA, initData *myData, int i, int DIM) {
    for (int t = 0; t < DIM; ++t) {
        // velocity matrix
        myDA->step[i + 1][t] = (myDA->s * myDA->sVector[i] + myDA->a * myDA->aVector[i] +
                                myDA->c * myDA->cVector[i] + myDA->f * myDA->fVector[i] +
                                myDA->e * myDA->eVector[i]) + myDA->w * myDA->step[i][t];

        // if the new position is outside the range of
        // the bounds, then make it equal to the bounds
        if (myDA->step[i + 1][t] > myData->max)
            myDA->step[i + 1][t] = myData->max;
        if (myDA->step[i + 1][t] < myData->min)
            myDA->step[i + 1][t] = myData->min;

        // position matrix
        myData->population[i + 1][t] = myData->population[i][t] + myDA->step[i + 1][t];

        // if the new population is outside the range of
        // the bounds, then make it equal to the bounds
        if (myData->population[i + 1][t] > myData->max)
            myData->population[i + 1][t] = myData->max;
        if (myData->population[i + 1][t] < myData->min)
            myData->population[i + 1][t] = myData->min;
    }
}

void distance(double *returnArr, double **a, double **b, int i, int j, int DIM) {
    for (int k = 0; k < DIM; ++k) {
        returnArr[i] = pow(sqrt(a[i][k] - b[j][k]), 2);
    }
}

// find the neighboring solutions
void findNeighbors(DA *myDA, initData *myData, int i, int DIM, int NS) {
    int index = 0;
    myDA->numNeighbors = 0;
    for (int k = 0; k < NS; ++k) {
        distance(myDA->o, myData->population, myData->population, i, k, DIM);
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

// if our vector distance is within the radius of our population
// then the vector is within the distance of our population
int lessR(DA *myDA, int DIM) {
    int counter = 0;
    for (int i = 0; i < DIM; ++i) {
        if (myDA->o[i] <= myDA->r && myDA->o[i] != 0) {
            counter++;
        }
    }
    if (counter == DIM) {
        return 1;
    }
    return 0;
}

int lessR2(DA *myDA, int DIM) {
    int counter = 0;
    for (int i = 0; i < DIM; ++i) {
        if (myDA->o[i] <= myDA->r) {
            counter++;
        }
    }
    if (counter == DIM) {
        return 1;
    }
    return 0;
}

void randomWalk(DA *myDA, initData *myData, int i, int j, int DIM) {
    for (int t = 0; t < DIM; ++t) {
        myData->population[i + 1][t] = myData->population[i][j] + levyFlight(DIM) * myData->population[i][j];
        myDA->step[i + 1][t] = 0;
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
    for (int i = 0; i <= DIM; i++) {
        DIM = DIM * i;
    }
    return DIM;
}
