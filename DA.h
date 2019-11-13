//
// Created by emyli on 10/16/2019.
//

#ifndef DA_495_DA_H
#define DA_495_DA_H

#include "Util.h"

/**
 * @brief
 *
 * Structure for the DA algorithm
 */
typedef struct _DA1 {
    double s; //!< separation weight (0.1)
    double a; //!< alignment weight (0.1)
    double c; //!< cohesion weight (0.7)
    double f; //!< food source weight (1)
    double e; //!< enemy weight (1)
    double w; //!< inertia weight (0.9-0.2)
    double r; //!< radius
    double weight; //!< a weight
    double food; //!< best solution
    double enemy; //! worst solution
    int foodPos; //!< best solution
    int enemyPos; //! worst solution
    double **step; //!< direction of the movement of the dragonflies
    double *sVector; //!< Separation array
    double *aVector; //!< Alignment array
    double *cVector; //!< Cohesion array
    double *fVector; //!< food source array
    double *eVector; //!< enemy array
    double *o; //!< distance array
    int numNeighbors; //!< number of neighbors within the population
    double **neighborsStep; //!< position of the step neighbors
    double **neighborsPop;//!< position of the population neighbors
} DA;

void readInput(initData *myData, int NS, int DIM, int iterations, int fitnessCounter);

void startDA(DA *myDA, initData *myData, int NS, int DIM, int iterations, int fitnessCounter, FILE *fileOut);

void updateWeights(DA *myDA, initData *myData, int iter, int maxIter);

void separation(DA *myDA, initData *myData, int NS, int t);

void alignment(DA *myDA, int NS, int t);

void cohesion(DA *myDA, initData *myData, int NS, int t);

void attraction(DA *myDA, initData *myData, int i, int DIM);

void distraction(DA *myDA, initData *myData, int i, int DIM);

void updateStepPosition(DA *myDA, initData *myData, int i, int DIM);

void randomWalk(DA *myDA, initData *myData, int i, int j, int DIM);

double levyFlight(int DIM);

int factorial(int DIM);

int lessR(DA *myDA, int DIM);

int lessR2(DA *myDA, int DIM);

void distance(DA *myDA, initData *myData, int i, int j, int DIM);

void findNeighbors(DA *myDA, initData *myData, int i, int DIM, int NS);

#endif //DA_495_DA_H
