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
    double *worstArr; //!< best solution
    double *bestArr; //! worst solution
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

/**
 * @brief
 * Function that reads from the DA input file to
 * initialize the structure variables for the DA algorithm
 * @param iterations number of iterations to run the algorithm for
 * @param fitnessCounter count how many times the fitness function was called
 * @param myData Data structure pointer for initializing the bounds and which function to run
 * @param NS number of solutions for a data type
 * @param DIM number of dimensions for a data type
 */
void readInput(initData *myData, int NS, int DIM, int iterations, int fitnessCounter);

/**
 * @brief
 * Function that runs the DA algorithm
 * @param fileOut output file for the fitness of the population
 * @param NS number of solutions for a data type
 * @param DIM number of dimensions for a data type
 * @param iterations number of iterations to run the algorithm for
 * @param myData Data structure pointer for initializing the bounds and which function to run
 * @param fitnessCounter count how many times the fitness function was called
 * @param myDA the DA structure
*/
void startDA(DA *myDA, initData *myData, int NS, int DIM, int iterations, int fitnessCounter, FILE *fileOut);

/**
 * @brief
 * Function that updates the weights for the DA algorithm
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param iter number of iterations for the algorithm
 * @param maxIter max number of iterations for the algorithm
 */
void updateWeights(DA *myDA, initData *myData, int iter, int maxIter);

/**
 * Function that calculates the separation factor
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void separation(DA *myDA, initData *myData, int NS, int t);

/**
 * Function that calculates the alignment factor
 * @param myDA the DA structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void alignment(DA *myDA, int NS, int t);

/**
 * Function that calculates the cohesion factor
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void cohesion(DA *myDA, initData *myData, int NS, int t);

/**
 * Function that calculates the attraction factor
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void attraction(DA *myDA, initData *myData, int i, int DIM);

/**
 * Function that calculates the distraction factor
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void distraction(DA *myDA, initData *myData, int i, int DIM);

/**
 * Function that updates the velocity vector and population
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void updateStepPosition(DA *myDA, initData *myData, int i, int DIM);

/**
 * Function that updates the velocity vector and population
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void updateStepPosition2(DA *myDA, initData *myData, int i, int DIM);

/**
 * Function that updates the population using Random Walk
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param DIM number of dimensions for a data type
 */
void randomWalk(DA *myDA, initData *myData, int i, int DIM);

/**
 * Function that calculates a part of the random walk equation
 * @param DIM number of dimensions for a data type
 */
double levyFlight(int DIM);

/**
 * Function that calculates the factorial of a number
 * @param DIM number of dimensions for a data type
 */
int factorial(int DIM);

/**
 * if our vector distance is within the radius of our population
 * then the vector is within the distance of our population
 * @param myDA the DA structure pointer
 * @param DIM number of dimensions for a data type
 */
int lessR(DA *myDA, int DIM);

/**
 * if our vector distance is within the radius of our population
 * then the vector is within the distance of our population
 * @param myDA the DA structure pointer
 * @param DIM number of dimensions for a data type
 */
int lessR2(DA *myDA, int DIM);

/**
 * if our vector distance is within the radius of our population
 * then the vector is within the distance of our population
 * @param myDA the DA structure pointer
 * @param DIM number of dimensions for a data type
 */
int greaterR3(DA *myDA, int DIM);

/**
 * Function that finds the distances between the dragonflies
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param j the next dragonfly
 * @param DIM number of dimensions for a data type
 */
void distance(DA *myDA, initData *myData, int i, int j, int DIM);


/**
 * @brief
 * Function that finds the neighboring solutions
 * @param myDA the DA structure pointer
 * @param myData Data structure pointer
 * @param i the current dragonfly
 * @param NS number of solutions for a data type
 * @param DIM number of dimensions for a data type
 */
void findNeighbors(DA *myDA, initData *myData, int i, int DIM, int NS);

#endif //DA_495_DA_H
