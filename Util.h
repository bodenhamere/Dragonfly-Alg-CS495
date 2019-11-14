//
// Created by emyli on 10/16/2019.
//

#ifndef DA_495_UTIL_H
#define DA_495_UTIL_H
/**
 * @file Util.h
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

/**
 * @brief
 *
 * Structure for the bounds, population, fitness and function number
 */
typedef struct _initData1 {
    double max; //! upper bound
    double min; //! lower bound
    double **population; //! random matrix
    double *fitness; //! fitness from matrix
    int functionNumber; //! number that corresponds to one of the 18 functions
} initData;

/**
 * sort the population in ascending order
 * @param myData the data structure
 * @param NS the number of solutions/size of a vector
 * @param DIM the size of the dimensions/size of a vector
 * @param b temporary vector to hold elements
 */
void sortAscendingOrder(initData *myData, int NS, int DIM, double *b);

/**
 * replace the elements of one array into another
 * \param a an array
 * \param b an array
 * \param NS the number of solutions/size of a vector
 * \return a array with array b's elements implemented
 */
double *replaceArray(double *a, const double *b, int NS);

/**
 * find the minimal value in a vector
 * \param a an array
 * \param NS the number of solutions/size of a vector
 */
void findBest(double* bestArr, const double *a,int NS);

/**
 * find the minimal value in a vector
 * \param a an array
 * \param NS the number of solutions/size of a vector
 */
void findWorst(double* worstArr, const double *a, int NS);


/**
 * copy a arrays contents to another
 * \param a a double array
 * \param b a double array
 * \param DIM the number of dimensions/size of a vector
 * \param NS the number of solutions/size of a vector
 * \return a double array with double array b's contents placed into
 */
double ** copyDbl(double **a, double **b, int NS, int DIM);

/**
 * copy a arrays contents to another
 * \param a an array
 * \param b an array
 * \param DIM the number of dimensions/size of a vector
 * \return a array with array b's contents placed into
 */
double * copySingle(double *a, const double *b, int DIM);

/**
 * print out a matrix with dimensions
 * printed out horizontally to a file
 * \param file file to print to
 * \param a an array
 * \param DIM the number of dimensions/size of a vector
 * \param NS the number of solutions/size of a vector
 */
void printDblDim(FILE *file, double **a, int DIM, int NS);

/**
 * print out a vector to a file
 * \param file a file to print to
 * \param a an array
 * \param DIM the number of dimensions/size of a vector
 */
void printSingle(FILE *file, double *a, int DIM);

#endif //DA_495_UTIL_H
