#ifndef DA_495_ARRAYMEM_H
#define DA_495_ARRAYMEM_H
/**
 * @file ArrayMem.h
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

/**
 * @brief
 * Generates a double pointer array that has pointer to a single pointer array
 * The arrays are initialized to 0 by using calloc()
 * \param col the size of the columns of the array
 * \param row the size of the rows of the array
 * \return double pointer array
 */
double **createDblArray(int col, int row);

/**
 * Generates a single pointer array
 * The array are initialized to 0 by using calloc()
 * \param n the size of the rows of the array
 * \return single pointer array
 */
double * singleArray(int n);

/**
 * Fills in a double pointer array with the MersenneTwister random numbers from a specified range
 * \param arr double pointer array
 * \param row the size of the rows of the array
 * \param col the size of the columns of the array
 * \param min the minimum size of the random numbers
 * \param max the maximum size of the random numbers
 * \return double pointer array
 */
double **fillIn(double **arr, int row, int col, double min, double max);

/**
 * Free the memory used for a double pointer array
 * \param matrix double pointer array
 * \param row the size of the rows of the array
 */
void freeMem(int row, double **matrix);

#endif //DA_495_ARRAYMEM_H
