#ifndef DA_495_FUNCTIONS_H
#define DA_495_FUNCTIONS_H
/**
 * @file Functions.h
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
 * calculate the Schwefel optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double schwefel (double *array, int n);

/**
 * @brief
 * calculate the DeJong 1 optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double deJong (double *array, int n);

/**
 * @brief
 * calculate the Rosenbrock's Saddle optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double rosenbrock (double *array, int n);

/**
 * @brief
 * calculate the Rastrigin optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double rastrigin (double *array, int n);

/**
 * @brief
 * calculate the Griewangk optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the product and summation
 */
double griewangk (double *array, int n);

/**
 * @brief
 * calculate the Sine Envelope Sine Wave optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double sineEnvelope (double *array, int n);

/**
 * @brief
 * calculate the Stretch V Sine Wave optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double sineStretched (double *array, int n);

/**
 * @brief
 * calculate the Ackley One optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double ackley1 (double *array, int n);

/**
 * @brief
 * calculate the Ackley Two optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double ackley2 (double *array, int n);

/**
 * @brief
 * calculate the Egg Holder optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double eggHolder (double *array, int n);

/**
 * @brief
 * calculate the Rana optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double rana (double *array, int n);

/**
 * @brief
 * calculate the Pathological optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double pathological (double *array, int n);

/**
 * @brief
 * calculate the Michalewicz optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double michalewicz (double *array, int n);

/**
 * @brief
 * calculate the Master's Cosine Wave optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double masters (double *array, int n);

/**
 * @brief
 * calculate the Quartic optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double quartic (double *array, int n);

/**
 * @brief
 * calculate the Levy optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double levy (double *array, int n);

/**
 * @brief
 * calculate the Step optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double step (double *array, int n);

/**
 * @brief
 * calculate the Alpine optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double alpine (double *array, int n);


#endif //DA_495_FUNCTIONS_H
