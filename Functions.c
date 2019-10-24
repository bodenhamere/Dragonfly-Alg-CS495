/**
 * @file Functions.c
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
#include <math.h>

/**
 * calculate the Schwefel optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double schwefel(double *array, int n) {
    double sum = 0.0;

    for(int i = 0; i < n; i++) {
        sum += (array[i] * -1) * sin(sqrt(fabs(array[i])));
    }

    return sum = (418.9829 * n) - sum;
}

/**
 * calculate the DeJong 1 optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double deJong(double *array, int n) {
    double sum = 0.0;

    for(int i = 0; i < n; i++) {
        sum += pow(array[i], 2);
    }

    return sum;
}

/**
 * calculate the Rosenbrock's Saddle optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double rosenbrock(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        sum += (100 * pow((pow(array[i], 2) - array[i+1]), 2)) + pow((1 - array[i]), 2);
    }

    return sum;
}

/**
 * calculate the Rastrigin optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double rastrigin(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n; i++) {
        sum += pow(array[i], 2) - (10 * cos(2 * M_PI * array[i]));
    }

    return 10 * n * sum;
}

/**
 * calculate the Griewangk optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the product and summation
 */
double griewangk(double *array, int n) {
    double sum = 0.0;
    double prod = 1;

    for(int i = 1; i < n; i++) {
        sum += (pow(array[i], 2) / 4000);
        prod *= cos(array[i] / sqrt(i));
    }

    return 1 + sum - prod;
}

/**
 * calculate the Sine Envelope Sine Wave optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double sineEnvelope(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        sum += 0.5 + (sin(pow(pow(array[i], 2) + pow(array[i+1], 2) - 0.5, 2))/
                      pow(1 + (0.001 * (pow(array[i], 2) + pow(array[i+1], 2))), 2));
    }

    return -1 * sum;
}

/**
 * calculate the Stretch V Sine Wave optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double sineStretched(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        sum += pow(pow(array[i], 2) + pow(array[i+1], 2), 0.25) *
               sin(pow(50 * pow(pow(array[i], 2) + pow(array[i+1], 2), 0.1), 2)) + 1;
    }

    return sum;
}

/**
 * calculate the Ackley One optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double ackley1(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        sum += ((1 / pow(M_E, 0.2)) * sqrt(pow(array[i], 2) + pow(array[i+1], 2))) +
               (3 * (cos(2 * array[i]) + sin(2 * array[i+1])));
    }

    return sum;
}

/**
 * calculate the Ackley Two optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double ackley2(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        sum += 20 + M_E - (20 / pow(M_E, 0.2 * sqrt((pow(array[i], 2) + pow(array[i+1], 2)) / 2))) -
               pow(M_E, 0.5 * (cos(2 * M_PI * array[i]) + cos(2 * M_PI * array[i+1])));
    }

    return sum;
}

/**
 * calculate the Egg Holder optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double eggHolder(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        sum += (-1 * array[i] * sin(sqrt(fabs(array[i] - array[i+1] - 47)))) -
               ((array[i+1] + 47) * sin(sqrt(fabs(array[i+1] + 47 + (array[i] / 2)))));
    }

    return sum;
}

/**
 * calculate the Rana optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double rana(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        sum += (array[i] * sin(sqrt(fabs(array[i+1] - array[i] + 1))) * cos(sqrt(fabs(array[i+1] + array[i] + 1)))) +
               ((array[i+1] + 1) * cos(sqrt(fabs(array[i+1] - array[i] + 1))) * sin(array[i+1] + array[i] + 1));
    }

    return sum;
}

/**
 * calculate the Pathological optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double pathological(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        sum += 0.5 + ((sin(pow(sqrt((100 * pow(array[i], 2)) + pow(array[i+1], 2)), 2)) - 0.5) /
                      (1 + pow((.001 * (pow(array[i], 2) - (2 * array[i] * array[i+1]) + pow(array[i+1], 2))), 2)));
    }

    return sum;
}

/**
 * calculate the Michalewicz optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double michalewicz(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n-1; i++) {
        sum+= sin(array[i]) * pow(sin((i * pow(array[i], 2)) / M_PI), 20);
    }

    return -1 * sum;
}


/**
 * calculate the Master's Cosine Wave optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double masters(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        sum += pow(M_E, -0.125 * (pow(array[i], 2) + pow(array[i+1], 2) + (0.5 * array[i+1]* array[i]))) *
               cos(4 * sqrt(pow(array[i], 2) + pow(array[i+1], 2) + (.5 * array[i] * array[i+1])));
    }

    return -1 * sum;
}

/**
 * calculate the Quartic optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double quartic(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n; i++) {
        sum += i * pow(array[i], 4);
    }

    return sum;
}

/**
 * calculate the Levy optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double levy(double *array, int n) {
    double sum = 0.0;

    for(int i = 1; i < n - 1; i++) {
        double wi = 1 + ((array[i] - 1) / 4);
        double wn = 1 + ((array[n] - 1) / 4);

        sum += (pow((wi - 1), 2) * (1 + (10 * pow(sin((M_PI * wi) + 1), 2)))) +
               (pow(wn - 1, 2) * (1 + pow(sin(2 * M_PI * wn), 2)));
    }

    return pow(sin(M_PI * (1 + ((array[1] - 1) / 4))), 2) + sum;
}

/**
 * calculate the Step optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double step(double *array, int n) {
    double sum = 0.0;

    for(int i = 0; i < n - 1; i++) {
        sum += pow(fabs(array[i]) + 0.5, 2);
    }

    return sum;
}

/**
 * calculate the Alpine optimization function
 * based on random number inputs
 * \param array single array
 * \param n the size of the rows of the array
 * \return final calculated number from the summation
 */
double alpine(double *array, int n) {
    double sum = 0.0;

    for(int i = 0; i < n - 1; i++) {
        sum += fabs((array[i] * sin(array[i])) + (0.1 * array[i]));
    }

    return sum;
}
