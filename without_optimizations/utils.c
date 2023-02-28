#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "utils.h"


double generateRandomA(unsigned int i, unsigned int j, unsigned int k) {
    double invRandMax = 1.0 / (double) RAND_MAX;

    return ((i == j) ? (double) (k << 1) : 1.0) * (double) rand() * invRandMax;
}

double generateRandomB(unsigned int k) {
    double invRandMax = 1.0 / (double) RAND_MAX;

    return (double) (k << 2) * (double) rand() * invRandMax;
}

int getSparseMatrixSize(unsigned int k, unsigned int n) {
    if (k >= 2 * n)
        return 0;

    int p = floor(k / 2) - 1;

    return n + 2 * (n * (p + 1) - (p * (p + 1) / 2) - p - 1);
}

int max(int a, int b) {
    if (a > b)
        return a;

    return b;
}

int min(int a, int b) {
    if (a < b)
        return a;

    return b;
}

double multiplyVectors(double* v1, double* v2, int n) {
    double sum = 0;

    for (int i = 0; i < n; i++)
        sum += v1[i] * v2[i];

    return sum;
}

double sumOfSquares(double* x, int n) {
    double sum = 0;

    for (int i = 0; i < n; i++)
        sum += pow(x[i], 2);

    return sum;
}

int testAlloc(void* target, char* name) {
    if (target)
        return 1;

    fprintf(stderr, "ERROR: Could not alloc %s!\n", name);

    return 0;
}


double timestamp() {
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC_RAW, &tp);

    return (double) (tp.tv_sec + tp.tv_nsec * 1.0e-9);
}
