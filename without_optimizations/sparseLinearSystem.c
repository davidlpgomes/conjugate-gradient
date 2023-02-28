#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sparseLinearSystem.h"

#include "utils.h"

SparseLinearSystem* allocLinearSystem(unsigned int k, unsigned int n) {
    SparseLinearSystem* sls = malloc(sizeof(SparseLinearSystem));

    if (!testAlloc(sls, "Sparse Linear System"))
        exit(1);

    sls->k = k;
    sls->n = n;

    sls->A = malloc(n * sizeof(double*));

    if (!testAlloc(sls->A, "Coefficients A")) {
        freeLinearSystem(sls);
        exit(1);
    }

    sls->sizeA = getSparseMatrixSize(k, n);

    if (!sls->sizeA) {
        fprintf(stderr, "ERROR: Invalid number of diagonals to size n\n");
        freeLinearSystem(sls);
        exit(1);
    }

    sls->A[0] = malloc(sls->sizeA * sizeof(double));

    if (!testAlloc(sls->A[0], "Coefficients A vector")) {
        freeLinearSystem(sls);
        exit(1);
    }

    for (int i = 1; i < sls->n; i++)
        sls->A[i] = sls->A[i - 1] + getLineSize(sls, i - 1);

    sls->b = malloc(n * sizeof(double));

    if (!testAlloc(sls->b, "Solutions b")) {
        freeLinearSystem(sls);
        exit(1);
    }

    return sls;
}

void copyLinearSystem(SparseLinearSystem* source, SparseLinearSystem* destination) {
    if (source->n != destination->n || source->k != destination->k)
        return;

    memcpy(destination->b, source->b, destination->n);
    memcpy(destination->A[0], source->A[0], destination->sizeA);

    return;
}

void fprintLinearSystem(FILE* stream, SparseLinearSystem* sls) {
    for (int i = 0; i < sls->n; i++) {
        fprintf(stream, "\n  ");

        for (int j = 0; j < sls->n; j++)
            fprintf(stream, "%10g", getMatrixValueByPos(sls, i, j));

        fprintf(stream, "   |   %g", sls->b[i]);
    }

    fprintf(stream, "\n\n");

    return;
}

void fprintVector(FILE* stream, double* vector, unsigned int n) {
    for (int i = 0; i < n; i++)
        fprintf(stream, "%.15g ", vector[i]);

    fprintf(stream, "\n");

    return;
}

void freeLinearSystem(SparseLinearSystem* sls) {
    if (!sls)
        return;

    if (sls->A[0])
        free(sls->A[0]);

    if (sls->A)
        free(sls->A);

    if (sls->b)
        free(sls->b);

    free(sls);

    return;
}

int getBeginIndexLine(SparseLinearSystem* sls, int i) {
    return max(0, i - sls->k / 2);
}

int getEndIndexLine(SparseLinearSystem* sls, int i) {
    return min(sls->n - 1, i + sls->k / 2);
}

int getLineSize(SparseLinearSystem* sls, int i) {
    return getEndIndexLine(sls, i) - getBeginIndexLine(sls, i) + 1;
}

void getMatrixTranspose(SparseLinearSystem* matrix, SparseLinearSystem* transpose) {
    for (int i = 0; i < matrix->n; i++)
        for (int j = getBeginIndexLine(matrix, i); j <= getEndIndexLine(matrix, i); j++)
            setMatrixValueByPos(transpose, j, i, getMatrixValueByPos(matrix, i, j));

    return;
}

double getMatrixValueByPos(SparseLinearSystem* sls, int i, int j) {
    if (j < getBeginIndexLine(sls, i) || j > getEndIndexLine(sls, i))
        return 0.0;

    return sls->A[i][j - getBeginIndexLine(sls, i)];
}

void multiplyTranposeWithA(SparseLinearSystem* matrix, SparseLinearSystem* transpose) {
    SparseLinearSystem* temp = allocLinearSystem(2 * matrix->k - 1, matrix->n);
    populateLinearSystem(temp);

    double value;

    for (int i = 0; i < matrix->n; i++) {
        for (int j = getBeginIndexLine(temp, i); j <= getEndIndexLine(temp, i); j++) {
            value = 0;

            for (int k = getBeginIndexLine(matrix, i); k <= getEndIndexLine(matrix, i); k++)
                value += getMatrixValueByPos(transpose, i, k) * getMatrixValueByPos(matrix, k, j);

            setMatrixValueByPos(temp, i, j, value);
        }
    }

    matrix->k = temp->k;
    matrix->sizeA = temp->sizeA;

    free(matrix->A[0]);
    free(matrix->A);
    matrix->A = temp->A;

    free(temp->b);
    free(temp);

    return;
}

void multiplyTranposeWithB(SparseLinearSystem* matrix, SparseLinearSystem* transpose) {
    memcpy(transpose->b, matrix->b, matrix->n * sizeof(double));

    for (int i = 0; i < matrix->n; i++) {
        matrix->b[i] = 0;

        for (int j = getBeginIndexLine(matrix, i); j <= getEndIndexLine(matrix, j); j++)
            matrix->b[i] += getMatrixValueByPos(transpose, i, j) * transpose->b[j];
    }

    return;
}

double normL2(double* x, int n) {
    double norma = 0;

    for (int i = 0; i < n; i++)
        norma += pow(x[i], 2);

    return sqrt(norma);
}

void populateLinearSystem(SparseLinearSystem* sls) {
    for (int i = 0; i < sls->n; i++)
        for (int j = 0; j < sls->n; j++)
            setMatrixValueByPos(sls, i, j, generateRandomA(i, j, sls->k));

    for (int i = 0; i < sls->n; i++)
        sls->b[i] = generateRandomB(sls->k);

    return;
}

void setMatrixValueByPos(SparseLinearSystem* sls, int i, int j, double value) {
    if (j < getBeginIndexLine(sls, i) || j > getEndIndexLine(sls, i))
        return;

    sls->A[i][j - getBeginIndexLine(sls, i)] = value;

    return;
}
