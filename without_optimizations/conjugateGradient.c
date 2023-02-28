#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <likwid.h>

#include "conjugateGradient.h"
#include "utils.h"

double calculateAlpha(SparseLinearSystem* sls, double* r, double* p) {
    double divisor = 0, sum = 0;

    for (int i = 0; i < sls->n; i++) {
        sum = 0;

        for (int j = 0; j < sls->n; j++)
            sum += p[j] * getMatrixValueByPos(sls, j, i);

        divisor += sum * p[i];
    }

    return sumOfSquares(r, sls->n) / divisor;
}

double calculateAlphaJacobi(SparseLinearSystem* sls, double* r, double* p, double* z) {
    double divisor = 0, sum = 0;

    for (int i = 0; i < sls->n; i++) {
        sum = 0;

        for (int j = 0; j < sls->n; j++)
            sum += p[j] * getMatrixValueByPos(sls, j, i);

        divisor += sum * p[i];
    }

    if (!divisor) return 0.0;

    return multiplyVectors(r, z, sls->n) / divisor;
}

double calculateBeta(SparseLinearSystem* sls, double* r, double* r1) {
    double div = sumOfSquares(r, sls->n);

    if (!div) return 0.0;

    return sumOfSquares(r1, sls->n) / div;
}

double calculateBetaJacobi(SparseLinearSystem* sls, double* r, double* r1, double* z, double* z1) {
    double div = multiplyVectors(r, z, sls->n);

    if (!div) return 0.0;

    return multiplyVectors(r1, z1, sls->n) / div;
}

void calculateM(SparseLinearSystem* sls, double* m) {
    for (int i = 0; i < sls->n; i++)
        m[i] = (double) 1 / getMatrixValueByPos(sls, i, i);

    return;
}

void calculateNextP(SparseLinearSystem* sls, double* p, double* x, double beta) {
    for (int i = 0; i < sls->n; i++)
        p[i] = x[i] + beta * p[i];

    return;
}

void calculateNextResidue(SparseLinearSystem* sls, double* r, double* r1, double* p, double alpha) {
    for (int i = 0; i < sls->n; i++) {
        r1[i] = 0;

        for (int j = 0; j < sls->n; j++)
            r1[i] += p[j] * getMatrixValueByPos(sls, j, i);

        r1[i] = r[i] - alpha * r1[i];
    }

    return;
}

double calculateNextX(SparseLinearSystem* sls, double* x, double* p, double alpha) {
    double oldX = x[0];
    x[0] += p[0] * alpha;

    double maxErr = fabs(x[0] - oldX) / fabs(x[0]), err;

    for (int i = 1; i < sls->n; i++) {
        oldX = x[i];
        x[i] += p[i] * alpha;

        err = fabs(x[i] - oldX) / fabs(x[i]);

        if (err > maxErr)
            err = maxErr;
    }

    return maxErr;
}

void calculateResidue(SparseLinearSystem* sls, double* x, double* r) {
    for (int i = 0; i < sls->n; i++) {
        r[i] = sls->b[i];

        int begI = getBeginIndexLine(sls, i);
        int endI = getEndIndexLine(sls, i);

        for (int j = begI; j <= endI; j++)
            r[i] -= sls->A[i][j - begI] * x[j];
    }

    return;
}

void calculateZ(double* z, double* m, double* r, int n) {
    for (int i = 0; i < n; i++)
        z[i] = m[i] * r[i];

    return;
}

void conjugateGradient(SparseLinearSystem* sls, double* x, FILE* stream, int maxItr, double err) {
    double start, timeRes, timeSum = 0.0;

    double* r = malloc(sls->n * sizeof(double));

    if (!testAlloc(r, "residue vector"))
        exit(1);

    calculateResidue(sls, x, r);

    double* r1 = malloc(sls->n * sizeof(double));

    if (!testAlloc(r1, "next residue vector"))
        exit(1);

    double* p = malloc(sls->n * sizeof(double));

    if (!testAlloc(p, "p vector"))
        exit(1);

    memcpy(p, r, sls->n * sizeof(double));

    double alpha, beta;
    double relativeErr = err + 1;

    int iterations = 0;
    while (iterations < maxItr && relativeErr > err) {
        start = timestamp();
        alpha = calculateAlpha(sls, r, p);

        relativeErr = calculateNextX(sls, x, p, alpha);

        calculateNextResidue(sls, r, r1, p, alpha);

        beta = calculateBeta(sls, r, r1);

        calculateNextP(sls, p, r1, beta);

        memcpy(r, r1, sls->n * sizeof(double));

        iterations++;

        timeSum += timestamp() - start;
        fprintf(stream, "# iter %d: %e\n", iterations, relativeErr);
    }

    start = timestamp();
    double resNorma = normL2(r, sls->n);
    timeRes = timestamp() - start;

    fprintf(stream, "# residuo: %.15g\n", resNorma);

    double itrTime = iterations ? (double) timeSum / iterations : 0;
    fprintf(stream, "# Tempo iter: %.15g\n", itrTime);
    fprintf(stream, "# Tempo residuo: %.15g\n#\n", timeRes);

    free(r);
    free(r1);
    free(p);

    return;
}

void conjugateGradientJacobi(SparseLinearSystem* sls, double* x, FILE* stream, int maxItr, double err) {
    double start, timePC, timeRes, timeSum = 0.0;

    double* m = malloc(sls->n * sizeof(double));

    if (!testAlloc(m, "M vector"))
        exit(1);

    start = timestamp();
    calculateM(sls, m);
    timePC = timestamp() - start;

    double* r = malloc(sls->n * sizeof(double));

    if (!testAlloc(r, "residue vector"))
        exit(1);

    double* r1 = malloc(sls->n * sizeof(double));

    if (!testAlloc(r1, "next residue vector"))
        exit(1);

    calculateResidue(sls, x, r);

    double* z = malloc(sls->n * sizeof(double));

    if (!testAlloc(z, "z vector"))
        exit(1);

    double* z1 = malloc(sls->n * sizeof(double));

    if (!testAlloc(z1, "next z vector"))
        exit(1);

    calculateZ(z, m, r, sls->n);

    double* p = malloc(sls->n * sizeof(double));

    if (!testAlloc(p, "p vector"))
        exit(1);

    memcpy(p, z, sls->n * sizeof(double));

    double alpha, beta;
    double relativeErr = err + 1;

    int iterations = 0;
    double startLikwid = timestamp();
    LIKWID_MARKER_START("OPT");
    while (iterations < maxItr && relativeErr > err) {
        start = timestamp();

        alpha = calculateAlphaJacobi(sls, r, p, z);

        relativeErr = calculateNextX(sls, x, p, alpha);

        calculateNextResidue(sls, r, r1, p, alpha);

        calculateZ(z1, m, r1, sls->n);

        beta = calculateBetaJacobi(sls, r, r1, z, z1);

        calculateNextP(sls, p, z1, beta);

        memcpy(r, r1, sls->n * sizeof(double));
        memcpy(z, z1, sls->n * sizeof(double));

        iterations++;

        timeSum += timestamp() - start;
        fprintf(stream, "# iter %d: %e\n", iterations, relativeErr);
    }
    LIKWID_MARKER_STOP("OPT");
    printf("TIME-OP1:%.15g\n", (timestamp() - startLikwid) / (iterations + 1));

    start = timestamp();
    double resNorma = normL2(r, sls->n);
    timeRes = timestamp() - start;

    fprintf(stream, "# residuo: %.15g\n", resNorma);
    fprintf(stream, "# Tempo PC: %.15g\n", timePC);

    double itrTime = iterations ? (double) timeSum / iterations : 0;
    fprintf(stream, "# Tempo iter: %.15g\n", itrTime);
    fprintf(stream, "# Tempo residuo: %.15g\n#\n", timeRes);

    free(r);
    free(r1);
    free(z);
    free(z1);
    free(m);
    free(p);

    return;
}
