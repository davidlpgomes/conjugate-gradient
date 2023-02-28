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
    double* sums = calloc(sls->n, sizeof(double));
    int beg, end, lineSize;

    for (int i = 0; i < sls->n; i++) {
        beg = getBeginIndexLine(sls, i);
        end = getEndIndexLine(sls, i);
        lineSize = end - beg + 1;

        for (int j = 0; j < lineSize - lineSize % 4; j += 4) {
            sums[j + beg] += p[i] * sls->A[i][j];
            sums[j + beg + 1] += p[i] * sls->A[i][j + 1];
            sums[j + beg + 2] += p[i] * sls->A[i][j + 2];
            sums[j + beg + 3] += p[i] * sls->A[i][j + 3];
        }

        for (int j = lineSize - lineSize % 4; j < lineSize; ++j)
            sums[j + beg] += p[i] * sls->A[i][j];
    }

    double divisor = multiplyVectors(sums, p, sls->n);
    free(sums);

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
    for (int i = 0; i < sls->n - sls->n % 4; i += 4) {
        p[i] = x[i] + beta * p[i];
        p[i + 1] = x[i + 1] + beta * p[i + 1];
        p[i + 2] = x[i + 2] + beta * p[i + 2];
        p[i + 3] = x[i + 3] + beta * p[i + 3];
    }

    for (int i = sls->n - sls->n % 4; i < sls->n; ++i)
        p[i] = x[i] + beta * p[i];

    return;
}

void calculateNextResidue(SparseLinearSystem* sls, double* r, double* r1, double* p, double alpha) {
    double* temp = calloc(sls->n, sizeof(double));
    int beg, end, lineSize;

    for (int i = 0; i < sls->n; i++) {
        beg = getBeginIndexLine(sls, i);
        end = getEndIndexLine(sls, i);
        lineSize = end - beg + 1;

        for (int j = 0; j < lineSize - lineSize % 4; j += 4) {
            temp[j + beg] += p[i] * sls->A[i][j];
            temp[j + beg + 1] += p[i] * sls->A[i][j + 1];
            temp[j + beg + 2] += p[i] * sls->A[i][j + 2];
            temp[j + beg + 3] += p[i] * sls->A[i][j + 3];
        }

        for (int j = lineSize - lineSize % 4; j < lineSize; ++j)
            temp[j + beg] += p[i] * sls->A[i][j];
    }

    for (int i = 0; i < sls->n - sls->n % 4; i += 4) {
        r1[i] = r[i] - alpha * temp[i];
        r1[i + 1] = r[i + 1] - alpha * temp[i + 1];
        r1[i + 2] = r[i + 2] - alpha * temp[i + 2];
        r1[i + 3] = r[i + 3] - alpha * temp[i + 3];
    }

    for (int i = sls->n - sls->n % 4; i < sls->n; ++i)
        r1[i] = r[i] - alpha * temp[i];

    free(temp);

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
        if (err > maxErr) maxErr = err;
    }

    return maxErr;
}

void calculateResidue(SparseLinearSystem* sls, double* x, double* r) {
    int begI, endI, size;
    double ri = 0;

    for (int i = 0; i != sls->halfK; ++i) {
        endI = i + sls->halfK;

        ri = 0;

        for (int j = 0; j <= endI; j++)
            ri += sls->A[i][j] * x[j];

        r[i] = sls->b[i] - ri;
    }

    for (int i = sls->halfK; i != sls->n - sls->halfK - 1; i++) {
        begI = i - sls->halfK;
        endI = i + sls->halfK;
        size = endI - begI + 1;

        ri = 0;

        for (int j = 0; j < size - size % 4; j += 4) {
            ri += sls->A[i][j] * x[j + begI];
            ri += sls->A[i][j + 1] * x[j + 1 + begI];
            ri += sls->A[i][j + 2] * x[j + 1 + begI];
            ri += sls->A[i][j + 3] * x[j + 2 + begI];
        }

        for (int j = size - size % 4; j < size; j++)
            ri += sls->A[i][j] * x[j + begI];

        r[i] = sls->b[i] - ri;
    }

    for (int i = sls->n - sls->halfK - 1; i != sls->n; i++) {
        begI = i - sls->halfK;

        ri = 0;

        for (int j = 0; j <= sls->n - 1 - begI; j++)
            ri += sls->A[i][j] * x[j + begI];

        r[i] = sls->b[i] - ri;
    }

    return;
}

void calculateZ(double* z, double* m, double* r, int n) {
    for (int i = 0; i < n - n % 4; i += 4) {
        z[i] = m[i] * r[i];
        z[i + 1] = m[i + 1] * r[i + 1];
        z[i + 2] = m[i + 2] * r[i + 2];
        z[i + 3] = m[i + 3] * r[i + 3];
    }

    for (int i = n - n % 4; i < n; ++i)
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
    LIKWID_MARKER_START("OP1");
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
    LIKWID_MARKER_STOP("OP1");
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
