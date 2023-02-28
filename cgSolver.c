#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <likwid.h>

#include "conjugateGradient.h"
#include "sparseLinearSystem.h"
#include "utils.h"

#define STR_LEN 1024

#define MIN_OPTS 5
#define MIN_N 11


/*!
    \brief Creates the SparseLinearSystem and solution's vector and runs the conjugate gradient method

    \param n dimension of the linear system
    \param k number of diagonals on the coefficients matrix
    \param i maximum iterations to stop the method
    \param e maximum approximate error to stop the method
    \param file stream to print informations
    \param function the conjugate gradient function to be used
*/
void runConjugateGradient(
    int n, int k, int i, double e, FILE* file, 
    void (function)(SparseLinearSystem*, double*, FILE*, int, double)
) {
    // Create Sparse Linear System
    SparseLinearSystem* sls = allocLinearSystem(k, n);
    populateLinearSystem(sls);

    // Create a copy of the original SLS
    SparseLinearSystem* originalSls = allocLinearSystem(k, n);
    copyLinearSystem(sls, originalSls);

    // Create transpose and multiply to the SLS
    SparseLinearSystem* transpose = allocLinearSystem(k, n);
    populateLinearSystem(transpose);

    getMatrixTranspose(sls, transpose);
    multiplyTranposeWithA(sls, transpose);
    multiplyTranposeWithB(sls, transpose);

    freeLinearSystem(transpose);

    // Create solutions vector
    double* x = calloc(n, sizeof(double));

    if (!testAlloc(x, "solutions vector"))
        exit(1);

    // Run conjugate gradient method
    function(sls, x, file, i, e);

    // Calculates the residue with original Linear System
    double* r = calloc(n, sizeof(double));
    if (!testAlloc(r, "residue vector")) exit(1);

    double start = timestamp();
    LIKWID_MARKER_START("OP2");
    calculateResidue(originalSls, x, r);
    LIKWID_MARKER_STOP("OP2");
    printf("TIME-OP2:%.15g\n", timestamp() - start);

    fprintf(file, "%d\n", sls->n);
    fprintVector(file, x, sls->n);

    // Free allocs
    free(r);
    free(x);
    freeLinearSystem(originalSls);
    freeLinearSystem(sls);

    return;
}

int main(int argc, char* argv[]) {
    /*
        Autor: David Lucas Pereira Gomes
        GRR: 20211757
        Curso: Ciência da Computação
    */
    LIKWID_MARKER_INIT;

    // Configures the random seed
    srand(20222);

    // Initiates the needed parameters
    int n = -1, k = -1, p = -1, i = -1;
    double e = -1;

    char* o = malloc(STR_LEN * sizeof(char));

    if (!testAlloc(o, "output file string"))
        return 1;

    int opt, optSum = 0;

    while ((opt = getopt(argc, argv, ":n:k:p:i:e:o:")) != -1) {
        if (opt == 'n') {
            n = atoi(optarg);

            if (n < MIN_N) {
                fprintf(stderr, "Dimensão da SL deve ser igual ou maior a %d!\n", MIN_N);
                return 1;
            }

            optSum++;
            continue;
        }

        if (opt == 'k') {
            k = atoi(optarg);

            if (k <= 1 || k % 2 == 0) {
                fprintf(stderr, "Número de diagonais de A deve ser ímpar e maior que 1!\n");
                return 1;
            }

            optSum++;
            continue;
        }

        if (opt == 'p') {
            p = atoi(optarg);
            optSum++;
            continue;
        }

        if (opt == 'i') {
            i = atoi(optarg);
            optSum++;
            continue;
        }

        if (opt == 'e') {
            e = atof(optarg);
            continue;
        }

        if (opt == 'o') {
            strcpy(o, optarg);
            optSum++;
            continue;
        }

        if (opt == ':')
            fprintf(stderr, "Todas as opções necessitam de um valor!\n");

        if (opt == '?')
            fprintf(stderr, "Opção %c desconhecida!\n", optopt);

        return 1;
    }

    if (optSum < MIN_OPTS) {
        fprintf(stderr, "Opções insuficientes!\n");
        return 1;
    }

    if (k > n) {
        fprintf(stderr, "O número de diagonais deve ser no máximo a dimensão do SL!\n");
        return 1;
    }

    FILE* file = fopen(o, "w");

    if (!file) {
        fprintf(stderr, "Erro ao abrir arquivo %s\n", o);
        return 1;
    }

    fprintf(file, "# dlpg21 David Lucas Pereira Gomes\n#\n");

    // Run the conjugate gradient method
    if (!p) runConjugateGradient(n, k, i, e, file, conjugateGradient);
    else runConjugateGradient(n, k, i, e, file, conjugateGradientJacobi);

    free(o);
    fclose(file);

    LIKWID_MARKER_CLOSE;

    return 0;
}
