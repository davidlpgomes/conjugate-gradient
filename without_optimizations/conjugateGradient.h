#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include <stdio.h>

#include "sparseLinearSystem.h"

/*!
    \brief Calculates the alpha used in the conjugate gradient method

    \param sls pointer to the SparseLinearSystem
    \param r pointer to the residue vector
    \param p pointer to the conjugate gradient method's p vector

    \returns The alpha calculated with r and p
*/
double calculateAlpha(SparseLinearSystem* sls, double* r, double* p);

/*!
    \brief Calculates the alpha used in the Jacobi preconditioned conjugate gradient method

    \param sls pointer to the SparseLinearSystem
    \param r pointer to the residue vector
    \param p pointer to the conjugate gradient method's p vector
    \param z pointer to the z vector

    \returns The alpha calculated with r, p and z
*/
double calculateAlphaJacobi(SparseLinearSystem* sls, double* r, double* p, double* z);

/*!
    \brief Calculates the beta used in the conjugate gradient method

    \param sls pointer to the SparseLinearSystem
    \param r pointer to the old residue vector
    \param r1 pointer to the new residue vector

    \returns The beta calculated with r and r1
*/
double calculateBeta(SparseLinearSystem* sls, double* r, double* r1);

/*!
    \brief Calculates the beta used in the Jacobi preconditioned conjugate gradient method

    \param sls pointer to the SparseLinearSystem
    \param r pointer to the old residue vector
    \param r1 pointer to the new residue vector
    \param z pointer to the old z vector
    \param z1 pointer to the new z vector

    \returns The beta calculated with r, r1, z and z1
*/
double calculateBetaJacobi(SparseLinearSystem* sls, double* r, double* r1, double* z, double* z1);

/*!
    \brief Calculates the inverse of the M matrix used in the Jacobi's preconditioning method

    \param sls pointer to the SparseLinearSystem
    \param m pointer to the m vector
*/
void calculateM(SparseLinearSystem* sls, double* m);

/*!
    \brief Calculates the new conjugate gradient method's p vector

    \param sls pointer to the SparseLinear System
    \param p pointer to the conjugate gradient method's p vector
    \param x pointer to the solutions vector
    \param beta conjugate gradient method's beta value
*/
void calculateNextP(SparseLinearSystem* sls, double* p, double* x, double beta);

/*!
    \brief Calculates the new residue of the SparseLinearSystem

    \param sls pointer to the SparseLinearSystem
    \param r pointer to the old residue vector
    \param r1 pointer to the new residue vector
    \param p pointer to the conjugate gradient method's p vector
    \param alpha conjugate gradient method's alpha value
*/
void calculateNextResidue(SparseLinearSystem* sls, double* r, double* r1, double* p, double alpha);

/*!
    \brief Calculates the new solutions vector

    \param sls pointer to the SparseLinearSystem
    \param x pointer to the solutions vector
    \param p pointer to the conjugate gradient method's p vector
    \param alpha conjugate gradient method's alpha value

    \returns The biggest approximate error of old x and new x
*/
double calculateNextX(SparseLinearSystem* sls, double* x, double* p, double alpha);

/*!
    \brief Calculates the residual of the SparseLinearSystem

    \param sls pointer to the SparseLinearSystem
    \param x pointer to the solutions vector
    \param r pointer to the residue vector
*/
void calculateResidue(SparseLinearSystem* sls, double* x, double* r);

/*!
    \brief Calculates the z vector used in the Jacobi's preconditioning method

    \param z pointer to the z vector
    \param m pointer to the M vector
    \param r pointer to the residual vector
    \param n size of the vectors
*/
void calculateZ(double* z, double* m, double* r, int n);

/*!
    \brief Solves the SparseLinearSystem using the Conjugate Gradient method without preconditioning

    \param sls pointer to the SparseLinearSystem
    \param x pointer to the solutions vector
    \param stream pointer to the stream to print informations
    \param maxItr maximum iterations to stop the method
    \param err maximum approximate error to stop the method
*/
void conjugateGradient(SparseLinearSystem* sls, double* x, FILE* stream, int maxItr, double err);

/*!
    \brief Solves the SparseLinearSystem using the Conjugate Gradient method with the Jacobi's preconditioning

    \param sls pointer to the SparseLinearSystem
    \param x pointer to the solutions vector
    \param stream pointer to the stream to print informations
    \param maxItr maximum iterations to stop the method
    \param err maximum approximate error to stop the method
*/
void conjugateGradientJacobi(SparseLinearSystem* sls, double* x, FILE* stream, int maxItr, double err);

#endif
