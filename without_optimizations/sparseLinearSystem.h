#ifndef SPARSE_LINEAR_SYSTEM_H
#define SPARSE_LINEAR_SYSTEM_H

#include <stdio.h>

/*!
    \brief Sparse linear system's struct

    \details The coefficients matrix is alloc as a unique vector, which has all coefficients.
    The lines of the matrix are stored sequentially, from the first to the last one.
    Each line is stored without the zeros from empty diagonals, thus, some lines
    have different sizes. For example, the first line has only the first element of the
    main diagonal and the higher diagonals. To access the beggining of each line easily,
    the param A receives a vector of pointers to the beggining of each line, and the vector
    with all the elements is allocated individually on A[0].

    \param A pointer to the vector of pointers to the diagonals' vectors of the coefficients matrix
    \param b pointer to the independent terms vector
    \param k number of diagonals
    \param n dimension of the linear system
    \param sizeA size of the diagonals vectors (sequence of diagonals)
*/
typedef struct {
    double** A;
    double* b;

    unsigned int k;
    unsigned int n;
    unsigned int sizeA;
} SparseLinearSystem;

/*!
    \brief Alloc the SparseLinearSystem struct using malloc

    \param k number of diagonals
    \param n dimension of the linear system

    \returns pointer to the SparseLinearSystem
*/
SparseLinearSystem* allocLinearSystem(unsigned int k, unsigned int n);

/*!
    \brief Copy the data from the source SLS to the destination SLS

    \param source SLS to copy data from
    \param destination SLS to write data to
*/
void copyLinearSystem(SparseLinearSystem* source, SparseLinearSystem* destination);

/*!
    \brief Prints SparseLinearSystem's A and b to stream

    \param stream pointer to the stream
    \param sls pointer to the SparseLinearSystem
*/
void fprintLinearSystem(FILE* stream, SparseLinearSystem* sls);

/*!
    \brief Prints vector to the stream

    \param stream pointer to the stream
    \param vector pointer to the vector
    \param n size of the vector
*/
void fprintVector(FILE* stream, double* vector, unsigned int n);

/*!
    \brief Free the SparseLinearSystem

    \param sls pointer to the SparseLinearSystem
*/
void freeLinearSystem(SparseLinearSystem* sls);

/*!
    \brief Given a line i, calculates the first j that has a value different from zero

    \param sls pointer to the SparseLinearSystem
    \param i index of the line
*/
int getBeginIndexLine(SparseLinearSystem* sls, int i);

/*!
    \brief Given a line i, calculates the last j that has a value different from zero

    \param sls pointer to the SparseLinearSystem
    \param i index of the line
*/
int getEndIndexLine(SparseLinearSystem* sls, int i);

/*!
    \brief Given a line i, calculates the size of the vector with values different from zero

    \param sls pointer to the SparseLinearSystem
    \param i index of the line
*/
int getLineSize(SparseLinearSystem* sls, int i);

/*!
    \brief Calculates the transpose of the SparseLinearSystems's coefficients matrix

    \param matrix pointer to the original SparseLinearSystem
    \param transpose pointer to the SparseLinearSystem to save the transpose matrix
*/
void getMatrixTranspose(SparseLinearSystem* matrix, SparseLinearSystem* transpose);

/*!
    \brief Given a line i and a column j of the A matrix, gets the value from the SparseLinearSystem's coefficients vector

    \param sls pointer to the SparseLinearSystem
    \param i index of the line
    \param j index of the column

    \returns The value on i,j from the SparseLinearSystem's coefficients vector
*/
double getMatrixValueByPos(SparseLinearSystem* sls, int i, int j);

/*!
    \brief Multiplies the transpose matrix to the SparseLinearSystem's coefficients matrix

    \param matrix pointer to the SparseLinearSystem
    \param transpose pointer to the transpose SparseLinearSystem
*/
void multiplyTranposeWithA(SparseLinearSystem* matrix, SparseLinearSystem* transpose);

/*!
    \brief Multiplies the transpose matrix to the SparseLinearSystem's b vector

    \param matrix pointer to the SparseLinearSystem
    \param transpose pointer to the transpose SparseLinearSystem
*/
void multiplyTranposeWithB(SparseLinearSystem* matrix, SparseLinearSystem* transpose);

/*!
    \brief Calculates the L2 norm of a vector

    \param x pointer to the vector
    \param n size of the vector

    \returns The L2 norm of x (square root of the sum of squares)
*/
double normL2(double* x, int n);

/*!
    \brief Populates the SparseLinearSystem's A and b with random values

    \param sls pointer to the SparseLinearSystem
*/
void populateLinearSystem(SparseLinearSystem* sls);

/*!
    \brief Given a line i and a column j of the A matrix, sets the value on the SparseLinearSystem's coefficients vector

    \param sls pointer to the SparseLinearSystem
    \param i index of the line
    \param j index of the column
    \param value value to be set on i,j
*/
void setMatrixValueByPos(SparseLinearSystem* sls, int i, int j, double value);

#endif
