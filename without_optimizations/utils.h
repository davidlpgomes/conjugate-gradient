#ifndef UTILS_H
#define UTILS_H

/*!
    \brief Generates the coefficient of a k-diagonal linear system

    \param i line of the matrix
    \param j column of the matrix
    \param k number of diagonals of the matrix

    \returns Random coefficient of the linear system's A
*/
double generateRandomA(unsigned int i, unsigned int j, unsigned int k);

/*!
    \brief Generates the independent term of a k-diagonal lineaer system

    \param k number of diagonals of the matrix

    \returns Random independent term of the linear system's b
*/
double generateRandomB(unsigned int k);

/*!
    \brief Calculates the size of the matrix's coefficients vector (sequence of diagonals)

    \param k total number of diagonals
    \param n dimension of the matrix

    \return Size of the vector, returns 0 if k is invalid to the given n
*/
int getSparseMatrixSize(unsigned int k, unsigned int n);

/*!
    \brief Get the maximum from the given values

    \param a first value
    \param b second value

    \returns The greatest value
*/
int max(int a, int b);

/*!
    \brief Get the minimum from the given values

    \param a first value
    \param b second value

    \returns The smallest value
*/
int min(int a, int b);

/*!
    \brief Calculates the sum of the multiplication of two vectors

    \param v1 pointer to the first vector
    \param v2 pointer to the second vector
    \param n size of the vectors

    \returns The sum of the multiplication of the elements of the vectors of equal positions
*/
double multiplyVectors(double* v, double* t, int n);

/*!
    \brief Calculates the sum of squares of the vector's elements

    \param x pointer to the vector
    \param n size of the vector

    \returns The sum of squares of the vector's elements
*/
double sumOfSquares(double* x, int n);

/*!
    \brief Test if target is alloc

    \param target pointer to target
    \param name name of the target

    \return 0 if success, -1 if failure
*/
int testAlloc(void* target, char* name);

/*!
    \brief Get the current timestamp

    \return Current timestamp
*/
double timestamp();

#endif
