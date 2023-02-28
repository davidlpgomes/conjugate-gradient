# Conjugate Gradient

An implementation of the conjugate gradient method for solving of sparse linear systems. It solves the systems with or without the Jacobi's preconditioning.

## Build
Run `make`.

## Usage
``` 
usage: ./cgSolver -n N -k K -p P -i I [-e E]

options:
  -n N                  N = size of the linear system (N x N). It must be bigger than 11
  -k K                  K = number of diagonals, it must be odd, at least 1 and at most N
  -p P                  P = 0 for conjugate gradient without preconditioning, any other number for Jacobi preconditioner
  -i I                  I = max iterations
  -e E                  E = minimum error to finish iterations (optional)
```
