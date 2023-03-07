# Conjugate Gradient

An implementation of the conjugate gradient method for solving sparse linear system. It solves it with or without the Jacobi preconditioner. The `without_optimizations` directory has the code without the optimizations implemented on the main code, such as loop unroll, which allows the compiler to use SIMD instruction and the XMM/YMM registers.

## Build
Run `make`.

## Usage
``` 
usage: ./cgSolver -n N -k K -p P -i I [-e E] -o O

options:
  -n N         N = size of the linear system (N x N). It must be bigger than 11
  -k K         K = number of diagonals, it must be odd, at least 1 and at most N
  -p P         P = 0 for conjugate gradient without preconditioning, any other number for Jacobi preconditioner
  -i I         I = max iterations
  -e E         E = minimum error to finish iterations (optional)
  -o O         O = output file
```

## Testing
The `run` shell script runs the algorithm on multiple sizes of linear system, while calculates the time spent, miss ratio on the L2 cache, bandwidth of the L3 cache, double precision MFLOP/s and AVX MFLOP/s, using [LIKWID](https://github.com/RRZE-HPC/likwid). Then, saves the results on CSVs, which can be used on the `plot.py` script to compare the algorithm performance with and without optimizations.
