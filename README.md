# distributed_cg
A conjugate gradient linear solver built on MPI for solving
Poisson's Equation in 2D.

This code was written as a starting point for Problem Set 4 in
MPCS 51087 at the University of Chicago.

## Compilation
Compilation settings can be configured at the top of the included makefile.
Settings include turning on/off MPI.

$> make

## Runtime Settings
- The physical domain size (n) can be set at runtime using the first
command line argument of the program.

- Parallelism can be enabled/disabled by using the second command
line argument of the program, specifying either 'serial' or 'parallel'
execution.

## Running
example:

$> ./cg 75 serial

will run the solver for a 75 x 75 physical domain size in serial.
