# distributed_cg
A conjugate gradient linear solver built on MPI for solving
Poisson's Equation in 2D.

This code was completed for Problem Set 4 in
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

$> ./cg 75 serial_sparse

will run the solver for a 75 x 75 physical domain size in serial.

# Analysis

## Serial Experimentation
Run on 1 thread and node on Midway

    Problem size        mode                  Runtime (seconds)          Memory Usage (MB)
    ---------------     ---------------     -------------------     ----------------------
    75 x 75             `serial_sparse`                    0.02                       0.21
    75 x 75             `serial_dense`                     9.19                     241.61
    10 000 x 10 000     `serial_sparse`                     ---                   3 814.70
    10 000 x 10 000     `serial_dense`                      ---      17 592 186 029 767.20

![75x75 Solution](picture.png)

## Parallel Plot

## Strong Scaling study

## Weak Scaling Study
Based on `sinteractive` experimentation, a problem size of 600 take ~11.95 seconds for 1 rank.
