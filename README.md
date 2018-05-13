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

Based on a problem size of 1600 x 1600. All converge in 4580 iterations.

 Ranks  Memory / Node (MB)  Runtime (seconds)
------ ------------------- ------------------
    1               48.85              245.62
    2               48.85              133.19
    4               24.44               74.11
    8               12.23               37.25
   16                6.13               23.28
   32
   64



## Weak Scaling Study
Based on `sinteractive` experimentation, a problem size of 600 take ~11.95 seconds for 1 rank.

Problem size    Ranks   Runtime (seconds)   Memory / Node (MB)  Iterations
------------   ------  ------------------  ------------------- -----------
600 x 600           1               11.95                13.74       1 825
850 x 850           2               17.93                13.79       2 523
1200 x 1200         4               27.59                13.75       3 550
1700 x 1700         8               47.09                13.81       4 609
2400 x 2400        16               86.63                13.77       6 478
3400 x 3400        32              144.23                13.83       9 128
4800 x 4800        64              213.99                13.81      12 566
