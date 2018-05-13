#include "cg_header.h"

int main(int argc, char* argv[])
{
    int nprocs = 1;
    int mype = 0;

    #ifdef MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    #endif

    // n is the physical domain size;
    int n = 1;

    // Read n from first command line argument
    if( argc != 3 )
        cli_error();
    else
        n = atoi(argv[1]);

    // Make sure n is reasonable
    assert(n > 0 && isfinite(n));

    if(mype==0)
        printf("Solving Poisson Equation on %d x %d domain...\n", n, n);

    // Use 2nd CLI argument to determine serial or parallel solve
    if (!strcmp(argv[2], "serial_dense"))
        run_dense(n);

    else if (!strcmp(argv[2], "serial_sparse"))
        run_sparse(n);

    else if (!strcmp(argv[2], "parallel"))
        #ifdef MPI
        run_parallel_sparse(n, mype, nprocs);
        #else
        printf("Fatal: MPI version not compiled\n");
        #endif

    else if(!mype)
        cli_error();

    #ifdef MPI
    MPI_Finalize();
    #endif

    return 0;
}

/////////////////////////////////////////////////////////////////
// Run "naive" CG solve, where A is fully stored in dense format
/////////////////////////////////////////////////////////////////
void run_dense(long n)
{
    printf("Running Dense CG solve...\n");

    // Dimension of operator matrix and vectors is n^2
    int N = n*n;

    // Allocate full A matrix and vectors
    double ** A = matrix( N );
    double *  x = (double*) calloc(N, sizeof(double));
    double *  b = (double*) calloc(N, sizeof(double));
    printf("Dense Memory  = %.2lf MB\n", (N*N+5*N)*sizeof(double)/1024.0/1024.0);

    // Compute elements of 'A' matrix (Poisson Operator)
    fill_A(A, N);

    // Compute elements of boundary condition vector 'b'
    fill_b(b, N);

    // Run Dense CG Solve
    double start = get_time();
    cg_dense(A, x, b, N);
    double stop = get_time();
    printf("Dense Runtime = %.2lf seconds\n", stop-start);

    // Save Solution Vector to File
    save_vector(x,N, "out/dense");

    // Free A matrix
    matrix_free(A);

    // Free vectors
    free(x);
    free(b);
}

/////////////////////////////////////////////////////////////////
// Run optimized CG solve, where A is assembled on the fly (OTF)
/////////////////////////////////////////////////////////////////
void run_sparse(long n)
{
    printf("Running Sparse CG solve...\n");

    double *x, *b, start, stop;
    int N = n*n; // Dimension of operator matrix and vectors is n^2

    // reset vectors
    x = calloc(N, sizeof(double));
    b = calloc(N, sizeof(double));
    printf("Sparse Memory = %.2lf MB\n", (5*N)*sizeof(double)/1024.0/1024.0);

    // Compute elements of boundary condition vector 'b'
    fill_b(b, N);

    // Run Sparse CG Solve
    start = get_time();
    cg_sparse_poisson(x, b, N);
    stop = get_time();
    printf("Sparse Runtime = %.2lf seconds\n", stop-start);

    // Save Solution Vector to File
    save_vector(x,N, "out/sparse");

    // Free vectors
    free(x);
    free(b);
}

/////////////////////////////////////////////////////////////////
// Run Domain Decomposed Parallel CG solve w/MPI
// where A is assembled on the fly (OTF)
/////////////////////////////////////////////////////////////////
#ifdef MPI
void run_parallel_sparse(long n, int mype, int nprocs)
{
    double *x, *b, start, stop;
    long N = n*n;
    long size = N / nprocs;

    if(!mype) printf("Running Parallel Sparse CG solve with %d nodes...\n", nprocs);

    if(N % nprocs){
        if(!mype) printf("Error: problem assumed divisible by number of ranks\n");
        exit(1);
    }

    x = calloc(size, sizeof(double));
    b = calloc(size, sizeof(double));

    if(!mype){
        double mem = (5 * size + 2 * n) * sizeof(double) / 1024.0 / 1024.0;
        printf("Parallel Sparse Memory = %.2lf MB per node, %.2lf MB total \n", mem, nprocs * mem);
    }

    parallel_fill_b(b, N, mype, nprocs);

    start = get_time();
    parallel_cg_sparse_poisson(x, b, N, mype, nprocs);
    stop = get_time();

    if(!mype) printf("Parallel Sparse Runtime = %.2lf seconds\n", stop-start);

    parallel_save_vector(x, N, "out/parallel", mype, nprocs);

    free(x);
    free(b);
}
#endif
