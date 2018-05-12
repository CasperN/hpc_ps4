#include "cg_header.h"

// Serial Conjugate Gradient Solver Function for Ax = b
// A must be symmetric and positive definite
// A = 2D operator matrix
// x = 1D solution vector
// b = 1D vector
// N = dimension
void cg_dense(double ** A, double * x, double * b, long N){
    long iter;
    double *r, *p, *z, rs_old;

    r = malloc(N * sizeof(double));
    p = malloc(N * sizeof(double));
    z = malloc(N * sizeof(double));
    matvec(r, A, x, N);                 // r = -A * x + b
    axpy(-1.0, r, 1.0, b, N);
    memcpy(p, r, N * sizeof(double));   // p = r;

    rs_old = dotp(r, r, N);             // rs_old = r' * r;

    for(iter = 0; iter < N; iter++ ){
        double alpha, rs_new;
        matvec(z, A, p, N);                 // z = A * p;
        alpha = rs_old / dotp(p, z, N);     // alpha = rs_old / (p' * z);
        axpy(1.0, x, alpha, p, N);          // x = x + alpha * p;
        axpy(1.0, r, -alpha, z, N);         // r = r - alpha * z;
        rs_new = dotp(r,r,N);               // rs_new = r * r
        if(sqrt(rs_new) < 1.0e-10) break;   // stopping condition
        axpy(rs_new/rs_old, p, 1.0, r, N);  // p = (rs_new / rs_old) * p + r;

        rs_old = rs_new;
    }
    printf("CG converged in %ld iterations.\n", iter);

    free(r);
    free(p);
    free(z);
}

// Serial Conjugate Gradient Solver Function for Ax = b
// where A is the 2D poisson matrix operator that
// is assembled on the fly (not stored)
// x = 1D solution vector
// b = 1D vector
// N = dimension
void cg_sparse_poisson(double * x, double * b, long N)
{

    double *r, *p, *z, rs_old;
    r = malloc(N * sizeof(double));
    p = malloc(N * sizeof(double));
    z = malloc(N * sizeof(double));

    matvec_OTF(r, x, N);
    axpy(-1.0, r, 1.0, b, N);           // r = -A*x + b
    memcpy(p, r, N * sizeof(double));   // p = r;
    rs_old = dotp(r, r, N);             // rs_old = r' * r;

    long iter;
    for(iter = 0; iter < N; iter++ ){
        double alpha, rs_new;
        matvec_OTF(z, p, N);                // z = A * p;
        alpha = rs_old / dotp(p, z, N);     // alpha = rs_old / (p' * z);
        axpy(1.0, x, alpha, p, N);          // x = x + alpha * p;
        axpy(1.0, r, -alpha, z, N);         // r = r - alpha * z;
        rs_new = dotp(r, r, N);             // rs_new = r * r
        if( sqrt(rs_new) < 1.0e-10 ) break; // stopping condition
        axpy(rs_new/rs_old, p, 1.0, r, N);  // p = (rs_new / rs_old) * p + r;

        // printf("% 8.3lf\t% 8.3lf\t% 8.3lf\n",alpha,rs_old,rs_new);
        rs_old = rs_new;
    }

    printf("CG converged in %ld iterations.\n", iter);
    free(r);
    free(p);
    free(z);
}

// General Matrix Vector Product for v = M * w
// v, w = 1D vectors
// M = 2D matrix
// N = dimension
void matvec( double * v, double ** M, double * w, long N )
{
    // Set solution vector to 0
    memset( v, 0, N*sizeof(double));

    for (long i = 0; i < N; i++)
        for (long j = 0; j < N; j++)
            v[i] += (M[i][j] * w[j]);
}

// Specific "On the fly" Matrix Vector Product for v = A * w
// where 'A' is the 2D Poisson operator matrix
// that is applied without explicit storage
// v, w = 1D vectors
// N = dimension
void matvec_OTF(double *v, double *w, long N){
    // Determine physical dimension
    long n = sqrt(N);

    // Set solution vector to 0
    memset(v, 0, N * sizeof(double));

    for(long i = 0; i < N; i++){
        long j = i % n;
        double far_left_diag, far_right_diag, main_diag, left_diag, right_diag;

        // (i/n,j) ~ physical (x,y) coordinates ~ (row, col) in solution matrix
        left_diag      = j != 0   ? -w[i-1] : 0.0;
        right_diag     = j != n-1 ? -w[i+1] : 0.0;
        far_left_diag  = i >= n   ? -w[i-n] : 0.0; // row > 0
        far_right_diag = i <  N-n ? -w[i+n] : 0.0; // row < n
        main_diag = 4 * w[i];

        v[i] = far_left_diag + left_diag + main_diag + right_diag + far_right_diag;
    }
}


// Dot product of c = a * b
// c = result scalar that's returned
// a, b = 1D Vectors
// N = dimension
double dotp( double * a, double * b, long N){
    double c = 0.0;
    for( long i = 0; i < N; i++ )
        c += a[i]*b[i];

    return c;
}


// Scale and add of two vectors (axpy)
// Solves w = alpha * w + beta * v (overwrites what's in w)
// alpha, beta = scalars
// w, v = 1D vectors
// N = dimension
void axpy( double alpha, double * w, double beta, double * v, long N){
    for( long i = 0; i < N; i++ )
        w[i] = alpha * w[i] + beta * v[i];
}


// Fills an Explicit Fully Stored 2D Poisson Operator Matrix 'A'
// A = 2D Matrix
// N = dimension
void fill_A(double ** A, long N){
    long n = sqrt(N);

    for( long i = 0; i < N; i++ ){
        for( long j = 0; j < N; j++ ){
            // Compute x-coordinate of index
            long x = j % n;

            if(i == j)
                A[i][j] = 4.0;          // main diagonal
            else if(i == j + 1){
                A[i][j] = -1.0;
                if(x == n-1)
                    A[i][j] = 0.0;      // left diagonal
            }
            else if(i == j - 1){
                A[i][j] = -1.0;
                if(x == 0)
                    A[i][j] = 0.0;      // right diagonal
            }
            else if(j+n == i)
                A[i][j] = -1.0;         // far left diagonal
            else if(j-n == i)
                A[i][j] = -1.0;         // far right diagonal
            else
                A[i][j] = 0.0;          // Otherwise, A is 0
        }
    }
}

// Sets a circular source term for the right hand side 'b' vector
// Takes as input 2D spatial domain indices
// i, j = indices
// n = physical dimension
double find_b(long i, long j, long n){
    double delta, x, y, radius;

    delta = 1.0 / (double) n;
    x = -.5 + delta + delta * j;
    y = -.5 + delta + delta * i;

    radius = 0.1; // Check if within a circle
    if(x * x + y * y < radius * radius)
        return delta * delta / 1.075271758e-02;
    else
        return 0.0;
}

// Fills a 1-D RHS vector specifying boundary conditions
// by calling the get_b method
// b = 1-D RHS vector
// N = dimension (length of b)
void fill_b(double * b, long N){
    long n = sqrt(N);
    for(long i = 0; i < n; i++ ){
        for(long j = 0; j < n; j++ ){
            long idx = i*n + j;
            b[idx] = find_b(i,j,n);
        }
    }
}
