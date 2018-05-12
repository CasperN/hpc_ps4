#include "cg_header.h"

// Parallel Distributed MPI Version of
// Conjugate Gradient Solver Function for Ax = b
// where A is the 2D poisson matrix operator that
// is assembled on the fly (not stored)
// x = 1D solution vector
// b = 1D vector
// N = dimension
void parallel_cg_sparse_poisson(double *x, double *b, long N, int mype, int nprocs){
    double *r, *p, *z, rs_old;
    long size = N / nprocs;
    long n = sqrt(N);

    z = malloc(size * sizeof(double));
    r = malloc(size * sizeof(double));
    p = malloc((n + size + n) * sizeof(double));
    p += n;
    // p has ghost cells, which are used in parallel_matvec_OTF
    // Data:     [left ghosts  |  actual vector    |   right ghosts)
    // Indices:  -n            0               size^         size+n^

    memcpy(p, x, size * sizeof(double));            // p = x (need ghost cells for A.x)
    parallel_matvec_OTF(r, p, N, mype, nprocs);     // r =  A.x
    parallel_axpy(-1, r, 1, b, N, mype, nprocs);    // r = -A.x + b
    memcpy(p, r, size * sizeof(double));            // p = r
    rs_old = parallel_dotp(r, r, N, mype, nprocs);  // rs_old = r.r

    long iter;
    for(iter=0; iter<N; iter++){
        double alpha, rs_new;
        parallel_matvec_OTF(z, p, N, mype, nprocs);             // z = A.p
        alpha = rs_old / parallel_dotp(p, z, N, mype, nprocs);  // alpha = rs_old /(p.z)
        parallel_axpy(1, x,  alpha, p, N, mype, nprocs);        // x = x + alpha * p
        parallel_axpy(1, r, -alpha, z, N, mype, nprocs);        // r = r - alpha * z
        rs_new = parallel_dotp(r, r, N, mype, nprocs);          // rs_new = r.r
        if (sqrt(rs_new) < 1.0e-10) break;                      // stopping condition
        parallel_axpy(rs_new/rs_old, p, 1, r, N, mype, nprocs); // p = (rs_new / rs_old) * p + r;
    }

    if(!mype) printf("CG converged in %ld iterations.\n", iter);
    free(p - n);
    free(r);
    free(z);
}

void exchange_ghosts(double *w, long n, long size, int mype, int nprocs){
    // `w` is a double vector of length size+2n.
    // It has two ghost cell buffers of len `n` surrounding data of len `size`.
    // Data:     [left ghosts  |  actual vector    |   right ghosts)
    // Indices:  -n            0               size^        size+ n^

    if(mype % 2){
        if(mype > 0){
            MPI_Sendrecv(
                w    , n, MPI_DOUBLE, mype - 1, 'l',        // send n L elems
                w - n, n, MPI_DOUBLE, mype - 1, 'r',        // recv into L ghost cells
                MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        }
        if (mype < nprocs - 1){                      // evens trade w/ R neighbor second
            MPI_Sendrecv(
                w + size - n , n, MPI_DOUBLE, mype + 1, 'r', // send n R elems
                w + size,      n, MPI_DOUBLE, mype + 1, 'l', // recv into R ghost cells
                MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        }
    } else {
        if (mype < nprocs - 1){                      // odds trade w/ R neighbor first
            MPI_Sendrecv(
                w + size - n, n, MPI_DOUBLE, mype + 1, 'r', // send n R elems
                w + size    , n, MPI_DOUBLE, mype + 1, 'l', // recv into R ghost cells
                MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        }
        if(mype > 0){                                // odds trade w/ L neighbor second
            MPI_Sendrecv(
                w,     n, MPI_DOUBLE, mype - 1, 'l',        // send n L elems
                w - n, n, MPI_DOUBLE, mype - 1, 'r',        // recv into L ghost cells
                MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        }
    }
}


// Parallel Distributed MPI Version of
// Specific "On the fly" Matrix Vector Product for v = A * w
// where 'A' is the 2D Poisson operator matrix
// that is applied without explicit storage
// v, w = 1D vectors
// N = dimension
void parallel_matvec_OTF(double *v, double *w, long N, int mype, int nprocs){

    long my_start, n, size;

    n = sqrt(N);
    size = N / nprocs;
    my_start = size * mype;
    memset(v, 0, size * sizeof(double));    // Set solution vector to 0

    exchange_ghosts(w, n, size, mype, nprocs);
    for(long i=0; i<size; i++){
        double far_left_diag, far_right_diag, main_diag, left_diag, right_diag;
        long I, j;

        I = my_start + i;
        j = I % n;
        // (I/n,j) ~ physical (x,y) coordinates ~ (row, col) in solution matrix
        left_diag      = j != 0   ? -w[i-1] : 0.0;
        right_diag     = j != n-1 ? -w[i+1] : 0.0;
        far_left_diag  = I >= n   ? -w[i-n] : 0.0; // row > 0
        far_right_diag = I <  N-n ? -w[i+n] : 0.0; // row < n
        main_diag = 4 * w[i];

        v[i] = far_left_diag + left_diag + main_diag + right_diag + far_right_diag;
    }
}

// Parallel Distributed MPI Version of
// Dot product of c = a * b
// c = result scalar that's returned
// a, b = 1D Vectors
// N = dimension
double parallel_dotp(double *a, double *b, long N, int mype, int nprocs){
    double global_sum, sum = 0;
    long size = N / nprocs;

    for(long i=0; i < size; i++)
        sum += a[i] * b[i];

    // printf("rank:%d sum:%lf\n", mype, sum);
    MPI_Reduce(&sum,  &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&global_sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // if(!mype) printf("global: %lf\n", global_sum);

    return global_sum;
}

// Parallel Distributed MPI Version of
// Scale and add of two vectors (axpy)
// Solves w = alpha * w + beta * v (overwrites what's in w)
// alpha, beta = scalars
// w, v = 1D vectors
// N = dimension
void parallel_axpy(double a, double *x, double b, double *y, long N, int mype, int nprocs){
    long size = N / nprocs;

    for(long i=0; i<size; i++){
        x[i] = a * x[i] + b * y[i];
    }

}


// Parallel Distributed MPI Version of
// Fills a 1-D RHS vector specifying boundary conditions
// by calling the get_b method
// b = 1-D RHS vector
// N = dimension (length of b)
void parallel_fill_b(double *b, long N, int mype, int nprocs){

    long size = N/nprocs, n = sqrt(N);
    
    for(long i=0; i<size; i++){
        int ii = mype * size + i;
        b[i] = find_b(ii/n, ii%n, n);
    }
}


void parallel_save_vector(double *vec, long N, char* fname, int mype, int nprocs){
    MPI_File f;
    long start, size;
    size = N / nprocs;
    start = size * mype;

    MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f);
    MPI_File_write_at_all(f, start * sizeof(double), vec, size, MPI_DOUBLE, MPI_STATUSES_IGNORE);
    MPI_File_close(&f);
}
