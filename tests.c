#include "cg_header.h"

#define N 200

void test_dot(int mype, int nprocs){
    double sres, pres, *sy, *py, *sx, *px;
    int size = N/nprocs;

    sx = malloc(sizeof(double) * N);
    sy = malloc(sizeof(double) * N);
    px = malloc(sizeof(double) * size);
    py = malloc(sizeof(double) * size);

    for(int rounds=0; rounds<N; rounds++){
        for(int i=0; i<N; i++){
            sx[i] = (double) rand() / RAND_MAX * 2 - 1;
            sy[i] = (double) rand() / RAND_MAX * 2 - 1;
        }
        for(int i=0; i<size; i++){
            int ii = mype * size + i;
            px[i] = sx[ii];
            py[i] = sy[ii];
        }
        sres = dotp(sx, sy, N);
        pres = parallel_dotp(px, py, N, mype, nprocs);
        assert(abs(sres - pres) < 1e-30);
    }
    free(sx); free(sy); free(py); free(px);
}

void test_axpy(int mype, int nprocs){

    double *sy, *py, *sx, *px, alpha, beta;
    int size = N/nprocs;

    sx = malloc(sizeof(double) * N);
    sy = malloc(sizeof(double) * N);
    px = malloc(sizeof(double) * size);
    py = malloc(sizeof(double) * size);

    for(int rounds=0; rounds<N; rounds++){
        alpha = rand();
        beta = rand();
        for(int i=0; i<N; i++){
            sx[i] = rand();
            sy[i] = rand();
        }
        for(int i=0; i<size; i++){
            int ii = mype * size + i;
            px[i] = sx[ii];
            py[i] = sy[ii];
        }
        axpy(alpha, sx, beta, sy, N);
        parallel_axpy(alpha, px, beta, py, N, mype, nprocs);
        for(int i=0; i<size; i++){
            int ii = mype * size + i;
            assert(abs(sx[ii] - px[i]) < 1e-30);
        }
    }
    free(sx); free(sy); free(py); free(px);
}


void test_matvecOTF(int mype, int nprocs){

    int n, size;
    double *px, *pres, *sx, *sres;

    size = N / nprocs;
    n = sqrt(N);

    sx   = malloc(sizeof(double) * N);
    sres = malloc(sizeof(double) * N);
    pres = malloc(size * sizeof(double));
    px = malloc(sizeof(double) * (n + size + n));
    px += n;

    for(int rounds=0; rounds<N; rounds++){

        for(int i=0; i<N; i++)
            sx[i] = (double) rand();

        for(int i=0; i<size; i++)
            px[i] =  sx[size * mype + i];

        matvec_OTF(sres, sx, N);
        parallel_matvec_OTF(pres, px, N, mype, nprocs);

        for(int pe=0; pe<nprocs; pe++)
            if(mype==pe)
                for(int i=0; i<size; i++)
                    assert(pres[i] == sres[size * mype + i]);
    }

    free(sx); free(sres); free(pres); free(px - n);
}

void test_ghost_cell_exchange(int mype, int nprocs){

    int n, size;
    double *sx, *px;

    n = sqrt(N);
    size = N / nprocs;
    sx = malloc(N * sizeof(double));
    px = calloc(n + size + n, sizeof(double));
    px += n;

    for(int rounds=0; rounds<N; rounds++){

        for(int i=0; i<N; i++)
            sx[i] = (double) rand();

        for(int i=0; i<size; i++)
            px[i] = sx[size * mype + i];

        exchange_ghosts(px, n, size, mype, nprocs);

        for(int i=-n; i<size + n; i++){
            int ii = i + size * mype;
            if(0 <= ii && ii < N)
                assert(sx[ii] == px[i]);
            else
                assert(px[i] == 0); // shouldn't have been touched since calloc
        }
    }
    free(px - n); free(sx);
}


void test_parallel_fill_b(int mype, int nprocs){
    double *sx, *px;
    int size = N/nprocs;
    sx = malloc(N * sizeof(double));
    px = malloc(size * sizeof(double));

    fill_b(sx, N);
    parallel_fill_b(px, N, mype, nprocs);
    for(int i=0; i<size; i++){
        int ii = size * mype;
        assert(abs(sx[ii] - px[i]) < 1e-30);
    }

}


int main(int argc, char const *argv[]) {

    int nprocs, mype;

    MPI_Init(&argc, (char***) &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    if (N % nprocs != 0) {
        if(!mype){
            printf("Expected ranks (%d) to divide %d.\n", nprocs, N);
        }
    } else {
        srand(2018);
        test_dot(mype, nprocs);
        test_axpy(mype, nprocs);
        test_ghost_cell_exchange(mype, nprocs);
        test_matvecOTF(mype, nprocs);
        test_parallel_fill_b(mype, nprocs);
    }

    MPI_Finalize();
    return 0;
}
