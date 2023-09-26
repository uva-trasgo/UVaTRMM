#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"
#include <assert.h>
#include "mkl.h"

void blacs_get_(int*, int*, int*);
void blacs_pinfo_(int*, int*);
void blacs_gridinit_(int*, char*, int*, int*);
void blacs_gridinfo_(int*, int*, int*, int*, int*);
void descinit_(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
void pdpotrf_(char*, int*, double*, int*, int*, int*, int*);
void pdtrmm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, int*, int*, double*, int*, int*, int*);
void blacs_gridexit_(int*);
int numroc_(int*, int*, int*, int*, int*);

int main(int argc, char *argv[]) {
    struct timeval start_time, end_time;
    int myrank_mpi, nprocs_mpi;
    MPI_Init( &argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
    int izero=0;
    int ione=1;
    int n = 1000;       // (Global) Matrix size
    int nprow = 2;   // Number of row procs
    int npcol = 2;   // Number of column procs
    int nb = 256;      // (Global) Block size
    char uplo='L';   // Matrix is lower triangular
    char layout='R'; // Block cyclic, Row major processor mapping

    if ( argc > 5 )
    {
    	printf("Usage: ./test matrix_size block_size nprocs_row nprocs_col\n");
	exit(-1);
    }

    if(argc > 1) {
        n = atoi(argv[1]);
    }
    if(argc > 2) {
        nb = atoi(argv[2]);
    }
    if(argc > 3) {
        nprow = atoi(argv[3]);
    }
    if(argc > 4) {
        npcol = atoi(argv[4]);
    }


    assert(nprow * npcol == nprocs_mpi);

    // Initialize BLACS
    int iam, nprocs;
    int zero = 0;
    int ictxt, myrow, mycol;
    blacs_pinfo_(&iam, &nprocs) ; // BLACS rank and world size
    blacs_get_(&zero, &zero, &ictxt ); // -> Create context
    blacs_gridinit_(&ictxt, &layout, &nprow, &npcol ); // Context -> Initialize the grid
    blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol ); // Context -> Context grid info (# procs row/col, current procs row/col)

    // Compute the size of the local matrices
    int mpA    = numroc_( &n, &nb, &myrow, &izero, &nprow ); // My proc -> row of local A
    int nqA    = numroc_( &n, &nb, &mycol, &izero, &npcol ); // My proc -> col of local A

    printf("Hi. Proc %d/%d for MPI, proc %d/%d for BLACS in position (%d,%d)/(%d,%d) with local matrix %dx%d, global matrix %d, block size %d\n",myrank_mpi,nprocs_mpi,iam,nprocs,myrow,mycol,nprow,npcol,mpA,nqA,n,nb);
    fflush(stdout);


    // Allocate and fill the matrices A and B
    // A[I,J] = (I == J ? 5*n : I+J)
    double *A, *B, contador = 0.0;
    A = (double *)calloc(mpA*nqA,sizeof(double)) ;
    B = (double *)calloc(mpA*nqA,sizeof(double)) ;
    if (A==NULL || B==NULL){ printf("Error of memory allocation A or B on proc %dx%d\n",myrow,mycol); exit(0); }
    int k = 0;
    for (int j = 0; j < nqA; j++) { // local col
        int l_j = j / nb; // which block
        int x_j = j % nb; // where within that block
        int J   = (l_j * npcol + mycol) * nb + x_j; // global col
        for (int i = 0; i < mpA; i++) { // local row
            int l_i = i / nb; // which block
            int x_i = i % nb; // where within that block
            int I   = (l_i * nprow + myrow) * nb + x_i; // global row
            assert(I < n);
            assert(J < n);
            if(I >= J) {
                A[k] = 9*iam + contador; contador++;
		B[k] = 1.0; //n*n;
            } else {
                A[k] = 0.0; //9*iam + contador; contador++;
		B[k] = 2.0;
            }
            //printf("A: %d %d -> %d %d -> %f\n", i, j, I, J, A[k]);
            k++;
        }
    }
    k = 0;
    for (int j = 0; j < nqA; j++) { // local col
        int l_j = j / nb; // which block
        int x_j = j % nb; // where within that block
        int J   = (l_j * npcol + mycol) * nb + x_j; // global col
        for (int i = 0; i < mpA; i++) { // local row
            int l_i = i / nb; // which block
            int x_i = i % nb; // where within that block
            int I   = (l_i * nprow + myrow) * nb + x_i; // global row
            assert(I < n);
            assert(J < n);
            //printf("B: %d %d -> %d %d -> %f\n", i, j, I, J, B[k]);
            k++;
        }
    }

    // Create descriptor
    int descA[20], descB[20];
    int info;
    int lddA = mpA > 1 ? mpA : 1;
    descinit_( descA, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &lddA, &info);
    descinit_( descB, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &lddA, &info);
    if(info != 0) {
        printf("Error in descinit, info = %d\n", info);
    }

    // Run dpotrf and time
    double MPIt1 = MPI_Wtime();
    gettimeofday(&start_time, NULL);
    printf("[%dx%d] Starting trmm\n", myrow, mycol);
    char SIDE = 'L', UPLO = 'L', TRANS='N', DIAG = 'N';
    double ALPHA = 1.0; /* Check MMTRIANG code */
    int IA = 1, JA = 1;
    int IB = 1, JB = 1;
    pdtrmm_( &SIDE, &UPLO, &TRANS, &DIAG,
               &n, &n,
               &ALPHA,
               A, &IA, &JA, descA,
               B, &IB, &JB, descB);
    MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&end_time, NULL);
    double MPIt2 = MPI_Wtime();
    printf("[%dx%d] Done, time %e s.\n", myrow, mycol, MPIt2 - MPIt1);
    double diff_sec    = (double) end_time.tv_sec - start_time.tv_sec;
    double diff_usec   = (double) end_time.tv_usec - start_time.tv_usec;
    double total_time  = (double) diff_sec + diff_usec/1000000.0;
    printf("Time: %lf\n", total_time);

    fflush(stdout);

/*
    printf("RESULT\n");
    k = 0;
    for (int j = 0; j < nqA && iam == 0; j++) { // local col
        int l_j = j / nb; // which block
        int x_j = j % nb; // where within that block
        int J   = (l_j * npcol + mycol) * nb + x_j; // global col
        for (int i = 0; i < mpA; i++) { // local row
            int l_i = i / nb; // which block
            int x_i = i % nb; // where within that block
            int I   = (l_i * nprow + myrow) * nb + x_i; // global row
            assert(I < n);
            assert(J < n);
            printf("B: %d %d -> %d %d -> %f\n", i, j, I, J, B[k]);
            k++;
	}
    }
*/

    free(A);
    free(B);

    // Exit and finalize
    blacs_gridexit_(&ictxt);
    MPI_Finalize();
    return 0;
}
