/***
 * Checkeo a lo bruto de soluciones
 * **/

/**
 * Ultima version original sin adaptaciones
 * mm_triangular_lapack.c
 * 	Multiplication of distributed square matrices, one of them lower triangular
 *
 * v4.0
 */
/*
 * mm_triangular.c
 * 	Multiplication of distributed square matrices, one of them lower triangular
 *
 * v2.3
 */
#define	VALUES_LIMIT	2341
#define DEBUG 0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "mkl.h"
int check;

/*#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>
#include "mkl.h"*/

/* Enum types for selector parameters */
typedef enum {
	INIT_VALUES, INIT_A_ID, INIT_B_ID
	} Params_init;
/* Main program */
int main( int argc, char *argv[] ) {
	/* 2. Read parameters */
	if ( argc < 5 ) {
		fprintf( stderr, "Usage: %s <file_prefix> <format_size> <format_dec> <num_procs>\n\n", argv[0] );
		exit( EXIT_FAILURE );
	}
	if ( strlen( argv[1] ) > 60 ) {
		fprintf( stderr, "Error: file_prefix should be 60 chars at most\n" );
		fprintf( stderr, "Usage: %s <file_prefix> <num_procs>\n\n", argv[0] );
		exit( EXIT_FAILURE );
	}
	int format1 = atoi( argv[2] );
	int format2 = atoi( argv[3] );
	int procs = atoi( argv[4] );

	// Select initialization scheme
	//int do_check = 0;
	//Params_init init_scheme;
	//if ( ! strcmp( argv[5], "values" ) )
	//	init_scheme = INIT_VALUES;
	//else if ( ! strcmp( argv[5], "values_check" ) ) {
    //            init_scheme = INIT_VALUES;
    //            //do_check = 1;
    //    }
	//else if ( ! strcmp( argv[5], "a_id" ) )
	//	init_scheme = INIT_A_ID;
	//else if ( ! strcmp( argv[5], "a_id_check" ) ) {
	//	init_scheme = INIT_A_ID;
	//	//do_check = 1;
	//}
	//else if ( ! strcmp( argv[5], "b_id" ) ){
	//	init_scheme = INIT_B_ID;
    //}
	//
	
	/* 2. Open first file to get matrix size */
	char file_name[64];
	char file_path[36]="./tests/outputs/";//./tests/outputs/\0";
	//sprintf( file_name, "%s.0.dat", argv[1] );
	sprintf( file_name, "%s%s.0.dat", file_path, argv[1] );
	FILE *file = fopen( file_name, "rb" );
	if ( file == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", file_name );
		fprintf( stderr, "Usage: %s <file_prefix> <num_procs>\n\n", argv[0] );
		//exit( EXIT_FAILURE );
		check=-2;
	}

	int sizes[2];
	int sizes_local[2];
	int ranges[4];
	fread(sizes, sizeof(int), 2, file);
	fclose( file );


    /* 3. Allocate and initialize matrix */
	#if DEBUG
	printf( "Reading matrix of %d x %d\n", sizes[0], sizes[1] );
	#endif 
    //El leido
	double *mat = (double *)calloc( (size_t)sizes[0] * sizes[1], sizeof(double) );

    int rows = sizes[0];
    int columns = sizes[1];
	/* 4. Initialize local parts of data structures */
	double *mat_A = (double *)calloc( rows*columns, sizeof( double ));
	double *mat_B = (double *)calloc( rows*columns,sizeof( double )); 
	double *mat_C = (double *)calloc( rows*columns,sizeof( double ));
	if ( mat_A == NULL || mat_B == NULL || mat_C == NULL ) {
		fprintf( stderr, "Error: Allocating memory for local parts of the matrices\n" );
        exit(EXIT_FAILURE);
	}
	#if DEBUG
    printf("--> Seting up A & B\n");
	#endif 
    fflush(stdout);
	long int i,j;
	for ( i=0; i<rows; i++ ) {
		for ( j=0; j<columns; j++ ) {
			mat_C[ i * columns + j ] = 0;
			if ( j >  i ) {
				mat_A[ i * columns + j ] = 0;
				mat_B[ i * columns + j ] = 0;
			} else {
				mat_A[ i * columns + j ] = ( ( (i * columns ) + j ) % VALUES_LIMIT + 1);
				mat_B[ i * columns + j ] = ( ( (i * columns ) + j ) % VALUES_LIMIT + 1);
			}
		}
	}
	#if DEBUB > 1
	for (i = 0; i < rows; i++){
		for (j = 0; j < columns; j++){
			printf("A[%ld][%ld] = %.0lf\n", i, j, mat_A[i*columns+j]);
		}
		printf("-- A[%ld] --\n", i);
	}


	for (i = 0; i < rows; i++){
		for (j = 0; j < columns; j++){
			printf("B[%ld][%ld] = %.0lf\n", i, j, mat_B[i*columns+j]);
		}
		printf("-- B[%ld] --\n", i);
	}
	#endif 



	/* 4. Read data */
	int proc;
	#if DEBUG
    printf("--> Reading Output\n");
    fflush(stdout);
	#endif 
	for( proc = 0; proc < procs; proc++ ) {
		// Open file
		//sprintf( file_name, "%s.%d.dat", argv[1], proc );
		sprintf( file_name, "%s%s.%d.dat", file_path, argv[1],proc );
		
		file = fopen( file_name, "rb" );
		if ( file == NULL ) {
			fprintf( stderr, "Error: File %s not found\n", file_name );
			fprintf( stderr, "Usage: %s <file_prefix> <num_procs>\n\n", argv[0] );
			exit( EXIT_FAILURE );
		}

		// Check sizes
		fread(sizes_local, sizeof(int), 2, file);
		if ( sizes[0] != sizes_local[0] || sizes[1] != sizes_local[1] ) {
			fprintf( stderr, "Error: Different matrix sizes in file %s, original (%d,%d), current (%d,%d)\n", file_name, sizes[0], sizes[1], sizes_local[0], sizes_local[1] );
			exit( EXIT_FAILURE );
		}

		// Read ranges
		/*
			local_rows.begin, 
			ranges[0]+local_rows.begin, 
			local_cols.begin, 
			local_cols.size };

		*/
		fread(ranges, sizeof(int), 4, file);
		if ( ranges[0] < 0 || ranges[2] < 0 || ranges[1] > sizes_local[0] || ranges[3] > sizes_local[1] ) {
			fprintf( stderr, "Error: Ranges outside the matrix sizes in file %s, matrix sizes (%d,%d), ranges (%d:%d,%d:%d)\n", file_name, sizes_local[0], sizes_local[1], ranges[0], ranges[1], ranges[2], ranges[3] );
			exit( EXIT_FAILURE );
		}

		// Read data
		int blocksize = ranges[3];
		for ( i=ranges[0]; i<ranges[1]; i++ ) {
			int result = fread( &mat[ i*sizes[1] + ranges[2] ], sizeof(double), blocksize, file);
			if ( result != blocksize ) {
				fprintf( stderr, "Error: Not enough data to read in file %s\n", file_name );
				exit( EXIT_FAILURE );
			}
		}

		fclose( file );
	}

	
	/* 5. Check result */
	char format[10];
    check = 1;
	sprintf( format, "%%%d.%dlf ", format1, format2 );
	char formatresult[23];
	char formatc[13];
	sprintf( formatresult, "Read Result: %%%d.%dlf ", format1, format2 );
	sprintf( formatc, "C: %%%d.%dlf ", format1, format2 );

//    int k;
	#if DEBUG
    printf("--> Making product\n");
    fflush(stdout);
	#endif 
    for (int a = 0; a < columns; a++) {
        // Dentro recorremos las filas de la primera (A)
        for (int i = 0; i < rows; i++) {
            mat_C[i*columns+a] = 0;
            // Y cada columna de la primera (A)
            for (int j = 0; j < columns; j++) {
                // Multiplicamos y sumamos resultado
                mat_C[i*columns+a] += mat_A[i*columns+j] * mat_B[j*columns+a];
            }
        }
    }

	#if DEBUG > 1
    printf("Matrix A: (%dx%d)\n", rows, columns);
    for (int row=0; row < rows; row++) {
		// Only lower triangular part
        fflush(stdout);
        for ( j=0; j < columns; j++ ) {
            printf( format, mat_A[ row*sizes[1]+j ] );
            fflush(stdout);
        }
        printf("\n");
	}
    printf("Matrix B: (%dx%d)\n", rows, columns);
    for (int row=0; row < rows; row++) {
		// Only lower triangular part
        fflush(stdout);
        for ( j=0; j < columns; j++ ) {
            printf( format, mat_B[ row*sizes[1]+j ] );
            fflush(stdout);
        }
        printf("\n");
	}
    printf("Matrix C: (%dx%d)\n", rows, columns);
    for (int row=0; row < rows; row++) {
		// Only lower triangular part
        fflush(stdout);
        for ( j=0; j < columns; j++ ) {
            printf( format, mat_C[ row*sizes[1]+j ] );
        }
        printf("\n");
	}
	printf("Matrix C_read: (%dx%d)\n", rows, columns);
    for (int row=0; row < rows; row++) {
		// Only lower triangular part
        fflush(stdout);
        for ( j=0; j < columns; j++ ) {
            printf( format, mat[ row*sizes[1]+j ] );
        }
        printf("\n");
	}
	#endif 
		
    #if DEBUG
    printf("--> Checking result\n");
    fflush(stdout);
	#endif 
    //Se pueden juntar y usar solo la matriz mat
	for ( i=0; i<sizes[0]; i++ ) {
		for ( j=0; j<sizes[1]; j++ ) {
            if (mat_C[i*columns+j] != mat[i*columns+j]){
			    printf("FAILED C[%ld][%ld]\n", i, j);
                printf(formatresult, mat[ i*sizes[1]+j ] );
				
                fflush(stdout);
                printf(formatc, mat_C[ i*sizes[1]+j ] );
                printf("\nFAILED C[%ld][%ld]\n", i, j);
				check = 0;
            }
		}
		#if DEBUG > 1
		printf("Row %ld check %d \n", i, check);
		#endif 
	}

	/* 6. End */
    printf("Check > %d\n", check);
    fflush(stdout);

    return 0;
}