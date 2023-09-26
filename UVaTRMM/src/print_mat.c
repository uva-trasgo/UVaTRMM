#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define DEBUG 0
int main( int argc, char *argv[] ) {
	/* 1. Read arguments */
	if ( argc < 4 ) {
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

	/* 2. Open first file to get matrix size */
	char file_name[64];
	sprintf( file_name, "%s.0.dat", argv[1] );
	FILE *file = fopen( file_name, "rb" );
	if ( file == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", file_name );
		fprintf( stderr, "Usage: %s <file_prefix> <num_procs>\n\n", argv[0] );
		exit( EXIT_FAILURE );
	}

	int sizes[2];
	int sizes_local[2];
	int ranges[4];
	fread(sizes, sizeof(int), 2, file);
	fclose( file );

	/* 3. Allocate and initialize matrix */
	#if DEBUG
	printf( "Reading matrix of %d x %d\n", sizes[0], sizes[1] );
	fflush(stdout);
	#endif 
	double *mat = (double *)calloc( (size_t)sizes[0] * sizes[1], sizeof(double) );

	/* 4. Read data */
	int i,j;
	int proc;
	for( proc = 0; proc < procs; proc++ ) {
		// Open file
		sprintf( file_name, "%s.%d.dat", argv[1], proc );
		#if DEBUG
		printf("File: %s\n", file_name);
		fflush(stdout);
		#endif 
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
	
	/* 5. Print data in text format */
	char format[10];
	sprintf( format, "%%%d.%dlf ", format1, format2 );
	for ( i=0; i<sizes[0]; i++ ) {
		for ( j=0; j<sizes[1]; j++ ) {
			printf( format, mat[ i*sizes[1]+j ] );
		}
		printf("\n");
	}

	/* 6. End */
	return 0;	
}

