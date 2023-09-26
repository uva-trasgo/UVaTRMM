
/***
 * Functions for range involving mkl types
 * **/
#ifndef MKL_RANGE
#define MKL_RANGE
#include "mkl_auxiliar.h"

/**************
 * STRUCTURES *
 * ************/
/**
 * @private Structure Range
 * Represents ranges of rows assigned to a process:
 * * Begin is the index of first row/columns
 * * Size is total number of rows/columns 
 * */
typedef struct {
	unsigned long int begin;
	unsigned long int size;
} Range;
/**
 * @private Structure Ranges
 * Represents a buffer of ranges 
 * */
typedef struct {
	unsigned int total;
	Range *array;
} Ranges;

int rangePrintf(const char *name, Range range);
int rangePrintf_rank_stage(const char *name, Range range, int mpi_rank, int stage);
Ranges rangesCalloc_abort(unsigned int total, char *name, MPI_Comm comm, int mpi_rank);
int rangesPrint_sizes(Ranges ranges, char *name);
#endif //MKL_RANGE
