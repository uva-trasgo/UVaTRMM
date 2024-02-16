#ifndef MATRIX_RANGE
	#define MATRIX_RANGE
	#include "mkl_range.h"

	/**************
 	* STRUCTURES *
 	* ************/
	typedef enum{
		RECTANGULAR, 
		TRIANGULAR_D,
		//TRIANGULAR_U
		//SYMETRIC
		//...
	} Matrix_type;
	/**
 		* @private Structure MatrixRange
	 	* Represents a matrix by ranges.
 		* Used for reducing functions parameters.
 	* */
	typedef struct {
		double *mat;
		Range rows;
		Range columns;
		Matrix_type type;
	} MatrixRange;

	/**
 		* @private Structure MatrixRanges
	 	* Represents a buffer of MatrixRange.
 		* Used for reducing functions parameters.
 	* */
	typedef struct {
		unsigned int total;
		MatrixRange *array;
	} MatrixRanges;

	/***
	 * FUNCTIONS
	 * **/
	MatrixRange matrix_range(Matrix_type type);
	MatrixRanges matrix_ranges(Matrix_type type, unsigned int total);
	void matrix_rangeMat_calloc_abort(MatrixRange *mat);
	void printf_matrix_ranges(double *mat, int row_begin, int row_end, int column_begin, int column_end, const char *name);
	void printf_matrix_ranges_rank_stage(double *mat, int row_begin, int row_end, int column_begin, int column_end, const char *name,int mpi_rank,int stage,int mat_block_before, _Bool check_triangular);
	void matrix_rangePrintf_rank_stage(MatrixRange mat,const char *name,int mpi_rank,int stage,int block_before);
#endif //MATRIX_RANGE
