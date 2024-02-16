/****************************************************************************************
 * @file mm_triangular.h     															*
 * @addtogroup mm_triangular															*
 * @brief Main macro variables and functions of mm mm_triangular and                    *
 * mm_triangular_triangular                                                             *
 * @author Arturo González Escribano													*
 * @author Rocío Carratalá Sáez															*
 * @author Maria Inmaculada Santamaria Valenzuela										*
 * @author Yuri Torres de la Sierra														*
 ****************************************************************************************/

#ifndef MM_TRIANGULAR
#define MM_TRIANGULAR
/**********************
 * @section Scalapack *
 * ********************/

#ifndef NTNDAT
    #define NTNDAT	10001 //¿Qué significa?
#endif 
#ifndef COMU 
    #define COMU(x) x
#endif 
#ifndef COME
#define COME(x) 
#endif 
#ifndef BIGBLOCK
#define BIGBLOCK 1
#endif 
#ifndef TAMA_BLOQUE
#define TAMA_BLOQUE 16
#endif 
#ifndef MAXINT
#define MAXINT 32768
#endif 
#ifndef CHECK
    #define CHECK 0
#endif
#ifndef CONC
    #define CONC 0
#endif 
#ifndef NORM 
    #define NORM 0
#endif 
#ifndef MININUMS
    #define MININUMS 0
#endif 

/***********************
 * @private @section Por borrar -> Con cuidado de dejar todo con la opción correcta
 * ********************/
    #define BUFFERED_B 1
    #define COPY_B 0
    #define USE_OLD_BALANCED 0
    #define STATUSSES 0

/***********************
 *  @section CONSTANTS *
 ***********************/
/**
 * @private VALUES_LIMIT
 * Constant to limit the values of 
 * initialization with consecutive values with 
 * cyclic scheme 
 * */
#ifndef VALUES_LIMIT
    #define	VALUES_LIMIT 2341
#endif


/**
 * @private BLOCK_LIMIT
 * Constant to limit the num of elements in communication blocks
 * initialization with consecutive values with 
 * cyclic scheme 
 * */
#ifndef BLOCK_LIMIT
    #define	BLOCK_LIMIT __INT_MAX__ ///2 //235929600 //134217728 //235929600
#endif 

/**
 * @private DEBUG
 * If DEBUG is 0, we only print times.
 * If DEBUG is 1, we print everything
 * */
#ifndef DEBUG
    #define DEBUG 0
#endif 

/**
 * @private STATS
 * If STATS is 0, don't calculate times.
 * If STATS is 1, calculate times
 * */
#ifndef STATS
    #define STATS 1
#endif 

/**
 * @private PRINT
 * If PRINT is 1, we print matrices
 * if PRINT is 0, we do not print them
 * */
#ifndef PRINT
#define PRINT 0
#endif 

/**
 * @private ERRROR_CHAR_SIZE
 * @brief Max size for error filename
 * */
#ifndef ERROR_CHAR_SIZE
    #define ERROR_CHAR_SIZE 300 
#endif

/**
 * @section HEADER FILES
 * */
//External 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include <sys/time.h>
#include "omp.h"

//Internal
#include "auxiliar.h"                   //Auxiliar functions for memmory and time handling
#include "mkl_auxiliar.h"               //Adds mkl and new types that makes its use easier
#include "mm_triangular_params.h"       //Types for programm parameters and functions for their lectures
#include "matrix_range.h"               //Definition of matrix type and basic functions associated to it
#include "map_functions.h"              //Matrix mapping function headers
#include "mm_triangular_datatypes.h"    //Datatypes for MPI sending-receive and building functions' headers (shape geberation)
#include "error.h"                      //Error handling datatypes and functions' headers
#include "mm_triangular_matrix_range.h" //MatrixRange functions that needs other libraries to be loaded
#include "mm_triangular_comms.h"        //Communication Functions' headers

/*****************************
 * @section GLOBAL VARIABLES *
 * ***************************/
//MPI
extern int mpi_rank, mpi_procs;     
extern int rows, columns;
extern unsigned int send_block_size_quotient;
//Auxiliar
extern int do_check;
//Parameters
extern Map_function f_mapping;              //Mapping function generic
extern Params_init init_scheme;             //Values generation selected option
extern Params_output output;                //Output generation selected option 
extern Params_comm_scheme comm_scheme;      //Communication Scheme selected
extern Params_comm_func comm_func;          //Communication pattern selected
extern Params_comm_shape comm_shape;        //Communication shape selected
extern Params_comm_mechanism comm_mechanism;//Communication mechanism selected (necessary for triang2type and triang2buff modes)

/**********************
 * Specific functions *
 * ********************/

/**
 * @private @brief printing function
 */
void print_output_results(int rows, int columns, Range local_cols, int mpi_rank,double *mat_C, int mpi_procs, char *func);
void print_output_results_generic(int rows, int columns,Range local_rows,Range local_cols, int mpi_rank,double *mat_C, int mpi_procs,char *func);
/**
 * @private @brief product for triangular*rectangular matrices.
 * @param[in]       mat_A           rectangular matrix
 * @param[in]       mat_B           rectangular matrix
 * @param[in out]   mat_C           rectangular matrix pointer
 * @param[in]       block_before_C  block before C in order to save product in the correct position
 * @param[in]       stage           Current For Product loop stage
 * @param[in]       mpi_rank        Current proccessor identificator
 * */
void mm_triangular_rectangular_product(MatrixRange mat_A, MatrixRange mat_B,MatrixRange *mat_C,size_t block_before_C,unsigned int stage, unsigned int mpi_rank);
/** 
 * @private @brief product for triangular*rectangular product when A is the moved matrix
 * @param[in]       mat_A           rectangular matrix
 * @param[in]       mat_B           rectangular matrix
 * @param[in out]   mat_C           rectangular matrix pointer
 * @param[in]       block_before_C  block before C in order to save product in the correct position
 * @param[in]       stage           Current For Product loop stage
 * @param[in]       mpi_rank        Current proccessor identificator
 */
void mm_triangular_triangular_product(MatrixRange mat_A, MatrixRange mat_B,MatrixRange *mat_C,size_t block_before_C, unsigned int stage, unsigned int mpi_rank);

/**
 * @private mm_rectangular_rectangular_product
 * @brief mat_A * mat_B = mat_C with all them rectangular matrices
 * @param[in]       mat_A           rectangular matrix
 * @param[in]       mat_B           rectangular matrix
 * @param[in out]   mat_C           rectangular matrix pointer
 * @param[in]       block_before_C  block before C in order to save product in the correct position
 * @param[in]       stage           Current For Product loop stage
 * @param[in]       mpi_rank        Current proccessor identificator
 * */
void mm_rectangular_rectangular_product(MatrixRange mat_A, MatrixRange mat_B,MatrixRange *mat_C,size_t block_before_C,unsigned int stage, unsigned int mpi_rank);
/**
 * @private get norm
 * @brief Gets infinity norm of a distributed matrix
 * @param[in]   myMat   Local part
 * @param[in]   mat     Global matrix (just for indexes - data is neither used or needed)
 * @param[in]   comm    MPI communicator
 * @param[out]  norm    Infinity Norm of mat
 * */
double matrix_rangeInfinityNorm(MatrixRange myMat, MatrixRange mat, MPI_Comm comm, int mpi_rank);

#endif
