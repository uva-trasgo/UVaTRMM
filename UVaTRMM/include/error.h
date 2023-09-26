/****************************************************************************************
 * @file auxiliar.h     																*
 * @addtogroup mm_triangular															*
 * @brief Auxiliar functions and macro                                                  *
 * @author Arturo González Escribano													*
 * @author Rocío Carratalá Sáez															*
 * @author Maria Inmaculada Santamaria Valenzuela										*
 * @author Yuri Torres de la Sierra														*
 ****************************************************************************************/
#ifndef ERROR_H
    #define ERROR_H
    #include "mkl_range.h"
    /***
     * @section TYPE DEFINITION
     * */
    /***
    * @private MM_TRIANGULAR_ERROR
    * @brief Codigos de error
    * */
    typedef enum {
        INVALID_PROCS_NUM,
        INVALID_ARGS,
        INVALID_COMMUNICATION_SHAPE,
        ERROR_OPENING_FILE,
        MEMORY_ERROR,
        INVALID_BLOCKS_NUM,
        INVALID_ROWS_BEFORE,
        EMPTY_BLOCK,
        HORIZONTAL_MATRIX,
        INVALID_RANGE,
        INVALID_ITERATION_RANGE,
        INVALID_RECTANGLE_NUM,
        INVALID_SIZE,
        INVALID_INVERSION_INDEXES,
        INVALID_INITIAL_SHAPE,
        TYPE_ARRAY_CREATION_ERROR,
        OTHER_ERROR

    } MM_TRIANGULAR_ERROR;
    /**
     * @section ERROR HANDLING FUNCTIONS
     * */
    /**
     * @private show_usage
     * @brief Shows a message about program parameters when they are invalid
     * */
    extern void show_usage( char *argv[], int mpi_rank, char *msg);
    //-- todo: ¿Debería estar en triangular.h?
    /**
     * @private check_results_aid
     * Function used to check that obtained result is correct when a is Identity matrix
     * */
    extern void check_results_aid(Range local_cols, int rows, double *mat_B, double *mat_C,int *check_result);
    //-- Deprecated: doesn't work :(
    void check_results_init_values(Range local_cols, Range current_range,int columns, double *mat_A, double *mat_B, double *mat_C,int *check_result,int rows, int mpi_rank);
    /***
     * @private error_mpi_abort 
     * @param[in] comm          MPI communicator
     * @param[in] error_value   Error numerical code   
     * @param[in] mpi_rank      Current proccessor's MPI numerical identificator
     * @param[in] more_info     Message of error
     * @brief Calling this function will abort comm with exit error_value and will save info of rank and error into an error_file
     * */
    void error_mpi_abort(MPI_Comm comm, MM_TRIANGULAR_ERROR error_value, int mpi_rank, char *more_info);
    char *error_num_name(MM_TRIANGULAR_ERROR error_value, char *info);
    void write_error_filename(char *error_filename, int mpi_rank );

#endif
