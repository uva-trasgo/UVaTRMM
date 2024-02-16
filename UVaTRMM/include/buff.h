/****************************************************************************************
 * @file auxiliar.h     																*
 * @addtogroup mm_triangular															*
 * @brief Auxiliar functions and macro                                                  *
 * @author Arturo González Escribano													*
 * @author Rocío Carratalá Sáez															*
 * @author Maria Inmaculada Santamaria Valenzuela										*
 * @author Yuri Torres de la Sierra														*
 ****************************************************************************************/
#ifndef BUFF_H
    #define BUFF_H
    #include "mkl_range.h"
    /**
     * @section STRUCT DEFINITIONS
     * */
    /**
    * @private Structure Buff
    * Simplifies the use of buffers
    * * Buffer is buffer itself
    * * Size is the number of elements in the buffer
    * */
    typedef struct {
        /*@{*/
    	double *buffer; /**< Double array */
    	int size;       /**< Current buffer size */
        /*@}*/
    } Buff;

    /**
     * @section FUNCTIONS
     * */
    /**
     * @private buff
     * @brief Struct Buff builder
     * @param[in] size      Size of the buff we want to create
     * @param[in] mpi_rank  Current proccessor identificator
     * @return Buff of size "size" or aborts in case of error
     * */
    Buff buff(int size, int mpi_rank);

    //-- Deprecated: copies buffer into array
    int buff__double_arrayCalloc(Buff buffer, double **array, int mpi_rank);
    //-- Deprecated: transforms buffer into matrix double array
    void buff__matrix_double_array(Range rows,int columns,double *mat, Buff buffer );
    //-- todo: Maybe Rows, Columns and Mat should be joined into just one parameter of MatrixRange type
    /**
     * @private set_triangular_buffer_data
     * @param[in] rows      Rows of mat
     * @param[in] columns   Columns of mat
     * @param[in] mat       Matrix given as double array
     * @param[in] mpi_rank  Current proccessor identificator
     * @return Buffer containing triangular under diagonal of mat
     * */
    extern Buff set_triangular_buffer_data(Range rows,int columns,double *mat, int mpi_rank);
#endif
