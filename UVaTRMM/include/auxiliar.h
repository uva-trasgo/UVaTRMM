/****************************************************************************************
 * @file auxiliar.h     																*
 * @addtogroup mm_triangular															*
 * @brief Auxiliar functions and macro                                                  *
 * @author Arturo González Escribano													*
 * @author Rocío Carratalá Sáez															*
 * @author Maria Inmaculada Santamaria Valenzuela										*
 * @author Yuri Torres de la Sierra														*
 ****************************************************************************************/

#ifndef INCLUDE_AUXILIAR
    #define INCLUDE_AUXILIAR
    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
    #include <math.h>
    #include <mpi.h>
    #include <sys/time.h>

    /****************************************
    * @section MACROS FOR BASIC OPERATIONS  *
    *****************************************/
   /**
    * @private MAX
    * @brief MAX(a,b) is the maximum of a and b
    * */
    #define MAX(a,b) a > b ? a : b;
    /**
    * @private MIN
    * @brief MIN(a,b) is the maximum of a and b
    * */
    #define MIN(a,b) a > b ? b : a;
    /****************************************
     * @section MEMORY HANDLING FUNCTIONS   *
     * **************************************/
    /***
     * @private calloc_abort
     * @brief Equivalent to calloc but adding abort in case memory allocation failed. 
     * @param[in] __nmemb   Number of elements
     * @param[in] __size    size of each element
     * @param[in] name      string identificator for created pointer
     * @param[in] comm      MPI communicator to abort in case of error
     * @param[in] mpi_rank  Current proccessor identificator
     * @return Void pointer of size __nmemb * __size
     * */
    void *calloc_abort(size_t __nmemb, size_t __size, const char *name, MPI_Comm comm, int mpi_rank);
    void freeNotNull(void *ptr);
    /***
     * @private malloc_abort
     * @brief Equivalent to malloc but adding abort in case memory allocation failed. 
     * @param[in] __nmemb   Number of elements
     * @param[in] __size    size of each element
     * @param[in] name      string identificator for created pointer
     * @param[in] comm      MPI communicator to abort in case of error
     * @param[in] mpi_rank  Current proccessor identificator
     * @return Void pointer of size __size
     * */
    void *malloc_abort(size_t __nmemb, size_t __size, const char *name, MPI_Comm comm, int mpi_rank);
    /***
     * @private realloc_abort
     * @private Equivalent to realloc but adding abort in case memory reallocation failed. 
     * @param[in out]   ptr       Pointer to realoc
     * @param[in]       __nmemb   Number of elements
     * @param[in]       __size    size of each element
     * @param[in]       name      string identificator for created pointer
     * @param[in]       comm      MPI communicator to abort in case of error
     * @param[in]       mpi_rank  Current proccessor identificator
     * @return ptr reallocated to new size __size
     * */
    void *realloc_abort(void *ptr, size_t __nmemb, size_t __size, const char *name, MPI_Comm comm, int mpi_rank);
    
    /**
     * @section Timming functions
     * */
    /**
     * @private set_ending_times
     * @brief Reduces times values in order to get final values for their analysis
     * */
    extern void set_ending_times(int mpi_rank,double total_t,double precomp_t,double comm1_t,double comm2_t,double mult_t,double for_t);
    /**
     * @private set_time 
     * @brief Auxiliar function for reducing the number of lines when reading time
     * @param[in] time  Initializated timeval struct
     * @return Milliseconds in double type since time was initialized untill current moment
     * */
    extern double set_time(struct timeval time);
#endif //INCLUDE_AUXILIAR
