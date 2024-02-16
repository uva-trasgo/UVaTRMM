#ifndef MKL_AUXILIAR
#define MKL_AUXILIAR
#include "mkl.h"


/**
 * MKL MACROS
 **/
/* MPI check error */
#define MPI_CHECK( a )	{ int error_code = a; if ( error_code != MPI_SUCCESS ) { int msg_len; char msg[1024]; MPI_Error_string( error_code, msg, &msg_len); fprintf( stderr, "\nRuntime error in line %d: Error message: %s\n\n", __LINE__, msg ); MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE ); } }

/**
 * @private Structure Requests
 * Simplifies the use of MPI_Request buffers
 * * arr is a MPI_Request buffer
 * * Size is the number of elements in the buffer
 * */
typedef struct {
	MPI_Request *arr;
	int 		size;
} Requests;

/**
 * @private Structure Requests
 * Simplifies the use of MPI_Request buffers
 * * arr is a MPI_Request buffer
 * * Size is the number of elements in the buffer
 * */
typedef struct {
	Requests *arr;
	int 		size;
} Requestss;



/**
 * @private Structure Statuses
 * Simplifies the use of MPI_Request buffers
 * * arr is a MPI_Request buffer
 * * Size is the number of elements in the buffer
 * */
typedef struct {
	MPI_Status  *arr;
	int 		size;
} Statuses;

typedef struct {
	Statuses *arr;
	int 		size;
} Statusess;

/**
 * Functions
 * */
Requests requestsCalloc_abort(size_t size); //, MPI_Comm comm, int mpi_rank);
Requests requestsMalloc_abort(size_t size); //, MPI_Comm comm, int mpi_rank);
Statuses statusesCalloc_abort(size_t size); //, MPI_Comm comm, int mpi_rank);
void requestsMakeNull(Requests *requests);

#endif //MKL_AUXILIAR
