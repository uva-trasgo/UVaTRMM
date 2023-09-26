#include "mm_triangular.h"

void mm_triangular_triangularMpi_commSplit(MPI_Comms *comms, MPI_Comm main_comm, unsigned int mpi_rank, unsigned int mpi_procs){
    _Bool color;
    unsigned int comm;
    if (mpi_procs > 1){
		comms -> buffer = calloc_abort(mpi_procs, sizeof(MPI_Comm ), "comms", MPI_COMM_WORLD, mpi_rank);
		for(comm = 0; comm < mpi_procs; comm ++){
			comms -> buffer[comm] = MPI_COMM_WORLD;
		}
	    #if DEBUG
	    printf("--> Communicators split\n");
	    fflush(stdout);
	    #endif
	    for ( comm=1; comm<mpi_procs; comm++ ) {
	    //for ( comm=0; comm<mpi_procs-1; comm++ ) {
		    color = (mpi_rank >= comm);
		    #if DEBUG 
		    printf("[%d | %d] color %d\n", mpi_rank, comm, color);
		    if (color==1){
    			printf("[%d] --  Sends to --> [%d] \n", mpi_rank, comm);
	    	}
		    fflush(stdout);
		    #endif 
		    MPI_CHECK( 
    			MPI_Comm_split( 
				    comms -> buffer[comm], 	//Old communicator to split
				    color, 			//New communicator that the calling process is to be assigned to
				    mpi_rank, 		//Index order to use in new communicator
				    &comms -> buffer[comm]   //New communicator handle
			    ) 
		    );//<<--Se queja aqui en dev
	    }
    } else {
		comms -> buffer=NULL;
    }
	int comm_size;
	for (int stage = 1; stage < mpi_procs; stage++){
		MPI_Comm_size( comms -> buffer[stage], &comm_size );
		printf("[%d | %d] comm size %d\n", mpi_rank, stage, comm_size);

	}
}
