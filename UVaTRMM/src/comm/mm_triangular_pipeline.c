/********************************************************************
 * @file mm_triangular_p2p										*
 * @addtogroup mm_triangular										*
 * @brief Definitions for mm_triangular program || Pipeline scheme.	*
 * @author Arturo González Escribano								*
 * @author Rocío Carratalá											*
 * @author Maria Inmaculada Santamaria Valenzuela					*
 * ******************************************************************/
#include "mm_triangular.h"
/***
 * @private comm_func_pipeline_recv
 * recv function for pipeline communication scheme
 * **/
void comm_func_pipeline_recv(
    MatrixRange *mat_A_remote, 
    int stage,
    Requests *requests, 
    int mpi_rank,
    int mpi_procs,
    Map_function f_mapping,
    int rows,
    MPI_Datatype *type_recv1,
    MPI_Datatype *type_recv2,
    Buff *triang_buff_recv,
    Params_comm_shape comm_shape,
    Params_comm_mechanism comm_mechanism,
    int comm_from
){
    #if DEBUG
    printf("[%d <- %d tg %d] recvd\n", mpi_rank, comm_from, stage+1);
    fflush(stdout);
    #endif
    requests->size = 2;
    requests->arr = (MPI_Request*) realloc(requests->arr, requests->size*sizeof(MPI_Request));
	requestsMakeNull(requests);

    set_type_recv(mpi_rank, mat_A_remote->rows, mat_A_remote -> columns.size, type_recv1, type_recv2, triang_buff_recv,comm_shape, comm_mechanism);
    switch (comm_shape)
    {
        case SHAPE_FULL:
            mpi_irecv_splitted(
                mat_A_remote -> mat,mat_A_remote -> rows.size * mat_A_remote -> columns.size, 
                comm_from, MPI_COMM_WORLD, requests,mpi_rank, MPI_DOUBLE,
                0, 1, comm_from //stage+1
            ); 
            #if DEBUG
            printf("[%d | %d] After Recv Full %d\n", mpi_rank, stage, requests->size);
            fflush(stdout);
            #endif 
        break;
	    case SHAPE_TRIANG_2:
		    mpi_irecv_splitted(
                mat_A_remote -> mat, 1, comm_from, MPI_COMM_WORLD, 
                requests, mpi_rank, MPI_DOUBLE, 0, 2, comm_from //stage+1
            );
            if ( comm_mechanism == MECH_BUFFER ) {
                mpi_irecv_splitted(
                    triang_buff_recv -> buffer,triang_buff_recv -> size, comm_from,
                    MPI_COMM_WORLD, requests, mpi_rank,MPI_DOUBLE, 1, 2, comm_from //stage+1
                );
            } else {
                mpi_irecv_splitted( 
                    mat_A_remote -> mat, 1, comm_from, MPI_COMM_WORLD, 
                    requests, mpi_rank, *type_recv2, 1, 2, comm_from //stage+1 
                );
		    }
        break;
        default:
            mpi_irecv_splitted(
                mat_A_remote -> mat, 1, comm_from, MPI_COMM_WORLD, 
                requests, mpi_rank,  *type_recv1, 0, 1, comm_from // stage+1
            );
        break;
	}
}


/***
 * @private comm_func_pipeline_send
 * send function for pipeline communication scheme
 * **/
void comm_func_pipeline_send(
	MatrixRange mat_A, 
	int stage,
	Requests *requests,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatype *type_send1,
	MPI_Datatype *type_send2,
	Buff *triang_buff_send,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism,
    int comm_to
){
    #if DEBUG
    printf("[%d -> %d tg %d] send_size (s%ld,b%ld)\n", mpi_rank, comm_to, stage+1, mat_A.rows.size, mat_A.rows.begin);
    fflush(stdout);
    #endif
    requests->size = 2;
	requests->arr = (MPI_Request*) realloc(requests->arr, requests->size*sizeof(MPI_Request));
	requestsMakeNull(requests);
    set_type_send(mpi_rank, type_send1,type_send2, triang_buff_send, &mat_A, comm_shape, comm_mechanism, FUNC_P2P);
    
	switch (comm_shape)
	{
	case SHAPE_FULL:
        #if DEBUG
        printf("[%d | %d] Before Send Full\n", mpi_rank, stage);
        fflush(stdout);
        #endif 
        mpi_isend_splitted(
            mat_A.mat, mat_A.rows.size * mat_A.columns.size, comm_to, MPI_COMM_WORLD, 
            requests, mpi_rank, MPI_DOUBLE, 0, 1, mpi_rank, 1 //stage+1, 1
        );
        #if DEBUG
        printf("[%d | %d] After Send Full %d\n", mpi_rank, stage, requests->size);
        fflush(stdout);
        #endif 
		break;
    case SHAPE_TRIANG_2:
        // Rectangle part
		mpi_isend_splitted(
            mat_A.mat, 1, comm_to, MPI_COMM_WORLD, 
            requests, mpi_rank,  *type_send1,0, 2, mpi_rank, 1 
        );
		// Triangle part
		if ( comm_mechanism == MECH_TYPE ) {
			mpi_isend_splitted( 
                mat_A.mat, 1, comm_to, MPI_COMM_WORLD, 
                requests, mpi_rank, *type_send2, 1, 2, mpi_rank, 1 //0
            );
		} else {
		    *triang_buff_send = set_triangular_buffer_data(mat_A.rows, mat_A.columns.size, mat_A.mat, mpi_rank);
			mpi_isend_splitted( 
                triang_buff_send -> buffer, triang_buff_send -> size, comm_to, MPI_COMM_WORLD, 
                requests, mpi_rank, MPI_DOUBLE, 1, 2, mpi_rank, 1
            );
		}
        break;
	default:
        mpi_isend_splitted(
            mat_A. mat, 1, comm_to, MPI_COMM_WORLD, 
            requests, mpi_rank, *type_send1, 0, 1, mpi_rank, 1 //stage+1, 0
        );
		break;
    }			
}



void comm_func_pipeline_recvArray(
    MatrixRange *mat_A_remote, 
    int stage,
    Requests *requests, 
    int mpi_rank,
    int mpi_procs,
    Map_function f_mapping,
    int rows,
    MPI_Datatypes *type_recv1,
    MPI_Datatypes *type_recv2,
    Buff *triang_buff_recv,
    Params_comm_shape comm_shape,
    Params_comm_mechanism comm_mechanism,
    int comm_from, unsigned int quotient, MPI_Comm comm
){
    #if DEBUG
    printf("[%d <- %d tg %d] recvd\n", mpi_rank, comm_from, stage+1);
    fflush(stdout);
    #endif
    requests->size = 2;
    requests->arr = (MPI_Request*) realloc(requests->arr, requests->size*sizeof(MPI_Request));
	requestsMakeNull(requests);

    set_type_recvArray(mpi_rank, mat_A_remote->rows, mat_A_remote -> columns.size, type_recv1, type_recv2, triang_buff_recv,comm_shape, comm_mechanism, quotient, comm, stage);
    switch (comm_shape)
    {
        case SHAPE_FULL:
            mpi_irecv_splitted(
                mat_A_remote -> mat,mat_A_remote -> rows.size * mat_A_remote -> columns.size, 
                comm_from, MPI_COMM_WORLD, requests,mpi_rank, MPI_DOUBLE,
                0, 1, comm_from //stage+1
            ); 
            #if DEBUG
            printf("[%d | %d] After Recv Full %d\n", mpi_rank, stage, requests->size);
            fflush(stdout);
            #endif 
        break;
	    case SHAPE_TRIANG_2:
		    mpi_irecv_array(
                mat_A_remote -> mat, comm_from, MPI_COMM_WORLD, 
                requests, mpi_rank, *type_recv1, 0, 2, comm_from //stage+1
            );
            if ( comm_mechanism == MECH_BUFFER ) {
                mpi_irecv_splitted(
                    triang_buff_recv -> buffer,triang_buff_recv -> size, comm_from,
                    MPI_COMM_WORLD, requests, mpi_rank,MPI_DOUBLE, 1, 2, comm_from //stage+1
                );
            } else {
                mpi_irecv_array( 
                    mat_A_remote -> mat, comm_from, MPI_COMM_WORLD, 
                    requests, mpi_rank, *type_recv2, 1, 2, comm_from //stage+1 
                );
		    }
        break;
        default:
            mpi_irecv_array(
                mat_A_remote -> mat,comm_from, MPI_COMM_WORLD,
                requests, mpi_rank,  *type_recv1, 0, 1, comm_from // stage+1
            );
        break;
	}
}


/***
 * @private comm_func_pipeline_send
 * send function for pipeline communication scheme
 * **/
void comm_func_pipeline_sendArray(
	MatrixRange mat_A, 
	int stage,
	Requests *requests,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatypes *type_send1,
	MPI_Datatypes *type_send2,
	Buff *triang_buff_send,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism,
    int comm_to,unsigned int quotient, MPI_Comm comm
){
    #if DEBUG
    printf("[%d -> %d tg %d] send_size (s%ld,b%ld)\n", mpi_rank, comm_to, stage+1, mat_A.rows.size, mat_A.rows.begin);
    fflush(stdout);
    #endif
    requests->size = 2;
	requests->arr = (MPI_Request*) realloc(requests->arr, requests->size*sizeof(MPI_Request));
	requestsMakeNull(requests);
    *triang_buff_send=set_type_sendArray(mpi_rank, type_send1,type_send2, triang_buff_send, &mat_A, comm_shape, comm_mechanism, FUNC_P2P, quotient, comm, stage);
    
	switch (comm_shape)
	{
	case SHAPE_FULL:
        #if DEBUG
        printf("[%d | %d] Before Send Full\n", mpi_rank, stage);
        fflush(stdout);
        #endif 
        mpi_isend_splitted(
            mat_A.mat, mat_A.rows.size * mat_A.columns.size, comm_to, MPI_COMM_WORLD, 
            requests, mpi_rank, MPI_DOUBLE, 0, 1, mpi_rank, 1 //stage+1, 1
        );
        #if DEBUG
        printf("[%d | %d] After Send Full %d\n", mpi_rank, stage, requests->size);
        fflush(stdout);
        #endif 
		break;
    case SHAPE_TRIANG_2:
        // Rectangle part
		mpi_isend_array(
            mat_A.mat, comm_to, MPI_COMM_WORLD,
            requests, mpi_rank,  *type_send1,0, 2, mpi_rank, 1 
        );
		// Triangle part
		if ( comm_mechanism == MECH_TYPE ) {
			mpi_isend_array( 
                mat_A.mat, comm_to, MPI_COMM_WORLD, 
                requests, mpi_rank, *type_send2, 1, 2, mpi_rank, 1 //0
            );
		} else {
		    *triang_buff_send = set_triangular_buffer_data(mat_A.rows, mat_A.columns.size, mat_A.mat, mpi_rank);
			mpi_isend_splitted( 
                triang_buff_send -> buffer, triang_buff_send -> size, comm_to, comm, 
                requests, mpi_rank, MPI_DOUBLE, 1, 2, mpi_rank, 1
            );
		}
        break;
	default:
        mpi_isend_array(
            mat_A. mat, comm_to, MPI_COMM_WORLD, 
            requests, mpi_rank, *type_send1, 0, 1, mpi_rank, 1//stage+1, 0
        );
		break;
    }			
}
