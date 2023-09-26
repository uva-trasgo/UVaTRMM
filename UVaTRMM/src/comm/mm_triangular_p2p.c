/********************************************************************
 * @file mm_triangular_p2p										*
 * @addtogroup mm_triangular										*
 * @brief Definitions for mm_triangular program || Skew scheme.		*
 * @author Arturo González Escribano								*
 * @author Rocío Carratalá											*
 * @author Maria Inmaculada Santamaria Valenzuela					*
 * ******************************************************************/

#include "mm_triangular.h"

/***
 * @private comm_func_p2p_recv
 * Receiving function for owner_skew mode
 * **/
void comm_func_p2p_recv(
    int comm_from,
	MatrixRange mat_A, 
	MatrixRange *current_A,
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
	Params_comm_mechanism comm_mechanism
){
	requests->size= 2;
	requests->arr = (MPI_Request*) realloc(requests->arr, requests->size*sizeof(MPI_Request));
	requestsMakeNull(requests);

    set_type_recv(mpi_rank, current_A->rows, mat_A.columns.size, type_recv1, type_recv2, triang_buff_recv, comm_shape,comm_mechanism);
	switch (comm_shape)
    {
	    case SHAPE_FULL:
    	    mpi_irecv_splitted(current_A -> mat, current_A -> rows.size * current_A -> columns.size, comm_from, MPI_COMM_WORLD, requests, mpi_rank, MPI_DOUBLE, 0, 2, stage+1);
        	break;
		case SHAPE_TRIANG_2: 
				mpi_irecv_splitted(current_A -> mat, 1, comm_from, MPI_COMM_WORLD, requests, mpi_rank, *type_recv1, 0, 2, stage+1);
			// Triangle part
			if ( comm_mechanism == MECH_BUFFER ) {
				mpi_irecv_splitted( triang_buff_recv -> buffer, triang_buff_recv -> size, comm_from,  MPI_COMM_WORLD, requests, mpi_rank, MPI_DOUBLE, 1, 2, stage+1 );
			} else {
				mpi_irecv_splitted( current_A -> mat, 1, comm_from, MPI_COMM_WORLD, requests, mpi_rank, *type_recv2, 1, 2, stage+1 );
			}
			break;
    	default:
        	mpi_irecv_splitted(current_A -> mat, 1, comm_from, MPI_COMM_WORLD, requests, mpi_rank, *type_recv1, 0, 2, stage+1);
        break;
	}
}

/***
 * @private comm_func_p2p_recv
 * Sending function for owner_skew mode
 * **/
void comm_func_p2p_send(
    int comm_to,
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
	Params_comm_mechanism comm_mechanism
){
	requests->size= 2;
	requests->arr = (MPI_Request*) realloc(requests->arr, requests->size*sizeof(MPI_Request));
	requestsMakeNull(requests);
	set_type_send(mpi_rank,type_send1, type_send2,triang_buff_send, &mat_A, comm_shape,comm_mechanism, FUNC_P2P);
	switch (comm_shape)
	{
	case SHAPE_FULL:
		mpi_isend_splitted(mat_A.mat, mat_A.rows.size * mat_A.columns.size, comm_to, MPI_COMM_WORLD, requests, mpi_rank,  MPI_DOUBLE, 0, 2, stage+1,1 );
		break;	
	case SHAPE_TRIANG_2:
		mpi_isend_splitted(mat_A.mat, 1, comm_to, MPI_COMM_WORLD, requests, mpi_rank, *type_send1, 0, 2, stage+1,1);
		if ( comm_mechanism == MECH_BUFFER ) {
			mpi_isend_splitted( triang_buff_send -> buffer,triang_buff_send -> size, comm_to, MPI_COMM_WORLD, requests, mpi_rank, MPI_DOUBLE, 1, 2, stage+1,1);
		} else {
			mpi_isend_splitted( mat_A.mat, 1, comm_to, MPI_COMM_WORLD, requests, mpi_rank, *type_send2, 1, 2, stage+1,1);
		}
		break;
	default:
		mpi_isend_splitted(mat_A.mat, 1, comm_to, MPI_COMM_WORLD, requests, mpi_rank,  *type_send1,0, 2, stage+1,1);
		break;
	}
}

void comm_func_p2p_sendArray(
    int comm_to,
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
	unsigned int quotient, 
	MPI_Comm comm 
){
	requests->size= 2;
	requests->arr = (MPI_Request*) realloc(requests->arr, requests->size*sizeof(MPI_Request));
	requestsMakeNull(requests);
	*triang_buff_send=set_type_sendArray(mpi_rank,type_send1, type_send2,triang_buff_send, &mat_A, comm_shape,comm_mechanism, FUNC_P2P, quotient, comm, stage);
	switch (comm_shape)
	{
	case SHAPE_FULL:
		mpi_isend_splitted(mat_A.mat, mat_A.rows.size * mat_A.columns.size, comm_to, MPI_COMM_WORLD, requests, mpi_rank,  MPI_DOUBLE, 0, 2, stage+1,1);
		break;	
	case SHAPE_TRIANG_2:
		mpi_isend_array(mat_A.mat, comm_to, MPI_COMM_WORLD, requests, mpi_rank, *type_send1, 0, 2, stage+1,  1);
		if ( comm_mechanism == MECH_BUFFER ) {
			mpi_isend_splitted( triang_buff_send -> buffer,triang_buff_send -> size, comm_to, MPI_COMM_WORLD, requests, mpi_rank, MPI_DOUBLE, 1, 2, stage+1,1);
		} else {
			mpi_isend_array( mat_A.mat,comm_to, MPI_COMM_WORLD, requests, mpi_rank, *type_send2, 1, 2, stage+1,1);
		}
		break;
	default:
		mpi_isend_array(mat_A.mat, comm_to, MPI_COMM_WORLD, requests, mpi_rank,  *type_send1,0, 2, stage+1,1);
		break;
	}
}

void comm_func_p2p_recvArray(
    int comm_from,
	MatrixRange mat_A, 
	MatrixRange *current_A,
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
	unsigned int quotient, 
	MPI_Comm comm 
){
	requests->size= 2;
	requests->arr = (MPI_Request*) realloc(requests->arr, requests->size*sizeof(MPI_Request));
	requestsMakeNull(requests);

    set_type_recvArray(mpi_rank, current_A->rows, mat_A.columns.size, type_recv1, type_recv2, triang_buff_recv, comm_shape,comm_mechanism, quotient, comm, stage);
	switch (comm_shape)
    {
	    case SHAPE_FULL:
    	    mpi_irecv_splitted(current_A -> mat, current_A -> rows.size * current_A -> columns.size, comm_from, MPI_COMM_WORLD, requests, mpi_rank, MPI_DOUBLE, 0, 2, stage+1);
        	break;
		case SHAPE_TRIANG_2: 
				mpi_irecv_array(current_A -> mat, comm_from, MPI_COMM_WORLD, requests, mpi_rank, *type_recv1, 0, 2, stage+1);
			// Triangle part
			if ( comm_mechanism == MECH_BUFFER ) {
				mpi_irecv_splitted( triang_buff_recv -> buffer, triang_buff_recv -> size, comm_from,  MPI_COMM_WORLD, requests, mpi_rank, MPI_DOUBLE, 1, 2, stage+1 );
			} else {
				mpi_irecv_array( current_A -> mat, comm_from, MPI_COMM_WORLD, requests, mpi_rank, *type_recv2, 1, 2, stage+1 );
			}
			break;
    	default:
        	mpi_irecv_array(current_A -> mat, comm_from, MPI_COMM_WORLD, requests, mpi_rank, *type_recv1, 0, 2, stage+1);
        break;
	}
}
