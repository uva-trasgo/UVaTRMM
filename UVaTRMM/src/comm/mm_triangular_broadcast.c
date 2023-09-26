/********************************************************************
 * @file mm_triangular_broadcast														*
 * @addtogroup mm_triangular															*
 * @brief Functions definitions for mm_triangular program needed for owner_broad scheme.*
 * @author Arturo González Escribano													*
 * @author Rocío Carratalá Sáez															*
 * @author Maria Inmaculada Santamaria Valenzuela										*
 * @author Yuri Torres de la Sierra														*
 * **************************************************************************************/
#include "mm_triangular.h"

/***
 * @private comm_func_bcast_send
 * Broadcast function for sending mat_A
 * Makes current_A = mat_A //--> ¿Debería? ¿No sobra?
 * */
void comm_func_bcast_send(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_request,
	Requests *broad_request2,
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
	*triang_buff_send = set_type_send(mpi_rank, type_send1, type_send2, triang_buff_send, &mat_A, comm_shape,comm_mechanism, FUNC_BROADCAST);
	#if DEBUG
	matrix_rangePrintf_rank_stage(mat_A, "Asending", mpi_rank, stage, 0);
	#endif 
	switch (comm_shape)
	{
	case SHAPE_FULL:
		mpi_broadcast_splitted(mat_A.mat, mat_A.rows.size * mat_A.columns.size, stage, MPI_COMM_WORLD, broad_request, mpi_rank, MPI_DOUBLE);
		break;
	case SHAPE_TRIANG_2:
		mpi_broadcast_splitted(mat_A.mat, 1, stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_send1);
		if ( comm_mechanism == MECH_BUFFER ) {
			mpi_broadcast_splitted(triang_buff_send -> buffer, triang_buff_send -> size, stage,  MPI_COMM_WORLD, broad_request2, mpi_rank, MPI_DOUBLE);
		} else {
			mpi_broadcast_splitted(mat_A.mat, 1, stage,  MPI_COMM_WORLD, broad_request2, mpi_rank, *type_send2);
		}
		break;
		default:
			mpi_broadcast_splitted(mat_A.mat, 1, stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_send1);
		break;
	}
	current_A -> mat = mat_A.mat;//realloc(mat_A, local_rows.size*columns*sizeof(double));
	current_A -> rows = mat_A.rows;
	current_A -> columns = mat_A.columns;
	#if DEBUG
	printf("| rank %d | SEND AFTER > | SIZE %ld | BUFFERSIZE %ld\n", stage, mat_A.rows.size, mat_A.rows.size*mat_A.columns.size*sizeof(double));
	printf_matrix_ranges_rank_stage(current_A ->  mat, current_A -> rows.begin, current_A -> rows.size+current_A -> rows.begin, 0, mat_A.columns.size,"CSA", mpi_rank, stage, 0,0);
	fflush(stdout);
	printf("After Broadcast shape full send | %d | %d\n", mpi_rank, stage);
	fflush(stdout);
	#endif 
}

/***
 * @private comm_func_bcast_recv
 * Broadcast function for receiving remote_A from mat_A and filling up current_A
 * **/
void comm_func_bcast_recv(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_request,
	Requests *broad_request2,
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
	//Allocate from remote part
	MatrixRange remote_A;
	remote_A.columns 	= mat_A.columns;
	remote_A.rows 		= f_mapping( rows, mpi_procs, stage, mpi_rank);
	matrix_rangeMat_calloc_abort(&remote_A, "remote_A", MPI_COMM_WORLD, mpi_rank);
	set_type_recv(mpi_rank, remote_A.rows, remote_A.columns.size, type_recv1, type_recv2, triang_buff_recv,comm_shape,comm_mechanism);
	switch (comm_shape){
		case SHAPE_FULL: 
			mpi_broadcast_splitted(
				remote_A.mat, 
				remote_A.rows.size * remote_A.columns.size, 
				stage, MPI_COMM_WORLD, broad_request, mpi_rank, MPI_DOUBLE
			);
		break;
		case SHAPE_TRIANG_2:
			mpi_broadcast_splitted(remote_A.mat, 1, stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_recv1);
			if ( comm_mechanism == MECH_BUFFER ) {
				mpi_broadcast_splitted(triang_buff_recv -> buffer, triang_buff_recv -> size, stage, MPI_COMM_WORLD, broad_request2, mpi_rank, MPI_DOUBLE);
			} else {
				mpi_broadcast_splitted(remote_A.mat, 1, stage, MPI_COMM_WORLD, broad_request2, mpi_rank, *type_recv2);
			}
		break;
		default:
			mpi_broadcast_splitted(remote_A.mat, 1,  stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_recv1);
		break;
	}
	#if DEBUG
	printf("%d | %d | Before wait\n", mpi_rank, stage);
	fflush(stdout);
	#endif
	/*MPI_CHECK( MPI_Wait(broad_request, MPI_STATUS_IGNORE ) );
	if (comm_shape == SHAPE_TRIANG_2){
		MPI_CHECK( MPI_Wait(broad_request2, MPI_STATUS_IGNORE ) );
	}
	#if DEBUG
	printf("%d | %d | After wait\n", mpi_rank, stage);
    fflush(stdout);
	#endif 
	if (comm_shape == SHAPE_TRIANG_2 && comm_mechanism == MECH_BUFFER){
		//set_matrix_data_from_buffer(remote_A.rows, remote_A.columns.size, remote_A.mat, *triang_buff_recv);
		buff__matrix_range(*triang_buff_recv, &remote_A);
	}*/
	current_A -> mat = remote_A.mat;
	current_A -> columns = remote_A.columns;
	current_A -> rows = remote_A.rows;
}
/***
 * @private comm_func_bcast
 * Send/Receive function for broadcast mode
 * **/
void comm_func_bcast(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_requests,
	Requests *broad_requests2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatype *type_send1,
	MPI_Datatype *type_send2,
	Buff *triang_buff_send,
	MPI_Datatype *type_recv1,
	MPI_Datatype *type_recv2,
	Buff *triang_buff_recv,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism
){
	if (stage == mpi_rank){
		comm_func_bcast_send(mat_A, current_A, stage, broad_requests, broad_requests2, mpi_rank, mpi_procs, f_mapping, rows, type_send1, type_send2, triang_buff_send, comm_shape, comm_mechanism);
	} else {
		comm_func_bcast_recv(mat_A, current_A, stage, broad_requests, broad_requests2, mpi_rank, mpi_procs, f_mapping, rows, type_recv1, type_recv2, triang_buff_recv, comm_shape, comm_mechanism);
	}
}

/***
 * @private comm_func_bcast_send
 * Broadcast function for sending mat_A
 * Makes current_A = mat_A //--> ¿Debería? ¿No sobra?
 * */
void comm_func_bcast_send_quotient(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_request,
	Requests *broad_request2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatype *type_send1,
	MPI_Datatype *type_send2,
	Buff *triang_buff_send,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism, unsigned int quotient
){
	*triang_buff_send = set_type_send(mpi_rank, type_send1, type_send2, triang_buff_send, &mat_A, comm_shape,comm_mechanism, FUNC_BROADCAST);
	#if DEBUG
	matrix_rangePrintf_rank_stage(mat_A, "Asending", mpi_rank, stage, 0);
	#endif 
	switch (comm_shape)
	{
	case SHAPE_FULL:
		mpi_broadcast_splitted_quotient(mat_A.mat, mat_A.rows.size * mat_A.columns.size, stage, MPI_COMM_WORLD, broad_request, mpi_rank, MPI_DOUBLE, quotient);
		break;
	case SHAPE_TRIANG_2:
		mpi_broadcast_splitted_quotient(mat_A.mat, 1, stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_send1, quotient);
		if ( comm_mechanism == MECH_BUFFER ) {
			mpi_broadcast_splitted_quotient(triang_buff_send -> buffer, triang_buff_send -> size, stage,  MPI_COMM_WORLD, broad_request2, mpi_rank, MPI_DOUBLE, quotient);
		} else {
			mpi_broadcast_splitted_quotient(mat_A.mat, 1, stage,  MPI_COMM_WORLD, broad_request2, mpi_rank, *type_send2,quotient);
		}
		break;
		default:
			mpi_broadcast_splitted_quotient(mat_A.mat, 1, stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_send1,quotient);
		break;
	}
	current_A -> mat = mat_A.mat;//realloc(mat_A, local_rows.size*columns*sizeof(double));
	current_A -> rows = mat_A.rows;
	current_A -> columns = mat_A.columns;
	#if DEBUG
	printf("| rank %d | SEND AFTER > | SIZE %ld | BUFFERSIZE %ld\n", stage, mat_A.rows.size, mat_A.rows.size*mat_A.columns.size*sizeof(double));
	printf_matrix_ranges_rank_stage(current_A ->  mat, current_A -> rows.begin, current_A -> rows.size+current_A -> rows.begin, 0, mat_A.columns.size,"CSA", mpi_rank, stage, 0,0);
	fflush(stdout);
	printf("After Broadcast shape full send | %d | %d\n", mpi_rank, stage);
	fflush(stdout);
	#endif 
}

/***
 * @private comm_func_bcast_recv
 * Broadcast function for receiving remote_A from mat_A and filling up current_A
 * **/
void comm_func_bcast_recv_quotient(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_request,
	Requests *broad_request2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatype *type_recv1,
	MPI_Datatype *type_recv2,
	Buff *triang_buff_recv,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism, unsigned int quotient
){
	//Allocate from remote part
	MatrixRange remote_A;
	remote_A.columns 	= mat_A.columns;
	remote_A.rows 		= f_mapping( rows, mpi_procs, stage, mpi_rank);
	matrix_rangeMat_calloc_abort(&remote_A, "remote_A", MPI_COMM_WORLD, mpi_rank);
	set_type_recv(mpi_rank, remote_A.rows, remote_A.columns.size, type_recv1, type_recv2, triang_buff_recv,comm_shape,comm_mechanism);
	switch (comm_shape){
		case SHAPE_FULL: 
			mpi_broadcast_splitted_quotient(
				remote_A.mat, 
				remote_A.rows.size * remote_A.columns.size, 
				stage, MPI_COMM_WORLD, broad_request, mpi_rank, MPI_DOUBLE, quotient
			);
		break;
		case SHAPE_TRIANG_2:
			mpi_broadcast_splitted_quotient(remote_A.mat, 1, stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_recv1, quotient);
			if ( comm_mechanism == MECH_BUFFER ) {
				mpi_broadcast_splitted_quotient(triang_buff_recv -> buffer, triang_buff_recv -> size, stage, MPI_COMM_WORLD, broad_request2, mpi_rank, MPI_DOUBLE, quotient);
			} else {
				mpi_broadcast_splitted_quotient(remote_A.mat, 1, stage, MPI_COMM_WORLD, broad_request2, mpi_rank, *type_recv2, quotient);
			}
		break;
		default:
			mpi_broadcast_splitted_quotient(remote_A.mat, 1,  stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_recv1, quotient);
		break;
	}
	#if DEBUG
	printf("%d | %d | Before wait\n", mpi_rank, stage);
	fflush(stdout);
	#endif
	/*MPI_CHECK( MPI_Wait(broad_request, MPI_STATUS_IGNORE ) );
	if (comm_shape == SHAPE_TRIANG_2){
		MPI_CHECK( MPI_Wait(broad_request2, MPI_STATUS_IGNORE ) );
	}
	#if DEBUG
	printf("%d | %d | After wait\n", mpi_rank, stage);
    fflush(stdout);
	#endif 
	if (comm_shape == SHAPE_TRIANG_2 && comm_mechanism == MECH_BUFFER){
		//set_matrix_data_from_buffer(remote_A.rows, remote_A.columns.size, remote_A.mat, *triang_buff_recv);
		buff__matrix_range(*triang_buff_recv, &remote_A);
	}*/
	current_A -> mat = remote_A.mat;
	current_A -> columns = remote_A.columns;
	current_A -> rows = remote_A.rows;
}
/***
 * @private comm_func_bcast
 * Send/Receive function for broadcast mode
 * **/
void comm_func_bcast_quotient(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_requests,
	Requests *broad_requests2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatype *type_send1,
	MPI_Datatype *type_send2,
	Buff *triang_buff_send,
	MPI_Datatype *type_recv1,
	MPI_Datatype *type_recv2,
	Buff *triang_buff_recv,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism, unsigned int quotient
){
	if (stage == mpi_rank){
		comm_func_bcast_send_quotient(mat_A, current_A, stage, broad_requests, broad_requests2, mpi_rank, mpi_procs, f_mapping, rows, type_send1, type_send2, triang_buff_send, comm_shape, comm_mechanism, quotient);
	} else {
		comm_func_bcast_recv_quotient(mat_A, current_A, stage, broad_requests, broad_requests2, mpi_rank, mpi_procs, f_mapping, rows, type_recv1, type_recv2, triang_buff_recv, comm_shape, comm_mechanism, quotient);
	}
}


/***
 * @private comm_func_bcast_send
 * Broadcast function for sending mat_A
 * Makes current_A = mat_A //--> ¿Debería? ¿No sobra?
 * */
void comm_func_bcast_sendArray(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_request,
	Requests *broad_request2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatypes *type_send1,
	MPI_Datatypes *type_send2,
	Buff *triang_buff_send,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism,
	int block_limit_quotient_denominator,
	MPI_Comm comm
){
	*triang_buff_send = set_type_sendArray(mpi_rank, type_send1, type_send2, triang_buff_send, &mat_A, comm_shape,comm_mechanism, FUNC_BROADCAST, block_limit_quotient_denominator, comm, stage);
	#if DEBUG
	matrix_rangePrintf_rank_stage(mat_A, "Asending", mpi_rank, stage, 0);
	#endif 
	switch (comm_shape)
	{
	case SHAPE_FULL:
		mpi_broadcast_splitted(mat_A.mat, mat_A.rows.size * mat_A.columns.size, stage, MPI_COMM_WORLD, broad_request, mpi_rank, MPI_DOUBLE);
		break;
	case SHAPE_TRIANG_2:
		mpi_broadcast_splitted_quotientArray(mat_A.mat, stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_send1);
		if ( comm_mechanism == MECH_BUFFER ) {
			mpi_broadcast_splitted(triang_buff_send -> buffer, triang_buff_send -> size, stage,  MPI_COMM_WORLD, broad_request2, mpi_rank, MPI_DOUBLE);
		} else {
			mpi_broadcast_splitted_quotientArray(mat_A.mat, stage,  MPI_COMM_WORLD, broad_request2, mpi_rank, *type_send2);
		}
		break;
		default:
			mpi_broadcast_splitted_quotientArray(mat_A.mat, stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_send1);
		break;
	}
	current_A -> mat = mat_A.mat;//realloc(mat_A, local_rows.size*columns*sizeof(double));
	current_A -> rows = mat_A.rows;
	current_A -> columns = mat_A.columns;
	#if DEBUG
	printf("| rank %d | SEND AFTER > | SIZE %ld | BUFFERSIZE %ld\n", stage, mat_A.rows.size, mat_A.rows.size*mat_A.columns.size*sizeof(double));
	printf_matrix_ranges_rank_stage(current_A ->  mat, current_A -> rows.begin, current_A -> rows.size+current_A -> rows.begin, 0, mat_A.columns.size,"CSA", mpi_rank, stage, 0,0);
	fflush(stdout);
	printf("After Broadcast shape full send | %d | %d\n", mpi_rank, stage);
	fflush(stdout);
	#endif 
}

/***
 * @private comm_func_bcast_recv
 * Broadcast function for receiving remote_A from mat_A and filling up current_A
 * **/
void comm_func_bcast_recvArray(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_request,
	Requests *broad_request2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatypes *type_recv1,
	MPI_Datatypes *type_recv2,
	Buff *triang_buff_recv,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism,
	int block_limit_quotient_denominator,
	MPI_Comm comm
){
	//Allocate from remote part
	MatrixRange remote_A;
	remote_A.columns 	= mat_A.columns;
	remote_A.rows 		= f_mapping( rows, mpi_procs, stage, mpi_rank);
	matrix_rangeMat_calloc_abort(&remote_A, "remote_A", MPI_COMM_WORLD, mpi_rank);
	set_type_recvArray(mpi_rank, remote_A.rows, remote_A.columns.size, type_recv1, type_recv2, triang_buff_recv,comm_shape,comm_mechanism, block_limit_quotient_denominator, comm, stage);
	switch (comm_shape){
		case SHAPE_FULL: 
			mpi_broadcast_splitted(
				remote_A.mat, 
				remote_A.rows.size * remote_A.columns.size, 
				stage, MPI_COMM_WORLD, broad_request, mpi_rank, MPI_DOUBLE
			);
		break;
		case SHAPE_TRIANG_2:
			mpi_broadcast_splitted_quotientArray(remote_A.mat, stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_recv1);
			if ( comm_mechanism == MECH_BUFFER ) {
				mpi_broadcast_splitted(triang_buff_recv -> buffer, triang_buff_recv -> size, stage, MPI_COMM_WORLD, broad_request2, mpi_rank, MPI_DOUBLE);
			} else {
				mpi_broadcast_splitted_quotientArray(remote_A.mat, stage, MPI_COMM_WORLD, broad_request2, mpi_rank, *type_recv2);
			}
		break;
		default:
			mpi_broadcast_splitted_quotientArray(remote_A.mat,  stage, MPI_COMM_WORLD, broad_request, mpi_rank, *type_recv1);
		break;
	}
	#if DEBUG
	printf("%d | %d | Before wait\n", mpi_rank, stage);
	fflush(stdout);
	#endif
	/*MPI_CHECK( MPI_Wait(broad_request, MPI_STATUS_IGNORE ) );
	if (comm_shape == SHAPE_TRIANG_2){
		MPI_CHECK( MPI_Wait(broad_request2, MPI_STATUS_IGNORE ) );
	}
	#if DEBUG
	printf("%d | %d | After wait\n", mpi_rank, stage);
    fflush(stdout);
	#endif 
	if (comm_shape == SHAPE_TRIANG_2 && comm_mechanism == MECH_BUFFER){
		//set_matrix_data_from_buffer(remote_A.rows, remote_A.columns.size, remote_A.mat, *triang_buff_recv);
		buff__matrix_range(*triang_buff_recv, &remote_A);
	}*/
	current_A -> mat = remote_A.mat;
	current_A -> columns = remote_A.columns;
	current_A -> rows = remote_A.rows;
}

void comm_func_bcast_quotientArray(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_requests,
	Requests *broad_requests2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatypes *type_send1,
	MPI_Datatypes *type_send2,
	Buff *triang_buff_send,
	MPI_Datatypes *type_recv1,
	MPI_Datatypes *type_recv2,
	Buff *triang_buff_recv,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism, unsigned int quotient, 
	MPI_Comm comm 
){
	if (stage == mpi_rank){
		comm_func_bcast_sendArray(mat_A, current_A, stage, broad_requests, broad_requests2, mpi_rank, mpi_procs, f_mapping, rows, type_send1, type_send2, triang_buff_send, comm_shape, comm_mechanism, quotient, comm);
	} else {
		comm_func_bcast_recvArray(mat_A, current_A, stage, broad_requests, broad_requests2, mpi_rank, mpi_procs, f_mapping, rows, type_recv1, type_recv2, triang_buff_recv, comm_shape, comm_mechanism, quotient, comm);
	}
}
