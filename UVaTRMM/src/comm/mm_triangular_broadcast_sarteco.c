/********************************************************************
 * @file mm_triangular_common										*
 * @addtogroup mm_triangular										*
 * @brief Common functions definitions for mm_triangular program.	*
 * @author Arturo González Escribano								*
 * @author Rocío Carratalá Sáez										*
 * @author Maria Inmaculada Santamaria Valenzuela					*
 * @author Yuri Torres de la Sierra									*
 * ******************************************************************/
#include "mm_triangular.h"
///@private Broadcast function for defined datatypes
void mpi_datatypeBcast_sarteco(
	double *buff, 		//Array 
	MPI_Datatype *type, //Type
	int proc, 			//Process sending/receiving
	Requests *requests,	//MPI_Request buffer
	MPI_Comm comm, 		//Communicator
	size_t block_before,
	int mpi_rank
){
	mpi_broadcast_splitted(buff+block_before, 1, proc, comm, requests, mpi_rank,  *type);
	/*MPI_CHECK(
		MPI_Ibcast(
			buff+block_before, 
			1, 
			*type, 
			proc, 
			comm, 
			&requests -> arr[proc]
		)
	);*/
	if(*type != MPI_DOUBLE){
		MPI_Type_free( type );
	}
}

///@private Broadcast function for defined datatypes
void mpi_datatypeBcast_sarteco_quotient(
	double *buff, 		//Array 
	MPI_Datatype *type, //Type
	int proc, 			//Process sending/receiving
	Requests *requests,	//MPI_Request buffer
	MPI_Comm comm, 		//Communicator
	size_t block_before,
	int mpi_rank, 
	unsigned int quotient
){
	mpi_broadcast_splitted_quotient(
		buff+block_before, 1, proc, comm, requests, mpi_rank,  *type, quotient
	);
	/*MPI_CHECK(
		MPI_Ibcast(
			buff+block_before, 
			1, 
			*type, 
			proc, 
			comm, 
			&requests -> arr[proc]
		)
	);*/
		if(*type != MPI_DOUBLE){
			MPI_Type_free( type);
		}
}


///@private Broadcast function for defined datatypes
void mpi_datatypeBcast_sarteco_quotientArray(
	double *buff, 		//Array 
	MPI_Datatypes *type, //Type
	int proc, 			//Process sending/receiving
	Requests *requests,	//MPI_Request buffer
	MPI_Comm comm, 		//Communicator
	size_t block_before,
	int mpi_rank
){
	#if DEBUG
	printf("[%d | %d] Just before broad splitted quotient array Block Before %ld\n", mpi_rank, proc, block_before);
	fflush(stdout);
	#endif 
	mpi_broadcast_splitted_quotientArray(
		buff+block_before, proc, comm, requests, mpi_rank,  *type
	);
	/*for (int i = 0; i < type->total;i++){
		if(type->arr[i] != MPI_DOUBLE && type->sizes[i] > 0){
			MPI_Type_free( &type -> arr[i]);
		}
	}
	free(type->arr);
	free(type->sizes);
	type->total=0;
	type->arr = NULL;*/
	#if DEBUG
	printf("[%d | %d] Just after broad splitted quotient array Block Before %ld --> \n", mpi_rank, proc, block_before);
	fflush(stdout);
	#endif 
}

///@private Sets up mat_A_current
void matrix_rangeSet_mat_A_current(
	MatrixRange *mat_A_current,
	MatrixRange mat_A_local,
	MatrixRange mat_A_remote,
	int 		stage,
	int 		mpi_procs,
	int 		mpi_rank,
	Map_function f_mapping,
	Ranges *mat_A_current_rows
){
	unsigned long int block_before_A = 0;
	//Set mat_A_current
    //mat_A_current_rows -> array[stage] = mat_A_local.rows;
	if (mpi_rank == stage){
        matrix_rangePoint_to(mat_A_current, mat_A_local);
    } else {
        mat_A_current_rows -> array[stage] = f_mapping(mat_A_remote.rows.size, mpi_procs, stage, mpi_rank);
		block_before_A = mat_A_current_rows -> array[stage].begin * mat_A_remote.columns.size;
        mat_A_current -> mat 	= mat_A_remote.mat + block_before_A;
		mat_A_current -> rows 	= mat_A_current_rows -> array[stage];
		mat_A_current -> columns=mat_A_remote.columns;
    } 
}



///@private Sending phase of Broadcast Sarteco option
void comm_func_bcast_sarteco_send(
	MatrixRange mat_A, 	///Local matrix
	//MatrixRange *current_A,
	int stage, 			///Current stage
	Requests *requests,	///Request arrays
	Requests *requests2,
	int mpi_rank,		///Current proc
	Map_function f_mapping,	
	int rows,
	MPI_Datatype *type_send1,
	MPI_Datatype *type_send2,
	Buff *triang_buff_send,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism){
	#if DEBUG
	printf("[%d | %d] --> bcast sarteco send\n", mpi_rank, stage);
    fflush(stdout);
	printf_matrix_ranges_rank_stage(
			mat_A.mat,
			mat_A.rows.begin,
			mat_A.rows.size+mat_A.rows.begin,
			mat_A.columns.begin,
			mat_A.columns.begin+mat_A.columns.size,
			"ALOCALSENT", mpi_rank, stage, 0, 0
	);
	printf("[%d] ", mpi_rank);
	rangePrintf("VAMOS QUE RANGESENT", mat_A.rows);
	fflush(stdout);
	#endif 
	switch (comm_shape)
	{
	case SHAPE_FULL: 
		mpi_broadcast_splitted(mat_A.mat, mat_A.columns.size*mat_A.rows.size, stage, MPI_COMM_WORLD, requests, mpi_rank, MPI_DOUBLE);
		break;
	case SHAPE_TRIANG_2:
		mpi_datatypeBcast_sarteco(mat_A.mat, type_send1, stage, requests, MPI_COMM_WORLD, 0, mpi_rank);
		if ( comm_mechanism == MECH_BUFFER ) {
			buffIbcast(triang_buff_send, stage, requests2, MPI_COMM_WORLD, mpi_rank);
		} else {
			mpi_datatypeBcast_sarteco(mat_A.mat, type_send1, stage, requests2, MPI_COMM_WORLD, 0,  mpi_rank);
		}
		break;
	default:
		mpi_datatypeBcast_sarteco(mat_A.mat, type_send1, stage, requests, MPI_COMM_WORLD, 0,  mpi_rank);
		break;
	}

	//current_A -> mat = mat_A.mat;//realloc(mat_A, local_rows.size*columns*sizeof(double));
	//current_A -> rows = mat_A.rows;
	//current_A -> columns = mat_A.columns;
	#if DEBUG
	printf("[%d | %d] bcast sarteco send -->\n", mpi_rank, stage);
	fflush(stdout);
	#endif 
}

///@private receiving phase of broadcast sarteco option
void comm_func_bcast_sarteco_recv(
	MatrixRange *remote_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *requests,
	Requests *requests2,
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
	#if DEBUG
	printf("[%d | %d] --> bcast sarteco recv\n", mpi_rank, stage);
    fflush(stdout);
	#endif 
	current_A -> rows = f_mapping( rows, mpi_procs, stage, mpi_rank);
	//#if DEBUG
	//printf("[%d]", mpi_rank);
	//rangePrintf("VAMOS QUE RANGERECVD", current_A -> rows);
	//fflush(stdout);
	//printf_matrix_ranges_rank_stage(remote_A -> mat, remote_A -> rows.begin, remote_A -> rows.begin+remote_A -> rows.size, remote_A -> columns.begin, remote_A -> columns.begin+remote_A -> columns.size, "ARECVD3",mpi_rank, stage, 0, 0);
	//fflush(stdout);
	//#endif 
	set_type_recv(mpi_rank, current_A -> rows, remote_A -> columns.size, type_recv1, type_recv2, triang_buff_recv, comm_shape,comm_mechanism);
	size_t block_before = current_A->rows.begin*remote_A -> columns.size;
	#if DEBUG
	printf("[%d] Block before = %ld\n", mpi_rank, block_before);
	#endif 
	
	//MPI_CHECK( MPI_Waitall(requests -> size, requests -> arr, MPI_STATUSES_IGNORE));//Le peta la cabeza si hay dos envíos a la vez sobre el mismo puntero
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (stage > 0){ //--Revisar que esto esté en le mismo punto que en el de Arturo
		MPI_CHECK( MPI_Wait(&requests -> arr[stage-1], MPI_STATUS_IGNORE ) );
	}
	switch (comm_shape){
		case SHAPE_FULL:
			/*MPI_CHECK(
			MPI_Ibcast(
				remote_A -> mat + block_before, 
				current_A -> rows.size*remote_A->columns.size,
				MPI_DOUBLE, stage, MPI_COMM_WORLD, &requests -> arr[stage] 
			)
			);*/
			mpi_broadcast_splitted(
				remote_A -> mat, remote_A -> columns.size*current_A -> rows.size, 
				stage, MPI_COMM_WORLD, requests, 
				mpi_rank, MPI_DOUBLE
			);
			break;
		case SHAPE_TRIANG_2:
			mpi_datatypeBcast_sarteco(remote_A -> mat, type_recv1, stage, requests, MPI_COMM_WORLD, block_before,  mpi_rank);
			if ( comm_mechanism == MECH_BUFFER ) {
				buffIbcast(triang_buff_recv, stage, requests2, MPI_COMM_WORLD, mpi_rank);
				MPI_CHECK( MPI_Wait(&requests2 -> arr[stage], MPI_STATUS_IGNORE ) );
				buff__matrix_range(*triang_buff_recv, remote_A);
			} else {
				mpi_datatypeBcast_sarteco(remote_A -> mat, type_recv2, stage, requests2, MPI_COMM_WORLD, block_before, mpi_rank);
			}	
		break;
		default:
			mpi_datatypeBcast_sarteco(remote_A -> mat, type_recv1, stage, requests, MPI_COMM_WORLD, remote_A -> rows.begin*remote_A -> columns.size, mpi_rank);
		break;
	}
	
	#if DEBUG
	printf("[%d | %d] bcast sarteco recv -->\n", mpi_rank, stage);
    fflush(stdout);
	#endif 
}


/***
 * @private broadcast function for sarteco option 
 * Depending on stage & requests calls receive or send function.
 * */
void comm_func_bcast_sarteco(
	MatrixRange mat_A_local,
	MatrixRange *mat_A_global,
	MatrixRange *mat_A_current,
	int stage,
	Requests *requests,
	Requests *requests2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	Type_function create_type,
	int rows,
	MPI_Datatype *type1,
	MPI_Datatype *type2,
	Buff *triang_buff,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism
) {
	unsigned long int block_before = 0;
	
	if (mpi_rank == stage) {
		mat_A_current -> mat = mat_A_local.mat;
		mat_A_current -> columns = mat_A_local.columns;
	} else {
		mat_A_current -> columns = mat_A_global -> columns;
		block_before = mat_A_current->rows.begin*mat_A_global -> columns.size;
		mat_A_current -> mat = mat_A_global -> mat + block_before;
	}
	#if DEBUG
	if (mpi_rank==stage) matrix_rangePrintf_rank_stage(*mat_A_current, "mat_A_currentSARTECO", mpi_rank, stage, 0);
	#endif 
	if (comm_shape != SHAPE_FULL){
		//matrix_rangePrintf_rank_stage(*mat_A_current, "matAC", mpi_rank, -1, 0);
		*type1 = create_type( 
			mat_A_current -> rows, 
			mat_A_current -> columns.size 
		);
		MPI_CHECK( MPI_Type_commit( type1 ) );
		/*MPI_CHECK( 
			MPI_Ibcast( 
				mat_A_current->mat, 1, 
				*type1, stage, 
				MPI_COMM_WORLD, 
				&requests -> arr[ stage ] 
			) 
		);*/
		mpi_broadcast_splitted(mat_A_current->mat, 1, stage,MPI_COMM_WORLD, requests, mpi_rank, *type1);
		MPI_Type_free(type1);
		if (comm_shape == SHAPE_TRIANG_2){
			if ( comm_mechanism == MECH_TYPE ) {
				*type2 = create_type_triangle(mat_A_current -> rows, mat_A_current -> columns.size);
				//MPI_CHECK( MPI_Ibcast( mat_A_current->mat, 1, *type2, stage, MPI_COMM_WORLD, &requests -> arr[ stage ] ) );
				mpi_datatypeBcast_sarteco( mat_A_current->mat, type2, stage, requests, MPI_COMM_WORLD, 0, mpi_rank );
				MPI_Type_free(type2);
			}	else { //Buffer
				*triang_buff = set_triangular_buffer_data(mat_A_current-> rows, mat_A_current->columns.size, mat_A_current->mat, mpi_rank);
				//MPI_CHECK(MPI_Ibcast(triang_buff -> buffer, triang_buff -> size, MPI_DOUBLE, stage, MPI_COMM_WORLD, &requests2 -> arr[stage]));
				buffIbcast(triang_buff,stage, requests, MPI_COMM_WORLD, mpi_rank);
			}
		}
	} else {
		//matrix_rangePrintf_rank_stage(*mat_A_current, "CURRENT_TEST", mpi_rank, stage, 0 );
		//printf("%d\n", requests-> arr[stage]);
		fflush(stdout);
		#if DEBUG
		printf("[%d | %d] SENDING SIZE %ld\n", mpi_rank, stage, mat_A_current -> rows.size*mat_A_current -> columns.size);
		#endif 
		/*MPI_CHECK( 
			MPI_Ibcast( 
				mat_A_current->mat,
				//mat_A_current -> rows.size*mat_A_current -> rows.size, 
				mat_A_current -> rows.size*mat_A_current -> columns.size, 
				MPI_DOUBLE, stage, 
				MPI_COMM_WORLD, 
				&requests -> arr[ stage ] 
			) 
		);*/
		mpi_broadcast_splitted(mat_A_current->mat, mat_A_current->rows.size*mat_A_current->columns.size,stage, MPI_COMM_WORLD, requests, mpi_rank, MPI_DOUBLE);
	}
}
