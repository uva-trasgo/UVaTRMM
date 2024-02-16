#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "mm_triangular.h"

///@private Main program
int main( int argc, char *argv[] ) {

	int exit_val = EXIT_SUCCESS;
	int mpi_rank, mpi_procs;
	unsigned int block_limit_denominator;

	//Matrices
	MatrixRange mat_A_global	= matrix_range(TRIANGULAR_D);
	//--Parts mpi_rank owns
	MatrixRange mat_A_local		= matrix_range(TRIANGULAR_D);
	//--Remote part (other proc owns) & current part (for operation)
	MatrixRanges mat_A_stages = {0, NULL};
	MatrixRange mat_B_local 	= matrix_range(RECTANGULAR);
	MatrixRange mat_C_local 	= matrix_range(RECTANGULAR);
	size_t block_before_C;	
	
	//Program parameters
	int do_check = 0; 
	int check_result = 1;
	int check_result_reduced = 1;
	Params_init init_scheme;
	Params_output output;
	Map_function f_mapping;
	Type_function create_type;
	Params_comm_scheme comm_scheme;
	Params_comm_func comm_func = FUNC_BROADCAST;
	Params_comm_shape comm_shape = SHAPE_FULL;
	Params_comm_mechanism comm_mechanism = MECH_BUFFER;
	int stage; //For operation stages
	
	//Time vars for gettimeofday
	struct timeval total_time, precomp_time, comm1_time, comm2_time, mult_time, for_time, aux_time;
	double total_t 		= 0.0;
	double mult_t 		= 0.0;
	_Bool check_comm1	= 0;
	
	//Communication buffers and data types
	//-- Type<...>1 is rectangle/combined
	//-- Type<...>2 is triangular part for non-buffered TRIANG_2 option
	MPI_Datatypes types_send1, types_send2, types_recv1, types_recv2;
	//-- Buffers for TRIANG2 buffered option
	Buff triang_buff_send = {NULL, 0};
	Buff triang_buff_recv = {NULL, 0};
	
	//Communication control
	//-- Main type Ibcast wait request
	Requestss requestss = {NULL, 0}; 	
	Requestss requestss2 = {NULL, 0}; 	
	//-- Communicators indexes
	int remote_origin; 
	int comm_from, comm_to;
	
	///////////////////////////
	/// STARTING PROGRAM 	///
	///MPI INITIALIZATION 	///
	///////////////////////////
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );	
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );//Total num of procs
	
	///////////////////////
	/// READ PARAMETERS ///
	///////////////////////
	if ( argc != 4 ) show_usage( argv, mpi_rank, "Number of arguments" );
	mat_A_global.rows.size = atoi( argv[1] );
	mat_A_global.columns.size = mat_A_global.rows.size;
	set_initial_scheme(argc, argv, mpi_rank, &init_scheme, &do_check);
	set_output_result(argc, argv, mpi_rank, &output);
	set_mapping_function(argc, argv, mpi_rank, &f_mapping);
	set_communication_scheme(argc, argv, mpi_rank, &comm_scheme, &comm_func);
	set_communication_shape(argc, argv, mpi_rank, &comm_shape, &comm_mechanism, &create_type);
	//-- Sanity check: Avoid idle processes to simplify code
	if ( mat_A_global.rows.size < mpi_procs ) show_usage( argv, mpi_rank, "Number of rows should be equal or greater than number of processors" );

	//////////////////////////////
	///	Partition and mapping  ///
	/// A -> local_rows x cols ///
	/// B -> rows x cols	   ///
	/// C -> rows x cols	   ///
	//////////////////////////////	
	//-- Set local matrix rows
	mat_A_local.rows = f_mapping( mat_A_global.rows.size, mpi_procs, mpi_rank , mpi_rank);
	//-- Set local matrix columns
	mat_B_local.columns = map_regular( mat_A_global.columns.size, mpi_procs, mpi_rank, mpi_rank);
	////////////////////////////
	///	Initialize matrices ///
	//////////////////////////	
	//-- Set sizes
	mat_A_local.columns.size 	= mat_A_global.columns.size;
	mat_A_local.columns.begin	= 0;
	mat_B_local.rows 			= mat_A_global.rows;
	mat_C_local.rows 			= mat_A_global.rows;
	mat_C_local.columns 		= mat_B_local.columns;
	//-- Alloc array memmory
	matrix_rangeMat_calloc_abort(&mat_A_local);
	matrix_rangeMat_calloc_abort(&mat_B_local);
	matrix_rangeMat_calloc_abort(&mat_C_local);
	//-- Initialize A & B values
	matrix_rangeSetA(&mat_A_local, init_scheme, mat_A_global.columns.size);
	matrix_rangeSetB(&mat_B_local, init_scheme, mat_A_global.columns.size);

    	block_limit_denominator=0;
	//////
	// Timestamp: 
	// start total & precompilation timming
	/////
	MPI_Barrier(MPI_COMM_WORLD);
	///////////////////////////////
	///	Send Derived datatypes ///
	/////////////////////////////
	//-- Derived data types for sending rectangles can be precalculated
	if (comm_scheme == COMM_OWNER && comm_func != FUNC_BROADCAST_SARTECO && comm_func != FUNC_BROADCAST){
		triang_buff_send = set_type_sendArray(mpi_rank, &types_send1, &types_send2, &triang_buff_send, &mat_A_local, comm_shape, comm_mechanism, comm_func, block_limit_denominator, MPI_COMM_WORLD, -1);
	}
	///////////////////////////////
	///	Send Derived datatypes ///
	/////////////////////////////
	requestss.size	= mpi_procs;
	requestss2.size	= mpi_procs;
	requestss.arr = (Requests*) calloc(requestss.size, sizeof(Requests));
	requestss.arr[0].size = 1;
	
	for (unsigned int i = 0; i < requestss.size; i++){
		requestss.arr[i] = requestsCalloc_abort(requestss.arr[0].size); //, MPI_COMM_WORLD, mpi_rank);
	}

	if (comm_shape == SHAPE_TRIANG_2 || comm_func == FUNC_P2P){
		requestss2.arr = (Requests*) calloc(requestss2.size, sizeof(Requests));
		for (unsigned int i = 0; i < requestss2.size; i++){
			requestss2.arr[0].size = 1;
			requestss2.arr[i] = requestsCalloc_abort(requestss2.arr[0].size); //, MPI_COMM_WORLD, mpi_rank);
		}
	}

	///////////////////////////////////////
	///	Initial setup 		   			///
	/// OWNER_BROAD_SARTECO: 			///
	/// Send A local and gets A global	///
	/// OTHER SCHEMES:					///
	/// Set current A to A local 		///
	///////////////////////////////////////
	#if STATS
	MPI_Barrier(MPI_COMM_WORLD);
	gettimeofday(&total_time, NULL);
	#endif 
	

	// RCSF if (comm_func == FUNC_BROADCAST_SARTECO && mpi_procs > 1){
	if (comm_func == FUNC_BROADCAST_SARTECO && mpi_procs >= 1){
		//-- --> Communication 1 Phase 1: send-receive
		mat_A_stages = matrix_ranges(TRIANGULAR_D, (unsigned int) mpi_procs);
		for (stage = 0; stage < mpi_procs; stage++){
			MPI_Barrier(MPI_COMM_WORLD);
			requestss.size=0;
			requestss2.size=0;
			if (mpi_rank == stage) {
				matrix_rangePoint_to(&mat_A_stages.array[stage], mat_A_local);
				triang_buff_send=set_type_sendArray(mpi_rank, &types_send1, &types_send2, &triang_buff_send, &mat_A_local, comm_shape, comm_mechanism, comm_func, block_limit_denominator, MPI_COMM_WORLD, stage);
			} else {
				mat_A_stages.array[stage].rows = f_mapping( 
					mat_A_global.rows.size, mpi_procs, stage, mpi_rank 
				);
				mat_A_stages.array[stage].columns = mat_A_global.columns;
				matrix_rangeMat_calloc_abort(&mat_A_stages.array[stage]);
				set_type_recvArray(
					mpi_rank, mat_A_stages.array[stage].rows,
					mat_A_global.columns.size, &types_send1, &types_send2, 
					&triang_buff_send, 
					comm_shape, comm_mechanism, block_limit_denominator, MPI_COMM_WORLD, stage
				);
			} 
			requestss.arr[stage].size = types_send1.total;
			size_t block_id, current_req=0;
        		for (block_id = 0; block_id < types_send1.total; block_id++){
                		if (types_send1.sizes[block_id] > 0){

					// v3
					int par = ((stage % 2) == 0) ? 1 : 0;
					int destiny;
					if (mpi_rank == stage && mpi_procs > 1 ) {
						destiny = mpi_rank+1;
						if (destiny >= mpi_procs) {
                                                        destiny = mpi_rank % 2 == 0 ? 1 : 0;
                                                }
						MPI_Send(
                                                                mat_A_stages.array[stage].mat,
                                                                1,
                                                                types_send1.arr[block_id],
                                                                destiny, //(mpi_rank+1) % mpi_procs,
                                                                0,
                                                                MPI_COMM_WORLD
                                                );
						destiny = mpi_rank + 2;
						if (destiny >= mpi_procs) {
							destiny = mpi_rank % 2 == 0 ? 0 : 1;
						}
						while (destiny != mpi_rank) {
							MPI_Isend(
                                                                mat_A_stages.array[stage].mat,
                                                                1,
                                                                types_send1.arr[block_id],
                                                                destiny, //(mpi_rank + isend) % mpi_procs,
                                                                0,
                                                                MPI_COMM_WORLD,
                                                                &requestss.arr[stage].arr[current_req]
                                                        );
							destiny += 2;
							if (destiny >= mpi_procs) {
                                                        	destiny = mpi_rank % 2 == 0 ? 0 : 1;
                                                	}
						}
					} else { if (mpi_procs > 1){
						int next_sender = stage + 1;
						if (next_sender >= mpi_procs) {
                                                	next_sender = stage % 2 == 0 ? 1 : 0;
                                                }
						if (mpi_rank == next_sender) {
							MPI_Recv(       
                                                                mat_A_stages.array[stage].mat,
                                                                1,
                                                                types_send1.arr[block_id],
                                                                stage,
                                                                0,
                                                                MPI_COMM_WORLD, MPI_STATUS_IGNORE
                                                	);
							destiny = mpi_rank + 2;
                                                	if (destiny >= mpi_procs) {
                                                        	destiny = mpi_rank % 2 == 0 ? 0 : 1;
                                                	}
                                                	while (destiny != mpi_rank) {
								MPI_Isend(
                                                                	mat_A_stages.array[stage].mat,
                                                        	        1,
                                                                	types_send1.arr[block_id],
                                                                	destiny, //(mpi_rank + isend) % mpi_procs,
                                                                	0,
                                                                	MPI_COMM_WORLD,
                                                               		&requestss.arr[stage].arr[current_req]
                                                        	);
								destiny += 2;
                                                        	if (destiny >= mpi_procs) {
                                                                	destiny = mpi_rank % 2 == 0 ? 0 : 1;
                                                        	}
							}
                                                } else {
							if ((mpi_rank % 2) == 0){
								if ( par ){
									MPI_Irecv(
                                                                		mat_A_stages.array[stage].mat,
                                                                		1,
                                                                		types_send1.arr[block_id],
                                                                		stage,
                                                                		0,
                                                                		MPI_COMM_WORLD,
                                                                		&requestss.arr[stage].arr[current_req]
                                                        		);
								} else {
									MPI_Irecv(
                                                                                mat_A_stages.array[stage].mat,
                                                                                1,
                                                                                types_send1.arr[block_id],
                                                                                next_sender,
                                                                                0,
                                                                                MPI_COMM_WORLD,
                                                                                &requestss.arr[stage].arr[current_req]
                                                                        );
								}
							} else {
								if ( par ){
                                                                        MPI_Irecv(
                                                                                mat_A_stages.array[stage].mat,
                                                                                1,
                                                                                types_send1.arr[block_id],
                                                                                next_sender,
                                                                                0,
                                                                                MPI_COMM_WORLD,
                                                                                &requestss.arr[stage].arr[current_req]
                                                                        );
                                                                } else {
                                                                        MPI_Irecv(
                                                                                mat_A_stages.array[stage].mat,
                                                                                1,
                                                                                types_send1.arr[block_id],
                                                                                stage,
                                                                                0,
                                                                                MPI_COMM_WORLD,
                                                                                &requestss.arr[stage].arr[current_req]
                                                                        );
                                                                }
							}
						}
					}}
                		} else {
                        		requestss.arr[stage].size = MAX(0,requestss.arr[stage].size-1);
                		}
   			}
		} 
	}
	
	///////////////////////////////////////
	///	PRODUCT LOOP 		   			///
	/// Multiplication in num procs. 	///
	/// stages in [0, .., num_procs-1] 	///
	///////////////////////////////////////
	// Aux variables
	long int index1, index2, c, f;
	size_t block_id, current_req, index, triangle_last_pos, mat_A_triangle_rows, mat_B_rectangle_columns, space_left_to_triangle, rows_begin_B, space_above_triangle;

	// Array to store the mult_t values
	struct timeval mult_times[mpi_procs*2];
	MPI_Barrier(MPI_COMM_WORLD);
	gettimeofday(&total_time, NULL); // Time start RCS
	for ( stage = 0; stage < mpi_procs; stage++ ) {
		//Communication 1 (Phase 2 for bcast sarteco: adjust current)
		//-- Set mat_A_current
		
		if (requestss.arr[stage].size > 0 && mpi_procs > 1){
			MPI_Waitall(requestss.arr[stage].size, requestss.arr[stage].arr, MPI_STATUSES_IGNORE);
		}
		
		#if STATS
	   	gettimeofday(&mult_times[stage * 2], NULL);
		#endif 

		/// 5.2. Compute matrix multiplication with current parts 
		block_before_C = mat_A_stages.array[stage].rows.begin * mat_C_local.columns.size;
		
        	double mat_B_copy[mat_B_local.columns.size * mat_A_stages.array[stage].rows.size];
        	index1, index2;
        	if (mpi_procs > 1){
                	for (f=0; f < mat_A_stages.array[stage].rows.size; f++){
				index1 = f * mat_B_local.columns.size;
                                index2 = (f + mat_A_stages.array[stage].rows.begin) * mat_B_local.columns.size;
				for (c=0; c < mat_B_local.columns.size; c++){
					mat_B_copy[index1 + c] = mat_B_local.mat[index2 + c];
                        	}
                	}
        	}

        	// START mm_triangular_rectangular_productRECTANGLE(mat_A_stages.array[stage], mat_B_local, mat_C_local, block_before_C, stage, mpi_rank);
        	cblas_dgemm (   
                	CblasRowMajor, CblasNoTrans, CblasNoTrans,
                	mat_A_stages.array[stage].rows.size, 
                	mat_B_local.columns.size, 
                	mat_A_stages.array[stage].rows.begin+1,
                	1.0, 
                	mat_A_stages.array[stage].mat,
                	mat_A_stages.array[stage].columns.size,
                	mat_B_local.mat,
                	mat_B_local.columns.size,
                	0.0, 
                	mat_C_local.mat+block_before_C,
                	mat_C_local.columns.size
        	);
		// END mm_triangular_rectangular_productRECTANGLE

		//START mm_triangular_rectangular_productTRIANGLE(mat_A_stages.array[stage], mat_B_local, mat_C_local, block_before_C, stage, mpi_rank);
		
		index;
        	triangle_last_pos = mat_A_stages.array[stage].rows.begin+mat_A_stages.array[stage].rows.size;
        	mat_A_triangle_rows     = MAX(0, mat_A_stages.array[stage].rows.size-1);
        	mat_A_triangle_rows     = MIN(mat_A_triangle_rows, mat_A_stages.array[stage].columns.size);
        	mat_A_triangle_rows     = MIN(mat_A_triangle_rows, mat_B_local.rows.size);
        	mat_B_rectangle_columns = MIN(mat_B_local.rows.begin+mat_B_local.rows.size, mat_B_local.columns.size);
        	space_left_to_triangle  = MAX(0,mat_A_stages.array[stage].rows.begin+mat_A_stages.array[stage].rows.size - mat_A_triangle_rows);
        	rows_begin_B            = space_left_to_triangle;
        	space_above_triangle    = MAX(0,mat_A_stages.array[stage].rows.size-mat_A_triangle_rows);
        	cblas_dtrmm (
                	CblasRowMajor, 
                	CblasLeft, 
                	CblasLower, 
                	CblasNoTrans, 
                	CblasNonUnit,
               		mat_A_triangle_rows,
                	mat_B_rectangle_columns, 
                	1.0,
                	mat_A_stages.array[stage].mat + space_above_triangle * mat_A_stages.array[stage].columns.size + space_left_to_triangle,
                	mat_A_stages.array[stage].columns.size,
                	mat_B_local.mat + mat_B_local.columns.size*rows_begin_B,
                	mat_B_local.columns.size
        	);
        	for ( f = rows_begin_B; f < rows_begin_B + mat_A_triangle_rows; f++ ){
                	index = f * mat_B_local.columns.size;
			for ( c = 0; c < mat_B_rectangle_columns; c++ ){
                        	mat_C_local.mat[index+c] += mat_B_local.mat[index+c];
                	}
        	}

        	if (mpi_procs > 1){
                	for (f=0; f < mat_A_stages.array[stage].rows.size; f++){
				index1 = f * mat_B_local.columns.size;
				index2 = (f + mat_A_stages.array[stage].rows.begin) * mat_B_local.columns.size;
                        	for (c=0; c < mat_B_local.columns.size; c++){
                                	mat_B_local.mat[index2 + c] = mat_B_copy[index1 + c];
                        	}
        		}
		}
		
		// END mm_triangular_rectangular_productTRIANGLE

		#if STATS
		gettimeofday(&mult_times[stage * 2 + 1], NULL);
		#endif 
	}
	MPI_Barrier(MPI_COMM_WORLD);
	gettimeofday(&aux_time, NULL);
        
	// Set final times (total and mult) of each MPI process
	total_t = (double) aux_time.tv_sec - total_time.tv_sec + (aux_time.tv_usec - total_time.tv_usec)/1000000.0;
        for ( stage = 0; stage < mpi_procs; stage++ ) {
		mult_t += (double) mult_times[stage * 2 + 1].tv_sec - mult_times[stage * 2].tv_sec + 
				(mult_times[stage * 2 + 1].tv_usec - mult_times[stage * 2].tv_usec)/1000000.0;
    	}

    	// 6. Print final execution times (total and product time)
	double total_t_red, mult_t_red;
	MPI_Reduce(&total_t, &total_t_red, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&mult_t, &mult_t_red, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	if ( mpi_rank == 0 ) {
    		printf( "\nTotal t: %lf\nTotal mpis: %d\nMult t: %lf\n", total_t_red, mpi_procs, mult_t_red);
	}

	// 7. Print norm of the result matrix and end execution
	matrix_rangeInfinityNorm(mat_C_local, mat_A_global, MPI_COMM_WORLD, mpi_rank);

	/// 9. End 
	free(mat_A_global.mat); free( mat_B_local.mat);	free( mat_C_local.mat);
	MPI_Finalize();
	return 0;
}
