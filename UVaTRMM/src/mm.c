/****************************************************************************************
 * @file mm_triangular.c																*
 * @addtogroup mm_triangular															*
 * @brief Multiplication of distributed square matrices, one of them lower triangular	*
 * @author Arturo González Escribano													*
 * @author Rocío Carratalá Sáez															*
 * @author Maria Inmaculada Santamaria Valenzuela										*
 * @author Yuri Torres de la Sierra														*
 * Goal: 																				*
 * * Given A & B squared-matrices with A lower triangular or rectangular returns A*B	*
 * Compilation: 																		*
 * * 1) Go to ./triangular product														*
 * * 2) Initialize environment by ". ./setup.sh"										*
 * Execution parameters: 																*
 * * 1) Matrix size																		*
 * * 2) Values generation param															* 
 * * * a_id : generates A = Id, B = consecutive values									*
 * * * b_id : generates A = consecutive values, B = Id									*
 * * * values: generates A & B consecutive values										*
 * * Execution line example:															*
 * 	rm output*;srun -n 2 -w manticore valgrind --log-file="valgrind.txt" 				*
 * 	--track-origins=yes ./mm_triangular.x 10  a_id_check output reg owner_broad full  	*
 * 	> broadcast_test.txt																*
 * Test file: ./tests.sh																*
 ****************************************************************************************/
//#ifndef GLOBAL_VARIABLES
//#define GLOBAL_VARIABLES
//#endif 
#include "mm_triangular.h"

///@private Main program
int main( int argc, char *argv[] ) {
	
	char *file_path="./tests/error/";//./tests/outputs/\0";
	char *error_filename;
	
	//////////////////
	/// VARIABLES ///
	////////////////	
	//MPI
	int mpi_rank, mpi_procs;
	unsigned int block_limit_denominator=0;	
	//Matrices
	//-- Global matrices
	MatrixRange mat_A_global	= matrix_range(RECTANGULAR);
	//-- Parts mpi_rank owns
	MatrixRange mat_A_local		= matrix_range(RECTANGULAR);
	//-- Stage Matrices
	MatrixRanges mat_A_stages = {0, NULL};
	//MatrixRange mat_A_remote	= matrix_range(RECTANGULAR); 
	MatrixRange mat_B_local 	= matrix_range(RECTANGULAR);
	MatrixRange mat_C_local 	= matrix_range(RECTANGULAR);
	size_t block_before_C;
	
	//Program parameters
	int do_check = 0; 
	int check_result = 1;
	int check_result_reduced = 1;
	Params_init init_scheme;
	Params_output output;
	Map_function f_mapping = map_regular;
	//Type_function create_type;
	Params_comm_scheme comm_scheme;
	Params_comm_func comm_func = FUNC_BROADCAST;
	//Params_comm_shape comm_shape = SHAPE_FULL;
	//Params_comm_mechanism comm_mechanism = MECH_BUFFER;
	//Auxiliar vars
	int stage; //For operation stages
	//Time vars for gettimeofday
	#if STATS
	struct timeval total_time, precomp_time, comm1_time, comm2_time, mult_time, for_time;
	double total_t 		= 0.0;
	double precomp_t 	= 0.0;
	double comm1_t 		= 0.0;
	double comm2_t 		= 0.0;
	double mult_t 		= 0.0;
	double for_t 		= 0.0;
	_Bool check_comm1	= 0;
	#endif 

	//Communication control
	//-- Main type Ibcast wait request
	Requestss requestss = {NULL, 0}; 	
	Requestss requestss2 = {NULL, 0};
	#if STATUSSES 	
	Statusess statusess = {NULL, 0};
	Statusess statusess2 = {NULL, 0};
	#endif 
	//-- Communicators indexes
	int remote_origin; 
	int comm_from, comm_to;
	
	///////////////////////////
	/// STARTING PROGRAM 	///
	///	MPI INITIALIZATION 	///
	///////////////////////////
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );	
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );//Total num of procs
	
	block_limit_denominator=0;

	//Building error file_name
	error_filename = (char *) calloc_abort(
		ERROR_CHAR_SIZE, 
		sizeof(char), 
		"error_filename", 
		MPI_COMM_WORLD, mpi_rank
	);

	do {
	} while(
		sprintf( 
		error_filename, 
		"%s%s/%d_%s_%s_%s_%s_%s_%s_%d.dat", 
		file_path, "full", mpi_procs, 
		argv[1], argv[2], argv[3], argv[4], 
		argv[6], argv[5], mpi_rank
		) <= 0
	);
	#if DEBUG
	printf("Error file name: %s\n", error_filename); 
	fflush(stdout);
	#endif 
	write_error_filename(error_filename, mpi_rank);
	fflush(stdout);
	
	///////////////////////
	/// READ PARAMETERS ///
	///////////////////////
	//-- Check total input parameters
	if ( argc < 7 ) show_usage( argv, mpi_rank, "Number of arguments" );
	//-- Set up sizes
	mat_A_global.rows.size = atoi( argv[1] );
	mat_A_global.columns.size = mat_A_global.rows.size;

	#if DEBUG
	printf("[%d] Aglobal Rows %ld\n", mpi_rank, mat_A_global.rows.size);
	#endif 
	//-- Values generation & check options
	set_initial_scheme(argc, argv, mpi_rank, &init_scheme, &do_check);
	//-- Output option
	set_output_result(argc, argv, mpi_rank, &output);
	//-- Set communication scheme
	set_communication_scheme(argc, argv, mpi_rank, &comm_scheme, &comm_func);
	//-- Sanity check: Avoid idle processes to simplify code
	if ( mat_A_global.rows.size < mpi_procs ) show_usage( argv, mpi_rank, "Number of rows should be equal or greater than number of processors" );
	
	/////////////////////////////
	///	Partition and mapping ///
	/////////////////////////////	
	//-- Set local matrix rows
	mat_A_local.rows = f_mapping( mat_A_global.rows.size, mpi_procs, mpi_rank, mpi_rank );
	//-- Set local matrix columns
	mat_B_local.columns = map_regular( mat_A_global.columns.size, mpi_procs, mpi_rank, mpi_rank);
	
	////////////////////////////
	///	Initialize matrices ///
	//////////////////////////	
	//-- Set sizes
	mat_A_local.columns.size 	= mat_A_global.columns.size;
	mat_A_local.columns.begin	= 0;
	mat_B_local.rows  			= mat_A_global.rows;
	mat_C_local.rows 			= mat_A_global.rows;
	mat_C_local.columns 		= mat_B_local.columns;


	//-- Alloc array memmory
	matrix_rangeMat_calloc_abort(&mat_A_local, "mat_A", MPI_COMM_WORLD, mpi_rank);
	matrix_rangeMat_calloc_abort(&mat_B_local, "mat_B_local", MPI_COMM_WORLD, mpi_rank);
	matrix_rangeMat_calloc_abort(&mat_C_local, "mat_C_local", MPI_COMM_WORLD, mpi_rank);

	//-- Initialize A & B values
	matrix_rangeSetA(&mat_A_local, init_scheme, mat_A_global.columns.size);
	matrix_rangeSetB(&mat_B_local, init_scheme, mat_A_global.columns.size);
	
	#if DEBUG
	printf_matrix_ranges_rank_stage(mat_A_local.mat, mat_A_local.rows.begin, mat_A_local.rows.size+mat_A_local.rows.begin, mat_A_local.columns.begin, mat_A_local.columns.begin+mat_A_local.columns.size, "matA", mpi_rank,-1, 0,0);
	printf("[%d] ", mpi_rank);
	rangePrintf("local_cols", mat_B_local.columns);
	printf_matrix_ranges_rank_stage(mat_B_local.mat, mat_B_local.rows.begin, mat_B_local.rows.size+mat_B_local.rows.begin, mat_B_local.columns.begin, mat_B_local.columns.begin+mat_B_local.columns.size, "matB", mpi_rank,-1, 0,0);
	#endif 
	
	//////
	// Timestamp: 
	// start total & precompilation timming
	/////
	MPI_Barrier(MPI_COMM_WORLD);
	#if STATS
	//gettimeofday(&total_time, NULL);
	gettimeofday(&precomp_time, NULL);
	#endif 
	
	requestss.size 	= mpi_procs;
	requestss.arr = (Requests*) calloc_abort(requestss.size, sizeof(Requests), "requestss", MPI_COMM_WORLD, mpi_rank);
	requestss.arr[0].size = 1;
	#if STATUSSES
	statusess.size	= mpi_procs;
	statusess.arr = (Statuses*) calloc_abort(requestss.size, sizeof(Statuses), "statusess", MPI_COMM_WORLD, mpi_rank);
	#endif 
	for (unsigned int i = 0; i < requestss.size; i++){
		requestss.arr[i] = requestsCalloc_abort(requestss.arr[0].size, MPI_COMM_WORLD, mpi_rank);
		#if STATUSSES
		statusess.arr[i] = statusesCalloc_abort(requestss.arr[0].size, MPI_COMM_WORLD, mpi_rank);
		#endif 
	}
	if (comm_func != FUNC_BROADCAST && comm_func != FUNC_BROADCAST_SARTECO){
		requestss2.size = mpi_procs;
		requestss2.arr = (Requests*) calloc_abort(requestss.size, sizeof(Requests), "requestss", MPI_COMM_WORLD, mpi_rank);
		requestss2.arr[0].size = 1;
		#if STATUSSES
		statusess2.size = mpi_procs;
		statusess2.arr = (Statuses*) calloc_abort(requestss.size, sizeof(Statuses), "statusess", MPI_COMM_WORLD, mpi_rank);
		#endif 
		for (unsigned int i = 0; i < requestss2.size; i++){
			#if DEBUG
				printf("[%d | %d] INITIALIZING REQUESTS2ARR ", mpi_rank, i);
			#endif 
			requestss2.arr[i] = requestsCalloc_abort(requestss2.arr[0].size, MPI_COMM_WORLD, mpi_rank);
			#if STATUSSES
			statusess2.arr[i] = statusesCalloc_abort(requestss2.arr[0].size, MPI_COMM_WORLD, mpi_rank);
			#endif 
		}
	} 

	
	///////////////////////////////////////
	///	Initial setup 		   			///
	/// OWNER_BROAD_SARTECO: 			///
	/// Send A local and gets A global	///
	/// OTHER SCHEMES:					///
	/// Set current A to A local 		///
	///////////////////////////////////////
	
	if (comm_func == FUNC_BROADCAST_SARTECO && mpi_procs > 1){
		mat_A_stages = matrix_ranges(RECTANGULAR, (unsigned int) mpi_procs);
		//-- --> Communication 1 Phase 1: send-receive
		#if STATS
		gettimeofday(&comm1_time, NULL);
		check_comm1=1;
		#endif 
		for (stage = 0; stage < mpi_procs; stage++){
			if (mpi_rank == stage) {
				#if DEBUG
					printf("[%d | %d] Broadcast Send\n", mpi_rank, stage);
				#endif 
				//mat_A_stages.array[stage].mat = mat_A_local.mat;
				matrix_rangePoint_to(&mat_A_stages.array[stage], mat_A_local);
			} else {
				#if DEBUG
					printf("[%d | %d] Broadcast Receive\n", mpi_rank, stage);
				#endif 
				mat_A_stages.array[stage].rows = f_mapping( 
					mat_A_global.rows.size, mpi_procs, stage, mpi_rank 
				);
				mat_A_stages.array[stage].columns = mat_A_global.columns;
				matrix_rangeMat_calloc_abort(&mat_A_stages.array[stage], "matBstg", MPI_COMM_WORLD, mpi_rank);
			} 
			if (
					mat_A_stages.array[stage].columns.size > 0 
				&& 	mat_A_stages.array[stage].rows.size > 0
			){
				mpi_broadcast_splitted_quotient(
					mat_A_stages.array[stage].mat, 
					mat_A_stages.array[stage].columns.size*mat_A_stages.array[stage].rows.size, 
					stage, 
					MPI_COMM_WORLD, &requestss.arr[stage], mpi_rank,  MPI_DOUBLE, block_limit_denominator
				);
				#if STATUSSES
				statusess.arr[stage].size	= requestss.arr[stage].size;
				statusess.arr[stage].arr	= (MPI_Status *) realloc(statusess.arr[stage].arr, sizeof(MPI_Status)*2);
				#endif 
			}
			//matrix_rangeIbcast(&mat_A_stages.array[0], stage, &requests, MPI_COMM_WORLD,0, mpi_rank); 
		} //For send-receive
		#if STATS
		comm1_t += set_time(comm1_time);
		#endif 
	} else if (mpi_procs > 1) {
		//First current part is the original part of A owned by the proc.
		if (comm_func == FUNC_P2P && comm_scheme == COMM_OWNER){
			mat_A_stages = matrix_ranges(RECTANGULAR, 2);
			matrix_rangePoint_to(&mat_A_stages.array[1], mat_A_local);
		} else {
			mat_A_stages = matrix_ranges(RECTANGULAR, 1);
		}
		matrix_rangePoint_to(&mat_A_stages.array[0], mat_A_local);
	} else {
		if (comm_func == FUNC_P2P && comm_scheme == COMM_OWNER){
			mat_A_stages = matrix_ranges(RECTANGULAR, 2);
			matrix_rangePoint_to(&mat_A_stages.array[0], mat_A_local);
			matrix_rangePoint_to(&mat_A_stages.array[1], mat_A_local);
		} else {
			mat_A_stages = matrix_ranges(RECTANGULAR, 1);
			matrix_rangePoint_to(&mat_A_stages.array[0], mat_A_local);
		}
	}
	#if STATS
	precomp_t = set_time(precomp_time);
	#endif
	#if DEBUG
	printf("[%d] Before for stage\n", mpi_rank);
	fflush(stdout);
	#endif 

	///////////////////////////////////////
	///	PRODUCT LOOP 		   			///
	/// Multiplication in num procs. 	///
	/// stages in [0, .., num_procs-1] 	///
	///////////////////////////////////////
	#if STATS
	MPI_Barrier(MPI_COMM_WORLD);
	gettimeofday(&total_time, NULL);
    gettimeofday(&for_time, NULL);
	#endif 
	for ( stage = 0; stage < mpi_procs; stage++ ) {
		//Communication 1 (Phase 2 for bcast sarteco: adjust current)
		#if STATS
			gettimeofday(&comm1_time, NULL);
			check_comm1=1;
		#endif 
		if (comm_func == FUNC_BROADCAST_SARTECO && mpi_procs > 1){
			//-- Set mat_A_stages.array[0]
			#if DEBUG
			printf("[%d | %d] Before waitall\n", mpi_rank, stage);
			fflush(stdout);
			#endif 
			if (mpi_procs > 1){
				MPI_CHECK( MPI_Waitall(requestss.arr[stage].size, requestss.arr[stage].arr, MPI_STATUSES_IGNORE)) ;
			}
			//#if DEBUG
			printf("[%d | %d] After waitall\n", mpi_rank, stage);
			fflush(stdout);
			//matrix_rangePrintf_rank_stage(mat_A_stages.array[stage], "AStageAfterWaitall", mpi_rank, stage, 0);
			//#endif 
		} else if (mpi_procs > 1) {
			//-- Reset mat_A_remote data & requests array
			mat_A_stages.array[0].mat = NULL;
			requestsMakeNull(&requestss.arr[stage]);
			if (comm_func == FUNC_BROADCAST){
				if (mpi_rank == stage) {
					#if DEBUG
						printf("[%d | %d] Broadcast Send\n", mpi_rank, stage);
					#endif 
					matrix_rangePoint_to(&mat_A_stages.array[0], mat_A_local);
				} else {
					#if DEBUG
						printf("[%d | %d] Broadcast Receive\n", mpi_rank, stage);
					#endif 
					mat_A_stages.array[0].rows = f_mapping( 
						mat_A_global.rows.size, mpi_procs, stage, mpi_rank 
					);
					mat_A_stages.array[0].columns = mat_A_global.columns;
					if (mat_A_stages.array[0].mat != NULL){
						free(mat_A_stages.array[0].mat);
						mat_A_stages.array[0].mat=NULL;
					}
					matrix_rangeMat_calloc_abort(&mat_A_stages.array[0], "matBstg", MPI_COMM_WORLD, mpi_rank);
				} 

				if (mpi_procs > 1){
					mpi_broadcast_splitted_quotient(
						mat_A_stages.array[0].mat, 
						mat_A_stages.array[0].rows.size*mat_A_stages.array[0].columns.size,
						stage, MPI_COMM_WORLD, &requestss.arr[stage], 
						mpi_rank, MPI_DOUBLE, block_limit_denominator
					);	
					MPI_CHECK( MPI_Waitall(requestss.arr[stage].size, requestss.arr[stage].arr, MPI_STATUSES_IGNORE)) ; //Abajo peta
				}
				/*mpi_broadcast_splitted(
					mat_A_stages.array[0].mat, 
					mat_A_stages.array[0].rows.size*mat_A_stages.array[0].columns.size,
					stage, MPI_COMM_WORLD, &requestss.arr[stage], 
					mpi_rank, MPI_DOUBLE
				);*/
				#if STATUSSES
				statusess.arr[stage].size	= requestss.arr[stage].size;
				statusess.arr[stage].arr	= (MPI_Status *) realloc(statusess.arr[stage].arr, sizeof(MPI_Status));
				#endif 

			} else {///P2P
				/// Start asynchronous communication ////
				// Remote origin:
				// 	P2P mode: Skew mode, displace ramk with stage
				//remote_origin = ( comm_func == FUNC_P2P ) ? mpi_rank - ( stage + 1 ) : stage;
				remote_origin = mpi_rank - ( stage + 1 );
				if ( remote_origin < 0 ) {
					remote_origin = mpi_procs + remote_origin;
				}
				mat_A_stages.array[0].rows 		= f_mapping( mat_A_global.rows.size, mpi_procs, remote_origin, mpi_rank );
				mat_A_stages.array[0].columns 	= mat_A_global.columns;
				switch (comm_scheme)
				{
				case COMM_OWNER://Skew
					comm_from = remote_origin;
					comm_to = mpi_rank + ( stage + 1 );
					if ( comm_to > mpi_procs-1 ) comm_to -= mpi_procs;
					break;
				default://Comm_previous = Pipeline
					comm_from = (mpi_rank == 0) ? mpi_procs-1 : mpi_rank-1;
					comm_to = (mpi_rank == mpi_procs-1) ? 0 : mpi_rank+1;
					break;
				}
				/* Scheme "owner" skips allocation and comm. in the last stage */
				if ( comm_scheme != COMM_OWNER || stage < mpi_procs - 1 ) {
				//if ( stage < mpi_procs - 1 ) {
					if (comm_scheme == COMM_OWNER && mpi_procs > 1) {//OWNER SKEW
						requestsMakeNull(&requestss.arr[stage]);
						requestsMakeNull(&requestss2.arr[stage]);
						mpi_isend_splitted(
							mat_A_local.mat, mat_A_local.columns.size*mat_A_local.rows.size, 
							comm_to, MPI_COMM_WORLD, &requestss.arr[stage], mpi_rank, MPI_DOUBLE, 0, 1, mpi_rank, 1
						);
						if (mat_A_stages.array[0].mat != NULL){
							free(mat_A_stages.array[0].mat);
							mat_A_stages.array[0].mat = NULL;
						}
						mat_A_stages.array[0].rows 		= f_mapping( mat_A_global.rows.size, mpi_procs, remote_origin, mpi_rank );
						mat_A_stages.array[0].columns 	= mat_A_global.columns;
						matrix_rangeMat_calloc_abort(&mat_A_stages.array[0], "AStg", MPI_COMM_WORLD, mpi_rank);
						#if DEBUG
						printf("[ %d | %d ]requestss2.size=%d\n", mpi_rank, stage, requestss2.arr[stage].size);
						fflush(stdout);
						#endif 
						mpi_irecv_splitted(
							mat_A_stages.array[0].mat, 
							mat_A_stages.array[0].columns.size*mat_A_stages.array[0].rows.size, 
							comm_from, MPI_COMM_WORLD, &requestss2.arr[stage], 
							mpi_rank, MPI_DOUBLE, 0, 1, remote_origin // remote_origin //¿Será que este remote_origin debiera ser comm_from?
						);
						#if STATUSSES
						statusess.arr[stage].size=requestss.arr[stage].size;
						statusess.arr[stage].arr=(MPI_Status*) realloc(statusess.arr[stage].arr, statusess.arr[stage].size*sizeof(MPI_Status));
						statusess2.arr[stage].size=requestss2.arr[stage].size;
						statusess2.arr[stage].arr=(MPI_Status*) realloc(statusess2.arr[stage].arr, statusess2.arr[stage].size*sizeof(MPI_Status));
						#endif
					} else if (comm_scheme == COMM_PREVIOUS && mpi_procs > 1){//Pipeline
							///Allocate for remote part
							matrix_rangeMat_calloc_abort(&mat_A_stages.array[0], "mat_A_remote", MPI_COMM_WORLD, mpi_rank);
							///Send Local
							#if DEBUG
							rangePrintf_rank_stage("send_rows", mat_A_local.rows, mpi_rank, stage);
							#endif 
							mpi_isend_splitted(
								mat_A_local.mat, mat_A_local.columns.size*mat_A_local.rows.size, 
								comm_to, MPI_COMM_WORLD, &requestss.arr[stage], mpi_rank, MPI_DOUBLE, 0, 1, mpi_rank, 1
							);
							requestss2.arr[stage] = requestsCalloc_abort(1, MPI_COMM_WORLD, mpi_rank);
							//Allocate From remote part
							//mat_A_stages.array[0].rows 		= f_mapping( mat_A_global.rows.size, mpi_procs, remote_origin, mpi_rank );
							//mat_A_stages.array[0].columns 	= mat_A_global.columns;
							//matrix_rangeMat_calloc_abort(&mat_A_stages.array[0], "AStg", MPI_COMM_WORLD, mpi_rank);
							//#if DEBUG
							//matrix_rangePrintf_rank_stage(mat_A_stages.array[0], "BEFORECV", mpi_rank, stage, 0);
							//#endif 
							#if DEBUG
							rangePrintf_rank_stage("send_rows", mat_A_stages.array[0].rows, mpi_rank, stage);
							#endif 
							mpi_irecv_splitted(
								mat_A_stages.array[0].mat, 
								mat_A_stages.array[0].columns.size*mat_A_stages.array[0].rows.size, 
								comm_from, MPI_COMM_WORLD, &requestss2.arr[stage], 
								mpi_rank, MPI_DOUBLE, 0, 1, comm_from
							);
							#if STATUSSES
							statusess.arr[stage].size=requestss.arr[stage].size;
							statusess.arr[stage].arr=(MPI_Status*) realloc(statusess.arr[stage].arr, statusess.arr[stage].size*sizeof(MPI_Status));
							statusess2.arr[stage].size=requestss2.arr[stage].size;
							statusess2.arr[stage].arr=(MPI_Status*) realloc(statusess2.arr[stage].arr, statusess2.arr[stage].size*sizeof(MPI_Status)); 
							#endif 
					}
				} 
			}
		}
		
        #if DEBUG
        printf("[%d | %d] Before local product\n", mpi_rank, stage);
		fflush(stdout);
        #endif 
		#if STATS
		comm1_t += set_time(comm1_time);
		#endif 
		#if (DEBUG && STATS)
		printf("[%d | %d] comm1_t = %f\n", mpi_rank, stage, comm1_t);
		fflush(stdout);
		#endif 
		#if STATS
	    gettimeofday(&mult_time, NULL);
		#endif 
		/// 5.2. Compute matrix multiplication with current parts 
		// Only multiply the lower triangle to improve performance
		// and to avoid undesired initializations with 0s in 
		// previous+boxes and previous+triang combinations
		//	
		// Inverted j,k loops
		// TODO: Tiling
		//int i, j, k, row;
		
		if (comm_func == FUNC_BROADCAST_SARTECO){//BROAD_SARTECO
			#if DEBUG
			matrix_rangePrintf_rank_stage(mat_A_stages.array[stage], "AcurrentJustBeforeProductSARTECO", mpi_rank, stage, 0);
			printf("Just before rectangular product owner BROAD SARTECO\n");
			fflush(stdout);
			#endif 
			block_before_C = mat_A_stages.array[stage].rows.begin*mat_C_local.columns.size;
			mm_rectangular_rectangular_product(
				mat_A_stages.array[stage], mat_B_local, &mat_C_local, block_before_C, stage, mpi_rank
			);
			free(mat_A_stages.array[stage].mat);
			mat_A_stages.array[stage].mat = NULL;
		} else if (comm_func == FUNC_BROADCAST){//BROAD
			#if DEBUG
			matrix_rangePrintf_rank_stage(mat_A_stages.array[0], "AcurrentJustBeforeProductBROADNOSARTECO", mpi_rank, stage, 0);
			printf("Just before rectangular product owner BROAD\n");
			fflush(stdout);
			#endif 
			block_before_C = mat_A_stages.array[0].rows.begin*mat_C_local.columns.size;
			mm_rectangular_rectangular_product(
				mat_A_stages.array[0], mat_B_local, &mat_C_local, block_before_C, stage, mpi_rank
			);
		} else if (comm_scheme == COMM_OWNER) {//OWNER_SKEW
			#if DEBUG
			printf("Just before rectangular product owner skew\n");
			fflush(stdout);
			#endif 
			block_before_C = mat_A_stages.array[1].rows.begin*mat_C_local.columns.size;
			//block_before_C = mat_A_local.rows.begin*mat_C_local.columns.size;
			mm_rectangular_rectangular_product(
				mat_A_stages.array[1], mat_B_local, &mat_C_local, block_before_C, stage, mpi_rank
				//&mat_A_local, &mat_B_local, &mat_C_local, block_before_C, stage, mpi_rank
			);
			#if DEBUG
			printf("Just after product\n");
			fflush(stdout);
			#endif 
		} else {//PIPELINE
			#if DEBUG
			matrix_rangePrintf_rank_stage(mat_A_local, "AcurrentJustBeforeProductPIPELINE", mpi_rank, stage, 0);
			printf("Just before rectangular product owner pipeline\n");
			fflush(stdout);
			#endif 
			block_before_C = mat_A_local.rows.begin*mat_C_local.columns.size;
			mm_rectangular_rectangular_product(
				mat_A_local, mat_B_local, &mat_C_local, block_before_C, stage, mpi_rank
			);
		}

		/*#if DEBUG
		printf("[%d | %d] Just after cblas\n", mpi_rank, stage);
		fflush(stdout);
		#endif
		#if DEBUG
		printf("[%d | %d] Just after product\n", mpi_rank, stage);
		printf_matrix_ranges_rank_stage(mat_C_local.mat, mat_A_stages.array[0].rows.begin, mat_A_stages.array[0].rows.begin+mat_A_stages.array[0].rows.size, mat_C_local.columns.begin, mat_C_local.columns.begin+mat_C_local.columns.size, "mat_C_current",mpi_rank, stage, block_before_C,0);
		fflush(stdout);
		#endif 
		#if DEBUG
		printf("[%d | %d] After local product\n", mpi_rank, stage);
		printf_matrix_ranges_rank_stage(mat_C_local.mat, mat_A_stages.array[0].rows.begin, mat_A_stages.array[0].rows.begin+mat_A_stages.array[0].rows.size, mat_C_local.columns.begin, mat_C_local.columns.begin+mat_C_local.columns.size, "mat_C_ALP",mpi_rank, stage, block_before_C,0);
		fflush(stdout);
        #endif */
		#if STATS
		mult_t += set_time(mult_time);
		#endif 
		/*if (init_scheme == INIT_VALUES && do_check){
			check_results_init_values(mat_B_local.columns,mat_A_stages.array[0].rows, mat_A_global.columns.size, mat_A_stages.array[0].mat, mat_B_local.mat, mat_C_local.mat, &check_result,mat_A_global.rows.size);
		}*/
		/*if (comm_scheme == COMM_PREVIOUS ){
			//Set up next sending matrix
			mat_A_local.columns = mat_A_stages.array[0].columns;
			mat_A_local.rows = mat_A_stages.array[0].rows;
			mat_A_local.mat = mat_A_stages.array[0].mat;
			#if DEBUG
			printf_matrix_ranges_rank_stage(mat_A_stages.array[0].mat, mat_A_stages.array[0].rows.begin, mat_A_stages.array[0].rows.size+mat_A_stages.array[0].rows.begin, 0, mat_A_stages.array[0].columns.size, "matACurrentchanged", mpi_rank,stage, 0);
			printf_matrix_ranges_rank_stage(mat_A_local.mat, mat_A_local.rows.begin, mat_A_local.rows.size+mat_A_local.rows.begin, 0, mat_A_global.columns.size, "matAChanged", mpi_rank,stage, 0);
			fflush(stdout);
			#endif 
		}*/
		#if STATS
		gettimeofday(&comm2_time, NULL);
		#endif
		/// 5.3. End async. comm. and rotate A structures for next iteration */
		///-- Scheme "owner" skips allocation and comm. in the last stage */
		/*#if DEBUG
		printf_matrix_ranges_rank_stage(mat_C_local.mat, mat_A_stages.array[0].rows.begin, mat_A_stages.array[0].rows.begin+mat_A_stages.array[0].rows.size, mat_B_local.columns.begin, mat_B_local.columns.size+mat_B_local.columns.begin, "c_check_outside_beforewaits", mpi_rank, stage, 0,0);
		fflush(stdout);
		#endif */
		
		if (
			(comm_func != FUNC_BROADCAST) &&
			! ( 
						comm_scheme == COMM_OWNER 
					&& 	stage == mpi_procs - 1 
				) 
		) {
			#if DEBUG
			printf("[%d | %d] --> Else if !commowner stage=mpi_procs-1\n", mpi_rank, stage);
			//printf_matrix_ranges_rank_stage(mat_C_local.mat, mat_A_stages.array[0].rows.begin, mat_A_stages.array[0].rows.begin+mat_A_stages.array[0].rows.size, mat_B_local.columns.begin, mat_B_local.columns.size+mat_B_local.columns.begin, "c_check_outside_bwaitall", mpi_rank, stage,0,0);
			fflush(stdout);
			#endif 
			/*if (
				(comm_scheme != COMM_PREVIOUS
				&& comm_func !=FUNC_BROADCAST_SARTECO) 
				|| (comm_func == FUNC_P2P && comm_scheme ==COMM_OWNER)//Quitado probando  errores con MPI_WAIT
			){*/
			if (
				comm_func !=FUNC_BROADCAST_SARTECO && comm_func != FUNC_BROADCAST //De OWNER_BROAD he tenido que pasar el wait arriba en triangular e iba. Testeando
			){
				if (mpi_procs > 1){
					#if DEBUG 
					printf("[%d | %d] Before waitall requestss\n", mpi_rank, stage);
					fflush(stdout);
					#endif 
					MPI_CHECK( 
						MPI_Waitall(requestss.arr[stage].size, requestss.arr[stage].arr, MPI_STATUSES_IGNORE) ;
					);
					#if DEBUG 
					printf("[%d | %d] Before waitall requestss2 %d \n", mpi_rank, stage, requestss2.arr[stage].size);
					fflush(stdout);
					#endif 
					MPI_CHECK( 
						MPI_Waitall(requestss2.arr[stage].size, requestss2.arr[stage].arr, MPI_STATUSES_IGNORE) ;
					);
					//#if DEBUG 
					printf("[%d | %d] After waitalls\n", mpi_rank, stage);
					fflush(stdout);
					//matrix_rangePrintf_rank_stage(mat_A_stages.array[0], "STAGEAFTERWAITALLS", mpi_rank, stage, 0);
					fflush(stdout);
					//#endif 
					if (comm_func == FUNC_P2P){
						if (comm_scheme == COMM_PREVIOUS){ //Pipeline
							matrix_rangePoint_to(&mat_A_local, mat_A_stages.array[0]);
						} else {//Owner skew
							matrix_rangePoint_to(&mat_A_stages.array[1], mat_A_stages.array[0]);
							//matrix_rangePoint_to(&mat_A_stages.array[1], mat_A_stages.array[0]);
						}
					}
				}
			} 
			
			/* Scheme "owner" do not deallocate mat_A in first stage */
			/*if ( 
				! ( 
							comm_scheme == COMM_OWNER 
						&& 	stage == 0 
					) 
				&& comm_func != FUNC_BROADCAST_SARTECO 
			) {
				printf("[%d | %d] He entrado aqui\n", mpi_rank, stage);
				free( mat_A_stages.array[0].mat ); -- > Da problemas en pipeline
			}*/
			
		}
		#if STATS
		comm2_t += set_time(comm2_time);
		#endif 
		//}
	}
	#if DEBUG
	printf("[%d] After for stage loop\n", mpi_rank);
	fflush(stdout);
	#endif 
	#if STATS
	for_t = set_time(for_time);
	gettimeofday(&comm2_time, NULL);
	#endif 
	// Scheme "previous", mat_A has been sent and deallocated in the first stage,
	// and reallocated and received in the last stage
	// The mat_A pointer should be restored for the free at the end of the program
	if (comm_scheme == COMM_PREVIOUS ){
		//Set up next sending matrix
		matrix_rangePoint_to(&mat_A_local, mat_A_stages.array[0]);
		//mat_A_local.columns = mat_A_stages.array[0].columns;
		//mat_A_local.rows = mat_A_stages.array[0].rows;
		//mat_A_local.mat = mat_A_stages.array[0].mat;
		#if DEBUG
		//printf_matrix_ranges_rank_stage(mat_A_stages.array[0].mat, mat_A_stages.array[0].rows.begin, mat_A_stages.array[0].rows.size+mat_A_stages.array[0].rows.begin, 0, mat_A_stages.array[0].columns.size, "matACurrentchanged", mpi_rank,stage, 0,0);
		printf_matrix_ranges_rank_stage(mat_A_local.mat, mat_A_local.rows.begin, mat_A_local.rows.size+mat_A_local.rows.begin, 0, mat_A_global.columns.size, "matAChanged", mpi_rank,stage, 0,0);
		//fflush(stdout);
		#endif 
	}
	
	// Annadido post-navidad 
	// Precomputed derived data types can be freed now
	/*
	Quitado por segfaults
	if ( comm_scheme == COMM_OWNER ) {
		
		if ( comm_shape == SHAPE_BOXES || comm_shape == SHAPE_TRIANG_1){
			MPI_Type_free( &type_send1 );
		} else if ( comm_shape == SHAPE_TRIANG_2) {
			MPI_Type_free( &type_send1 );
			MPI_Type_free( &type_send2 );
		}
		if ( comm_func == FUNC_BROADCAST_SARTECO || comm_func == FUNC_BROADCAST) {
			free( mat_A_remote.mat);
		}
	}
	*/

	#if DEBUG
	printf("[%d] Before check results\n", mpi_rank);
	fflush(stdout);
	#endif 
	#if STATS
	comm2_t += set_time(comm2_time);
	#endif 
	MPI_Barrier(MPI_COMM_WORLD);
	#if STATS
	total_t = set_time(total_time);
	#endif 
    /// 6. Check results 
	///-- Check for values is done up
	if ( do_check ) {
		#if DEBUG
		if ( mpi_rank == 0 ) {
			printf("[%d] Starting results check...\n", mpi_rank);
			fflush(stdout);
		}
		#endif 
		// Check that values on C are the same as the values on B
		if (init_scheme == INIT_A_ID){
			check_result = 1;
			#if DEBUG
			printf_matrix_ranges_rank_stage(mat_C_local.mat, 0,mat_A_global.rows.size,mat_B_local.columns.begin, mat_B_local.columns.size+mat_B_local.columns.begin, "c_before_check_", mpi_rank, stage,0,0);
			fflush(stdout);
			#endif 
			check_results_aid(mat_B_local.columns, mat_A_global.rows.size, mat_B_local.mat, mat_C_local.mat, &check_result);
		} else if (init_scheme == INIT_VALUES){//Por comprobar que este bien puesto aqui -> Usar test exterior
			check_results_init_values(
				mat_B_local.columns,
				mat_A_local.rows, 
				mat_A_local.columns.size, 
				mat_A_local.mat, 
				mat_B_local.mat,
				mat_C_local.mat, &check_result,
				mat_B_local.rows.size, 
				mpi_rank
				);
		}
		#if DEBUG
		printf("\t[%d] Reduce the result...\n", mpi_rank);
		fflush(stdout);
		#endif
		// Reduce the result, rank 0 prints the global result
		check_result_reduced = 1;
		#if DEBUG
		printf("[%d | %d] Check result %d\n", mpi_rank, stage, check_result);
		#endif
		MPI_CHECK( MPI_Reduce( &check_result, &check_result_reduced, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD ) );
		//MPI_Barrier(MPI_COMM_WORLD); --> Innecesario porque ya va en el reduce
		if ( mpi_rank == 0 ) {
			if ( check_result_reduced ) {
				printf("\t\t >> Check: OK\n");
				fflush(stdout);
			} else {
				printf("\t\t >> Check: FAILED !!\n");
				fflush(stdout);
			}
		} 
		#if DEBUG
		printf("\t[%d] Result reduced\n", mpi_rank);
		#endif
	}
	#if DEBUG
	printf("[%d] After check results\n", mpi_rank);
	fflush(stdout);
	#endif 
	/* 7. Output result */
	//#if DEBUG
	//matrix_rangePrintf_rank_stage(mat_C_local, "resC", mpi_rank, -1, 0);
	//#endif 
	if ( output == OUTPUT_YES ) {
		print_output_results_generic(
			mat_A_global.rows.size, 
			mat_A_global.columns.size, 
			mat_C_local.rows,
			mat_C_local.columns,
			mpi_rank, 
			mat_C_local.mat, 
			mpi_procs, "full"
		);
		fflush(stdout);
	}

	
	

	#if STATS
	/// 8. Write clocks 
	comm1_t = check_comm1 ? comm1_t : 0.0;
	set_ending_times(mpi_rank, total_t, precomp_t, comm1_t, comm2_t, mult_t, for_t);
	#endif
	#if NORM 
	matrix_rangeInfinityNorm(mat_C_local, mat_A_global, MPI_COMM_WORLD, mpi_rank);
	#endif
	/// 9. End 
	free(mat_A_global.mat);
	//free( mat_A_local.mat);
	free( mat_B_local.mat);
	free( mat_C_local.mat);
	#if DEBUG
	printf("BEFORE MPI_FINALIZE\n");
	fflush(stdout);
	#endif 
	MPI_Finalize();
	#if DEBUG
	printf("AFTER MPI_FINALIZE\n");
	fflush(stdout);
	#endif 
	return EXIT_SUCCESS;
}
