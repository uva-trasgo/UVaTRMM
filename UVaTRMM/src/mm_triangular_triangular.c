/****************************************************************************************
 * @file mm_triangular_common															*
 * @addtogroup mm_triangular_triangular													*
 * @brief Main program of mm_triangular_triangular										*
 * @author Arturo González Escribano													*
 * @author Rocío Carratalá Sáez															*
 * @author Maria Inmaculada Santamaria Valenzuela										*
 * @author Yuri Torres de la Sierra														*
 * Goal: 																				*
 * * Given A & B squared-matrices with A & B lower triangular returns A*B				*
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
#define COMSPLITFUNC 1
#include "mm_triangular.h"
char *error_filename;
#define MULTIPLE_TRIANG_BUFF 0

/* Main program */
int main( int argc, char *argv[] ) {
	
    //signal(SIGTERM, manejador);
	char *file_path="./tests/error/";//./tests/outputs/\0";
	char *error_filename;
	
	//////////////////
	/// VARIABLES ///
	////////////////
	//MPI Variables
	//-- Current rank id
	int mpi_rank;
	//-- Total proccesors num
	int mpi_procs; 
	unsigned int block_limit_denominator;
	//-- Broadcast communicators array
	MPI_Comms comms;
	//Matrices
	//-- Global matrices
	MatrixRange mat_B_global 	= matrix_range(TRIANGULAR_D);
	//-- Local matrices
	MatrixRange mat_A_local 	= matrix_range(TRIANGULAR_D);
	MatrixRange mat_B_local 	= matrix_range(TRIANGULAR_D);
	MatrixRange mat_C_local 	= matrix_range(TRIANGULAR_D);
	#if BUFFERED_B
	MatrixRanges mat_B_stages; 
	#else 
	MatrixRange mat_B_stage	= matrix_range(TRIANGULAR_D);
	#endif 
	
	//Program parameters
	int do_check = 0;
	Params_init init_scheme;
	Params_output output;
	Map_function f_mapping;
	Type_function create_type;
	//Params_comm_scheme comm_scheme;
	//Params_comm_func comm_func = FUNC_BROADCAST;
	Params_comm_shape comm_shape;
	Params_comm_mechanism comm_mechanism = MECH_BUFFER;
	//Auxiliar vars
	//-- Operation stages (0 .. mpi_procs)
	unsigned int stage;
	size_t i, j;
	size_t index1;
	//size_t index2;
	int output_range[4];
	//-- Current communicator 
	
	_Bool Check_comm1 = 0;
	#if STATS
	//-- Time vars for gettimeoftoday
	struct timeval total_time, precomp_time, comm1_time, comm2_time, mult_time, for_time;//, aux_time;
	double total_t 		= 0.0;
	double precomp_t 	= 0.0;
	double comm1_t 		= 0.0;
	double comm2_t 		= 0.0;
	double mult_t 		= 0.0;
	double for_t 		= 0.0;
	#endif 
	//Communication buffers and data types
	//-- Type<...>1 is rectangle/combined
	//-- Type<...>2 is triangular part for non-buffered TRIANG_2 option
	//--- todo: Comprobar si con 1 solo basta
	MPI_Datatypes stage_type, stage_type2;
	//-- Buffer for triang2 buffered option
	#if MULTIPLE_TRIANG_BUFF
	Buff *triang_buffers;
	#else 
	Buff triang_buff = {NULL, 0};
	#endif 
	//Communication control
	//-- Wait requests
	Requestss requestss = {NULL, 0}; 	
	Requestss requestss2 = {NULL, 0}; 
	#if STATUSSES
	Statusess statusess = {NULL, 0};
	Statusess statusess2 = {NULL, 0};
	#endif 
	#if !BUFFERED_B
	size_t block_before_B = 0;
	#endif 

	///////////////////////////
	/// STARTING PROGRAM 	///
	///	MPI INITIALIZATION 	///
	///////////////////////////
	MPI_Init( &argc, &argv ); 
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
	//block_limit_denominator=8*mpi_procs;
	block_limit_denominator=0;
	#if BUFFERED_B
	mat_B_stages = matrix_ranges(TRIANGULAR_D, (unsigned int) mpi_procs);
	#endif 

	error_filename = (char *) calloc_abort(
		ERROR_CHAR_SIZE, 
		sizeof(char), 
		"error_filename", 
		MPI_COMM_WORLD, 
		mpi_rank
	);
	
	sprintf( 
		error_filename, 
		"%s%s/%d_%s_%s_%s_%s_%s_%s_%d.dat", 
		file_path, "triang_triang", mpi_procs, 
		argv[1], argv[2], argv[3], argv[4], 
		argv[6], argv[5], mpi_rank
	);
	#if DEBUG
	printf("Error file name: %s\n", error_filename); 
	fflush(stdout);
	#endif 
	write_error_filename(error_filename, mpi_rank);
	fflush(stdout);
	
	MPI_Barrier(MPI_COMM_WORLD);

	///////////////////////
	/// READ PARAMETERS ///
	///////////////////////
	if ( argc < 7 ) show_usage( argv, mpi_rank, "Number of arguments" );
	//-- Setup global ranges (rows and cols)
	mat_B_global.rows.size	 	= atoi( argv[1] );
	mat_B_global.columns.size	= mat_B_global.rows.size;
	mat_B_global.rows.begin		= 0;
	mat_B_global.columns.begin	= 0;
	//-- Setup communication & checking options
	set_initial_scheme(argc, argv, mpi_rank, &init_scheme, &do_check);
	set_output_result(argc, argv, mpi_rank, &output);
	set_mapping_function(argc, argv, mpi_rank, &f_mapping);
	/*if (mpi_procs >= floor(mat_B_global.rows.size/2))
	{
    	error_mpi_abort(MPI_COMM_WORLD, INVALID_PROCS_NUM, mpi_rank, "Too many procs for this size");
	}*/
	//comm_func=FUNC_BROADCAST_SARTECO;
	//comm_scheme=COMM_OWNER;
	set_communication_shape(argc, argv, mpi_rank, &comm_shape, &comm_mechanism, &create_type);
	//-- Sanity check: Avoid idle processes to simplify code */
	if ( mat_B_global.rows.size < mpi_procs ){
		show_usage( argv, mpi_rank, "Number of rows should be equal or greater than number of processors" );
	}
	///////////////////////////////
	///	Partition and mapping 	///
	/// A -> local_rows x cols ////
	/// B -> local_rows x cols ////
	/// C -> local_rows x cols ////
	///////////////////////////////
	// In this option all matrices
	// got same cols number (total)
	// and row number (mpi_rank's block)
	mat_A_local.rows 	= f_mapping( mat_B_global.rows.size, mpi_procs, mpi_rank, mpi_rank);
	mat_A_local.columns = mat_B_global.columns;
	
	mat_B_local.rows 	= mat_A_local.rows;
	//mat_B_local.columns = f_mapping(mat_B_global.columns.size, mpi_procs, mpi_rank);
	mat_B_local.columns = mat_B_global.columns;
	
	mat_C_local.rows 	= mat_A_local.rows;
	mat_C_local.columns = mat_B_local.columns;

	//mat_B_stage.rows=mat_B_local.rows;
	//mat_B_stage.columns=mat_B_local.columns;
	////////////////////////////
	///	Initialize matrices ///
	//////////////////////////
	//-- Allocation
	matrix_rangeMat_calloc_abort(&mat_A_local, 	"A_local", MPI_COMM_WORLD, mpi_rank);
	matrix_rangeMat_calloc_abort(&mat_B_local, 	"B_local", MPI_COMM_WORLD, mpi_rank);
	matrix_rangeMat_calloc_abort(&mat_C_local, 	"C_local", MPI_COMM_WORLD, mpi_rank);
	#if !BUFFERED_B
	matrix_rangeMat_calloc_abort(&mat_B_global, "Bglobal", MPI_COMM_WORLD, mpi_rank);
	#endif 
	//-- Initialization (A & B)
	matrix_rangeSetA(&mat_A_local, init_scheme, mat_B_global.columns.size);
	matrix_rangeSetB(&mat_B_local, init_scheme, mat_B_global.columns.size);
	#if DEBUG
	//matrix_rangePrintf_rank_stage(mat_A_local, "mat_A_local", mpi_rank, -1, 0);
	//matrix_rangePrintf_rank_stage(mat_B_local, "mat_B_local", mpi_rank, -1, 0);
	printf_matrix_ranges_rank_stage(mat_A_local.mat, mat_A_local.rows.begin, mat_A_local.rows.size+mat_A_local.rows.begin, mat_A_local.columns.begin, mat_A_local.columns.begin+mat_A_local.columns.size, "matA", mpi_rank,-1, 0,0);
	printf("[%d] ", mpi_rank);
	rangePrintf("local_rows", mat_B_local.rows);
	rangePrintf("local_cols", mat_B_local.columns);
	printf_matrix_ranges_rank_stage(mat_B_local.mat, mat_B_local.rows.begin, mat_B_local.rows.size+mat_B_local.rows.begin, mat_B_local.columns.begin, mat_B_local.columns.begin+mat_B_local.columns.size, "matB", mpi_rank,-1, 0,0);
	fflush(stdout);
	#endif 
	//////
	// Timestamp: 
	// start total & precompilation timming
	/////
	#if DEBUG
	printf("Before first barrier\n"); fflush(stdout);
	#endif 
	MPI_Barrier(MPI_COMM_WORLD);// <- se queja este barrier en dev
	#if DEBUG
	printf("After first barrier\n"); fflush(stdout);
	#endif 
	#if STATS
	//gettimeofday(&total_time, NULL);
	gettimeofday(&precomp_time, NULL);
	#endif 

	///////////////////////////////////////
	/// BUILD BROADCAST COMMUNICATORS 	///
	///////////////////////////////////////
	/// Build specific comms so			///
	/// mpi_rank sends to -> proc:		///
	/// 0 -> [1, ..., mpi_procs]		///
	/// 1 -> [2, ..., mpi_procs]		///
	///	...								///
	/// mpi_procs -> []					///
	///////////////////////////////////////	
	//-- Initialize to MPI_COM_WORLD
	#if DEBUG
	printf("Creating %d communicators\n", mpi_procs-1);
	fflush(stdout);
	#endif 
	#if COMMSPLITFUNC
	mm_triangular_triangularMpi_commSplit(&comms, MPI_COMM_WORLD, mpi_rank, mpi_procs);
	#else 
	    _Bool color;
    unsigned int comm;
    if (mpi_procs > 1){
		comms.buffer = calloc_abort(mpi_procs, sizeof(MPI_Comm ), "comms", MPI_COMM_WORLD, mpi_rank);
		for(comm = 0; comm < mpi_procs; comm ++){
			comms.buffer[comm] = MPI_COMM_WORLD;
		}
	    #if DEBUG
	    printf("--> Communicators split\n");
	    fflush(stdout);
	    #endif
	    for ( comm=1; comm<mpi_procs; comm++ ) { //Algo hay mal aqui
		//for ( comm=0; comm<mpi_procs; comm++ ) {
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
				    comms.buffer[comm], 	//Old communicator to split
				    color, 			//New communicator that the calling process is to be assigned to
				    mpi_rank, 		//Index order to use in new communicator
				    &comms.buffer[comm]   //New communicator handle
			    ) 
		    );//<<--Se queja aqui en dev
	    }
    } else {
		comms.buffer=NULL;
    }
	#if DEBUG
	int comm_size;
	for (int stage = 1; stage < mpi_procs; stage++){
		MPI_Comm_size( comms.buffer[stage], &comm_size );
		printf("[%d | %d] comm size %d\n", mpi_rank, stage, comm_size);
	}
	#endif 
	#endif 

	MPI_Barrier(MPI_COMM_WORLD); 
	//////////////////
	/// BROADCAST ////
	//////////////////
	//-- Allocate requests & statuses
	if (mpi_procs > 1){
	requestss.size	= mpi_procs;
	requestss2.size	= mpi_procs;
	requestss.arr = (Requests*) calloc_abort(requestss.size, sizeof(Requests), "requestss", MPI_COMM_WORLD, mpi_rank);
	#if STATUSSES
	statusess.size	= mpi_procs;
	statusess2.size	= mpi_procs;
	statusess.arr = (Statuses*) calloc_abort(requestss.size, sizeof(Statuses), "statusess", MPI_COMM_WORLD, mpi_rank);
	#endif 
	if (comm_shape==SHAPE_TRIANG_2){
		requestss2.arr = (Requests*) calloc_abort(requestss2.size, sizeof(Requests), "requestss2", MPI_COMM_WORLD, mpi_rank);
		#if STATUSSES
		statusess2.arr = (Statuses*) calloc_abort(requestss2.size, sizeof(Statuses), "statusess2", MPI_COMM_WORLD, mpi_rank);
		#endif 
	}
	requestss.arr[0].size = 1;
	for (unsigned int i=0; i < requestss.size; i++){
		requestss.arr[i] = requestsCalloc_abort(requestss.arr[0].size, MPI_COMM_WORLD, mpi_rank);
		#if STATUSSES
		statusess.arr[i] = statusesCalloc_abort(requestss.arr[0].size, MPI_COMM_WORLD, mpi_rank);
		#endif 
	}
	if (comm_shape == SHAPE_TRIANG_2){
		for (unsigned int i=0; i < requestss.size; i++){
			requestss2.arr[i]= requestsCalloc_abort(requestss.arr[0].size, MPI_COMM_WORLD, mpi_rank);
			#if STATUSSES
			statusess2.arr[i]= statusesCalloc_abort(requestss.arr[0].size, MPI_COMM_WORLD, mpi_rank);
			#endif 
		}
	}

	

	 
	for ( stage = 0; stage <= mpi_rank; stage++ ) {
		#if DEBUG
		printf("Bcast stage Loop\n");
		fflush(stdout);
		#endif 
		if (
				!(
						mpi_rank == stage 
					&& 	mpi_rank == mpi_procs-1
				)
		){
			//-- Set size &
			#if DEBUG
			printf("[%d | %d] Inside if\n", mpi_rank, stage);
			#endif 
			#if BUFFERED_B
			mat_B_stages.array[stage].rows = f_mapping( mat_B_global.rows.size, mpi_procs, stage, mpi_rank );
			mat_B_stages.array[stage].columns = mat_B_global.columns;
			#else 
			mat_B_stage.rows = f_mapping( mat_B_global.rows.size, mpi_procs, stage, mpi_rank );
			#endif 
			if (mpi_rank == stage){
				#if BUFFERED_B
				mat_B_stages.array[stage].mat = mat_B_local.mat;
				#else 
				mat_B_stage.mat = mat_B_local.mat;
				#endif
			} else {
				#if !BUFFERED_B
				block_before_B = mat_B_stage.rows.begin * mat_B_global.columns.size;
				mat_B_stage.mat = mat_B_global.mat+block_before_B;
				#endif 
			}
			//-- Send -> Receive
			//////
			/// Time: 
			/// Start | comm1_time
			/////
			#if STATS
			gettimeofday(&comm1_time, NULL);
			#endif
			Check_comm1=1;
			#if DEBUG
			printf("Send-receive\n");
			#endif 
			if (mpi_procs>1){
				#if DEBUG
				printf("[%d | %d] About to send/receive\n", mpi_rank, stage);
				fflush(stdout);
				#endif 
				if (mpi_rank == stage){
					#if BUFFERED_B
					matrix_rangePoint_to(&mat_B_stages.array[stage], mat_B_local);
					#else 
					matrix_rangePoint_to(&mat_B_stage, mat_B_local);
					#endif 
					triang_buff=set_type_sendArray(stage, &stage_type, &stage_type2, &triang_buff, &mat_B_local, comm_shape, comm_mechanism, FUNC_BROADCAST_SARTECO,  block_limit_denominator, MPI_COMM_WORLD, stage);
				} else {
					#if BUFFERED_B
					mat_B_stages.array[stage].columns = mat_B_local.columns;
					mat_B_stages.array[stage].rows=f_mapping( mat_B_global.rows.size, mpi_procs, stage, mpi_rank );
					matrix_rangeMat_calloc_abort(&mat_B_stages.array[stage], "matBstg", MPI_COMM_WORLD, mpi_rank);
					set_type_recvArray(mpi_rank, mat_B_stages.array[stage].rows, mat_B_global.columns.size, &stage_type, &stage_type2, &triang_buff, comm_shape, comm_mechanism, block_limit_denominator, MPI_COMM_WORLD, stage);
					#else 
					block_before_B 	= mat_B_stage.rows.begin * mat_B_global.columns.size;	
					mat_B_stage.mat = mat_B_global.mat + block_before_B;
					mat_B_stage.columns = mat_B_local.columns;
					mat_B_stage.rows=f_mapping( mat_B_global.rows.size, mpi_procs, stage, mpi_rank );
					set_type_recv(mpi_rank, mat_B_stage.rows, mat_B_global.columns.size, &stage_type, &stage_type2, &triang_buff, comm_shape, comm_mechanism);
					#endif
				}
				#if DEBUG
				printf("[%d | %d] Dtypes %d\n", mpi_rank, stage, stage_type.total);
				fflush(stdout);
				#if BUFFERED_B
				//matrix_rangePrintf_rank_stage(mat_B_stages.array[stage], "Bstg_presarteco", mpi_rank, stage,0);
				#else 
				//matrix_rangePrintf_rank_stage(mat_B_stage, "Bstg_presarteco", mpi_rank, stage,0);
				#endif 
				//matrix_rangePrintf_rank_stage(mat_B_local, "Bloc_presarteco", mpi_rank, stage,0);	
				#endif 
				if (comm_shape != SHAPE_FULL){
					#if BUFFERED_B
					#if DEBUG
					printf("[%d | %d] Bcast :)\n", mpi_rank, stage);
					#endif 
					//MPI_CHECK( MPI_Ibcast( mat_B_stages.array[stage].mat, 1, stage_type, 0, comms.buffer[stage], &requests.arr[stage]) );
					mpi_broadcast_splitted_quotientArray(
						mat_B_stages.array[stage].mat, 
						0, 
						comms.buffer[stage], &requestss.arr[stage],
						mpi_rank, 
						stage_type
					);
					//MPI_Type_free(&stage_type);
					//mpi_datatypesFree(&stage_type);
					if (comm_shape == SHAPE_TRIANG_2){
						if (comm_mechanism == MECH_TYPE){
							mpi_broadcast_splitted_quotientArray(mat_B_stages.array[stage].mat, 0, comms.buffer[stage], &requestss2.arr[stage], mpi_rank, stage_type2);
							//mpi_datatypesFree(&stage_type2);
						} else {//Buffer
							buffIbcast_quotient(&triang_buff, stage, &requestss2.arr[stage], MPI_COMM_WORLD, mpi_rank, block_limit_denominator);
						}
					}
					#if DEBUG
					printf("[%d | %d] Dtypes %d Requests %d After Send\n", mpi_rank, stage, stage_type.total, requestss.arr[stage].size);
					fflush(stdout);
					#endif 
				} else { //SHAPE_FULL
						#if DEBUG
						printf("[%d] send stg %d Ibcast\n", mpi_rank, stage);
						fflush(stdout);
						//matrix_rangePrintf_rank_stage(mat_B_stages.array[stage], "sentB", mpi_rank, stage, 0);
						fflush(stdout);
						#if STATUSSES
						for (int i = 0; i < requests.size; i++){
							printf("[%d | %d] req[%d] %d\n", mpi_rank, stage, i, requests.arr[i]);
						}
						#endif 
			
						#endif 
						int comm_size;
						if (comms.buffer == NULL){
							printf("WTF\n");
							MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
						}
						MPI_Comm_size( comms.buffer[stage], &comm_size );
						#if DEBUG
						printf("[%d | %d] currcomsize %d\n", mpi_rank, stage, comm_size);
						#endif 
						mpi_broadcast_splitted_quotient(
							mat_B_stages.array[stage].mat, 
							mat_B_stages.array[stage].rows.size*mat_B_stages.array[stage].columns.size,
							0,
							comms.buffer[stage], 
							&requestss.arr[stage],mpi_rank, MPI_DOUBLE, block_limit_denominator
						);
					} // If shape full / other
					#if (DEBUG && STATUSSES)
					statusess.arr[stage].size = requestss.arr[stage].size;
					statusess.arr[stage].arr = (MPI_Status*) realloc(statusess.arr[stage].arr, statusess.arr[stage].size*sizeof(MPI_Status));
					//if (comm_shape == SHAPE_TRIANG_2 || comm_func == FUNC_P2P){
					if (comm_shape == SHAPE_TRIANG_2){
						statusess2.arr[stage].size = requestss2.arr[stage].size;
						statusess2.arr[stage].arr = (MPI_Status*) realloc(statusess2.arr[stage].arr, statusess2.arr[stage].size*sizeof(MPI_Status));
					}
					#endif 
	
					#else
					MPI_CHECK( MPI_Ibcast( mat_B_stage.mat, 1, stage_type, 0, comms.buffer[stage], &requestss.arr[stage]););
					//MPI_Type_free(&stage_type);
					if (comm_shape == SHAPE_TRIANG_2){
						if (comm_mechanism == MECH_TYPE){
							MPI_CHECK( MPI_Ibcast( mat_B_stage.mat, 1, stage_type2, 0, comms.buffer[stage], &requests2.arr[ stage ] ););
							//MPI_Type_free(&stage_type2);
						} else {//Buffer
							MPI_CHECK(MPI_Ibcast(triang_buff.buffer, triang_buff.size, MPI_DOUBLE, 0, comms.buffer[stage], &requests2.arr[stage]););
						}
					}
					} else { //SHAPE_FULL
						#if DEBUG
						printf("[%d] send stg %d Ibcast\n", mpi_rank, stage);
						fflush(stdout);
						//matrix_rangePrintf_rank_stage(mat_B_stage, "sentB", mpi_rank, stage, 0);
						fflush(stdout);
						for (int i = 0; i < requests.size; i++){
							printf("[%d | %d] req[%d] %d\n", mpi_rank, stage, i, requests.arr[i]);
						}
						#endif 
						int comm_size;
						if (comms.buffer == NULL){
							printf("WTF\n");
							MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
						}
						MPI_Comm_size( comms.buffer[stage], &comm_size );
						#if DEBUG
						printf("[%d | %d] currcomsize %d\n", mpi_rank, stage, comm_size);
						matrix_rangeIbcast(&mat_B_stage, stage, &requests, comms.buffer[stage], 0, mpi_rank);
						#endif 
					} // If shape full / other
			#endif 
		}
		} //Avoid lat send
	} //For stage send - receive   
} else {
	matrix_rangePoint_to(&mat_B_stages.array[0], mat_B_local);
}
	#if DEBUG
	printf("[%d] After broad if\n", mpi_rank);
	fflush(stdout);
	#endif 
	//////
	/// Timestamp: 
	/// End		| precompilation
	/// Start 	|  product
	//////
	#if STATS
	precomp_t = set_time(precomp_time);
	MPI_Barrier(MPI_COMM_WORLD);
	gettimeofday(&total_time, NULL);
	gettimeofday(&for_time, NULL);
	#endif
	////////////////////
	/// PRODUCT LOOP ///
	////////////////////
	for ( stage = 0; stage <= mpi_rank; stage++ ) {
		/* Scheme "broadcast", wait for the end of the comm. corresponding to this stage */
		#if BUFFERED_B
		mat_B_stages.array[stage].rows = f_mapping( mat_B_global.rows.size, mpi_procs, stage, mpi_rank );
		#else 
		mat_B_stage.rows = f_mapping( mat_B_global.rows.size, mpi_procs, stage, mpi_rank );
		#endif
		if (
				!(
						mpi_rank == mpi_procs - 1 
					&& 	stage == mpi_procs-1
				) 
		){
			#if DEBUG
			printf("[%d | %d] Before Wait %d\n", mpi_rank, stage, requestss.arr[stage].size);
			fflush(stdout);
			#endif 
			if (requestss.arr[stage].size > 0){
				//mpi_wait_splitted(&requestss.arr[stage], &statusess.arr[stage]);
				MPI_CHECK( MPI_Waitall(requestss.arr[stage].size, requestss.arr[stage].arr, MPI_STATUSES_IGNORE)) ;
			}
			#if DEBUG
			printf("[%d | %d] After First Wait %d\n", mpi_rank, stage, requestss.arr[stage].size);
			fflush(stdout);
			#endif 
			if (comm_shape == SHAPE_TRIANG_2){
				if (requestss2.arr[stage].size > 0){
					MPI_CHECK( MPI_Waitall(requestss2.arr[stage].size, requestss2.arr[stage].arr, MPI_STATUSES_IGNORE)) ;
				}
				//mpi_wait_splitted(&requestss2.arr[stage], &statusess2.arr[stage]);
			}
			//MPI_CHECK( MPI_Wait( &requests.arr[ stage ], &statuses.arr[ stage ] ) );
			//if (comm_shape == SHAPE_TRIANG_2){
				//MPI_CHECK( MPI_Wait( &requests2.arr[ stage ], &statuses2.arr[ stage ] ) );
			//}
			#if DEBUG
			printf("[%d] After wait %d\n", mpi_rank, stage);
			fflush(stdout);
			#endif 
		}
		
		#if DEBUG > 1
		#if BUFFERED_B
		printf(
			"[%d | %d] remoteRange(%ld,%ld)\n", 
			mpi_rank, stage, 
			mat_B_stage.array[stage].rows.begin, 
			mat_B_stage.array[stage].rows.size
		);
		#else 
		printf(
			"[%d | %d] remoteRange(%ld,%ld)\n", 
			mpi_rank, stage, 
			mat_B_stage.rows.begin, 
			mat_B_stage.rows.size
		);
		#endif 
		fflush(stdout);
		//matrix_rangePrintf_rank_stage(
			mat_B_global, "BrB", 
			mpi_rank, stage, 0
		);
		fflush(stdout);
		#endif 
		#if BUFFERED_B
		if (mpi_rank == stage){
			matrix_rangePoint_to(&mat_B_stages.array[stage], mat_B_local);
		} 
		#else 
		if (mpi_rank == stage){
			matrix_rangePoint_to(&mat_B_stage, mat_B_local);
		} else {
			block_before_B 	= mat_B_stage.rows.begin * mat_B_global.columns.size;	
			mat_B_stage.mat = mat_B_global.mat + block_before_B;
			mat_B_stage.columns = mat_B_local.columns;
			//mat_B_stage.rows=f_mapping( mat_B_global.rows.size, mpi_procs, stage );
		}
		#endif 
		
		#if DEBUG > 1
		#if BUFFERED_B
		//matrix_rangePrintf_rank_stage(mat_B_stages.array[stage], "Bstg_postsarteco", mpi_rank, stage,0);
		//matrix_rangePrintf_rank_stage(mat_B_stage, "B", mpi_rank, stage, 0);
		for (i = 0; i < mat_B_stages.array[stage].rows.size; i ++){
			for (j = 0; j < mat_B_stages.array[stage].columns.size; j++){
				printf("[%d | stg %d] B[%ld][%ld] = %.3lf\n",
				 mpi_rank, stage, 
				 i+mat_B_stage.array[stage].rows.begin, j, 
				 mat_B_stage.array[stage].mat[i*mat_B_stage.array[stage].columns.size+j]);
			}
		}
		#else 
		//matrix_rangePrintf_rank_stage(mat_B_stage, "Bstg_postsarteco", mpi_rank, stage,0);
		//matrix_rangePrintf_rank_stage(mat_B_stage, "B", mpi_rank, stage, 0);
		for (i = 0; i < mat_B_stage.rows.size; i ++){
			for (j = 0; j < mat_B_stage.columns.size; j++){
				printf("[%d | stg %d] B[%ld][%ld] = %.3lf\n",
				 mpi_rank, stage, 
				 i+mat_B_stage.rows.begin, j, 
				 mat_B_stage.mat[i*mat_B_stage.columns.size+j]);
			}
		}
		#endif 
		fflush(stdout);
		#endif 
		///////////////////////////////////////////////////
		/// MATRIX MULTIPLICATION						///
		/// Only multiply the lower triangle to improve	///
		/// performance and to avoid 					///
		/// undesired initializations with 0s in 		///
		/// previous+boxes and previous+triang 			///
		/// combinations								///
		///////////////////////////////////////////////////
		/// timestamp 
		/// End 	| 	communication time
		/// Start	|	Multiplication time
		//////
		#if STATS
		comm1_t += set_time(comm1_time);
		gettimeofday(&mult_time, NULL);
		#endif 
		#if DEBUG > 1
			//matrix_rangePrintf_rank_stage(mat_A_local, "A_local just before product",mpi_rank, stage, 0);
		#endif 
		
		unsigned long int block_before_C = 0; //mat_A_local.rows.begin*mat_C_local.columns.size;
		#if DEBUG > 1
		//matrix_rangePrintf_rank_stage(mat_B_stage, "Bstg", mpi_rank,stage, 0);
		#endif 
		#if BUFFERED_B
		mm_triangular_triangular_product(
			mat_A_local, mat_B_stages.array[stage], 
			&mat_C_local, block_before_C, stage, mpi_rank
		);
		#else 
		mm_triangular_triangular_product(
			mat_A_local, mat_B_stage, &mat_C_local, 
			block_before_C, stage, mpi_rank
		);
		#endif 
		/* 5.3. End async. comm. and rotate A structures for next iteration */
		/* Scheme "owner" skips allocation and comm. in the last stage */
		if ( stage != mpi_procs - 1 ) {
			// Triangle part: buffer unmarshalling 
			if ( 
					comm_shape == SHAPE_TRIANG_2 
				&& 	comm_mechanism == MECH_BUFFER 
			) {
				index1 = 0;
				for ( 
					i=0; 
					i<mat_A_local.rows.size; 
					i++ 
				) {
					for ( j=mat_A_local.rows.begin; j<=i+mat_A_local.rows.begin; j++ ) {
						#if MULTIPLE_TRIANG_BUFF
						mat_B_stage.mat[ i*mat_B_global.columns.size+j ] = triang_buffers[stage].buffer[ index1 ];
						#else 
						#if BUFFERED_B
						mat_B_stages.array[stage].mat[ i*mat_B_global.columns.size+j ] = triang_buff.buffer[ index1 ];
						#else 
						mat_B_stage.mat[ i*mat_B_global.columns.size+j ] = triang_buff.buffer[ index1 ];
						#endif 
						#endif 
						index1++;
					}
				}
				// Free triang buffers
				//free( triang_buffers[stage].buffer ); ---- HERE
				free(triang_buff.buffer);
			}
		}
		#if BUFFERED_B
		//fprintf(stderr, "[%d | %d] I am using matrixB[%d] \n", mpi_rank, stage, stage);
		free(mat_B_stages.array[stage].mat);
		#endif 
	}

	///////
	/// Timestamp
	///	End 	| 	For time
	/// Start	| 	Comm2 - liberaciones
	//////
	#if STATS
	for_t = set_time(for_time);
	//aux_clock = MPI_Wtime();
	gettimeofday(&comm2_time, NULL);
	// Comm scheme "owner"
	// Precomputed derived data types can be freed now
	//free( mat_B_global.mat); --> que dice boxes que no le gusta este free
	comm2_t += set_time(comm2_time);
	#endif 
	// Time: Stop clock
	MPI_Barrier(MPI_COMM_WORLD);
	#if STATS
	total_t = set_time(total_time);
	#endif 

	/* 7. Output result */
	/*if ( output == OUTPUT_YES ) {
		print_output_results(
			mat_B_global.rows.size, 
			mat_B_global.columns.size, 
			mat_B_local.columns, 
			mpi_rank, 
			mat_C_local.mat, 
			mpi_procs, 
			"triang_triang");
		fflush(stdout);
	}*/
	if ( output == OUTPUT_YES ) {
		//print_output_results(mat_B_global.rows.size, mat_B_global.columns.size, mat_C_local.columns, mpi_rank, mat_C_local.mat, mpi_procs, "triang_triang");
		char file_name[63];
		sprintf( file_name, "./tests/outputs/output_triang_triang.%ld.%d.%d.dat", mat_B_global.rows.size, mpi_procs, mpi_rank );
		#if DEBUG
		printf("file_name = %s\n", file_name);
		#endif 
		//sprintf( file_name, "output.%d.%d.%d.dat", rows, mpi_procs, mpi_rank );
		FILE *file = fopen( file_name, "wb" );
		if ( file == NULL ) {
			char *info = calloc (50+strlen(file_name), sizeof(char));
			sprintf(info, "[triangular_triangular] Error: Opening the file %s\n", file_name );
			error_mpi_abort( MPI_COMM_WORLD, EXIT_FAILURE, mpi_rank, info);
		}
		fwrite( &mat_B_global.rows.size, sizeof(int), 1, file );
		fwrite( &mat_B_global.columns.size, sizeof(int), 1, file );
		output_range[0] = mat_C_local.rows.begin;
		output_range[1] = mat_C_local.rows.size+mat_C_local.rows.begin;
		output_range[2] = mat_C_local.columns.begin;
		output_range[3] = mat_C_local.columns.size +mat_C_local.columns.begin;
		fwrite( output_range, sizeof(int), 4, file );
		for ( i=0; i<mat_C_local.rows.size; i++ ) {
			fwrite( 
				&mat_C_local.mat[ i*mat_C_local.columns.size ], 
				sizeof(double), 
				mat_C_local.columns.size, file 
			);
		}
		fclose( file );
	}

	//#if DEBUG
	//matrix_rangePrintf_rank_stage(mat_C_local, "Local C", mpi_rank, -1, 0);
	fflush(stdout);
	//#endif
	/* 8. Write clocks */	
	comm1_t = Check_comm1 ? comm1_t : 0.0;
	//printf("checking error %.3lf\n", comm1_t);
	#if STATS
	set_ending_times(mpi_rank, total_t, precomp_t, comm1_t, comm2_t, mult_t, for_t);
	#endif 
	#if NORM 
	matrix_rangeInfinityNorm(mat_C_local, mat_B_global, MPI_COMM_WORLD, mpi_rank);
	#endif
	/* 9. End */
	MPI_Finalize();
	
	free( mat_A_local.mat );
//	free( mat_B_local.mat ); -- Este free se hace dentro del bucle :)
	free( mat_C_local.mat );
	return 0;
	
}	
