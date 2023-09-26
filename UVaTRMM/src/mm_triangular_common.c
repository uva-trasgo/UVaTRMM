/********************************************************************
 * @file mm_triangular_common										*
 * @addtogroup mm_triangular										*
 * @brief Common functions definitions for mm_triangular program.	*
 * @author Arturo González Escribano								*
 * @author Rocío Carratalá											*
 * @author Yuri Torres de la Sierra									*
 * @author Maria Inmaculada Santamaria Valenzuela					*
 * ******************************************************************/
#include <mm_triangular.h>

/************************************************
 * @section CONSTRUCTORS 						*
 * @brief Constructors for created structures.	*
 * **********************************************/
/***
 * @private buff
 * @brief Constructor
 * **/
Buff buff(int size, int mpi_rank){
	Buff buffer;
	buffer.size = size;
	buffer.buffer = calloc_abort(size, sizeof(double), "buff", MPI_COMM_WORLD, mpi_rank);
	return buffer;
}
/***
 * @private buff
 * @brief Constructor
 * **/
Sizes sizes(int size, int mpi_rank, MPI_Comm comm, char *name){
	Sizes buffer;
	buffer.size = size;
	buffer.buffer = calloc_abort(size, sizeof(size_t), name, comm, mpi_rank);
	return buffer;
}

/***
 * @private rangeSet
 * @brief Constructor
 * */
Range rangeSet(size_t begin, size_t size){
	Range range;
	range.begin = begin;
	range.size = size;
	return range;
}
/***
 * @private matrix_rangeCalloc
 * @brief Constructor
 * **/
MatrixRange matrix_rangeCalloc(size_t rows, size_t columns, const char *name, MPI_Comm comm, int mpi_rank){
	MatrixRange mat;
	mat.rows = rangeSet(0, rows);
	mat.columns = rangeSet(0, columns);
	matrix_rangeMat_calloc_abort(&mat, name, comm, mpi_rank);
	return mat;
}

/************************************************
 * @section MEMORY AUXILIAR FUNCTIONS			*
 * @brief Constructors for created structures.	*
 * **********************************************/

void freeNotNull(void *ptr){
	if (ptr != NULL){
		free(ptr);
	}
}

/***
 * @private calloc_abort
 * Soy un calloc pero con control de error:)
 * Adaptado a MPI
 * **/
void *calloc_abort(
	size_t __nmemb, 
	size_t __size, 
	const char *name, 
	MPI_Comm comm,
	int mpi_rank
){
	void *ptr = calloc(__nmemb, __size);
	fflush(stdin);
	fflush(stdout);
	
	char *info = (char*) calloc(ERROR_CHAR_SIZE, sizeof(char));
	if (ptr == NULL) {
		snprintf(info, ERROR_CHAR_SIZE, "%s%s", "Error Allocating memory for ", name);
		fprintf( stderr, "%s", info);
		fflush(stdout);
		error_mpi_abort(comm, MEMORY_ERROR, mpi_rank, info);
	}
	return ptr;
}

/***
 * @private calloc_abort
 * Soy un calloc pero con control de error:)
 * Adaptado a MPI
 * **/
void *malloc_abort(
	size_t __nmemb, 
	size_t __size, 
	const char *name, 
	MPI_Comm comm,
	int mpi_rank
){
	void *ptr = malloc(__nmemb* __size);
	fflush(stdin);
	fflush(stdout);
	
	char *info = (char*) calloc_abort(34+strlen(name), sizeof(char), "infostr", comm, mpi_rank);
	if (ptr == NULL) {
		sprintf(info, "%s%s", "Error Allocating memory for ", name);
		fprintf( stderr, "%s", info);
		fflush(stdout);
		error_mpi_abort(comm, MEMORY_ERROR, mpi_rank, info);
	}
	return ptr;
}

/***
 * @private calloc_abort
 * Soy un calloc pero con control de error:)
 * Adaptado a MPI
 * **/
void *realloc_abort(
	void *ptr,
	size_t __nmemb, 
	size_t __size, 
	const char *name, 
	MPI_Comm comm,
	int mpi_rank
){
	ptr = realloc(ptr, __nmemb* __size);
	fflush(stdin);
	fflush(stdout);
	
	char *info = (char*) calloc_abort(34+strlen(name), sizeof(char), "infostr", comm, mpi_rank);
	if (ptr == NULL) {
		sprintf(info, "%s%s", "Error Allocating memory for ", name);
		fprintf( stderr, "%s", info);
		fflush(stdout);
		error_mpi_abort(comm, MEMORY_ERROR, mpi_rank, info);
	}
	return ptr;
}



/***
 * @private Uses calloc for array with buffer's size
 * **/
int buff__double_arrayCalloc(Buff buffer, double **array, int mpi_rank){
	*array = calloc_abort(buffer.size, sizeof(double), "array", MPI_COMM_WORLD, mpi_rank);
	return buffer.size;
}

/***************************
 * @section BUFF FUNCTIONS *
 * *************************/
/***
 * @private buffCopy
 * Copies data of Orig into Dest
 * **/
void buffCopy(Buff Orig, Buff *Dest, int mpi_rank){
	int size = Orig.size;
	int i;
	*Dest = buff(size, mpi_rank);
	for(i = 0; i < size; i++){
		Dest -> buffer[i] = Orig.buffer[i];
	}
}
/***
 * @private buff_matrix_double_array
 * Copies buff data into mat
 * */
void buff__matrix_double_array(Range rows,int columns,double *mat,Buff buffer){
	int posa = 0;
	int pos = 0;
	int col, fil;
	int end = rows.size;
	int begin = rows.begin;
	for (fil = 0; fil < end; fil++){
		for (col = begin;col <= begin+fil; col++ ){
			posa =  fil*columns+col;
			mat[posa] = buffer.buffer[pos];
			pos++;
		}
	}
}

/***
 * @private buff__matrix_range
 * Copies buff data into mat
 * **/
void buff__matrix_range(Buff buffer, MatrixRange *mat){
	unsigned long int  posa = 0;
	unsigned long int  pos = 0;
	unsigned long int  col, fil;
	unsigned long int  end = mat -> rows.size;
	unsigned long int  begin = mat -> rows.begin;
	for (fil = 0; fil < end; fil++){
		for (col = begin; col <= begin+fil; col++ ){
	//for (fil = 1; fil < end; fil++){
	//	for (col = begin+1; col <= begin+fil; col++ ){
			posa =  fil*mat -> columns.size+col;
			mat -> mat[posa] = buffer.buffer[pos];
			pos++;
		}
	}
}

/***
 * @private set_triangular_buffer_data
 * fits triangular part of mat rows block 
 * into a buffer
 * **/
Buff set_triangular_buffer_data_deprecated_intento_corregir_superposicion(Range rows,int columns,double *mat, int mpi_rank){	
	Buff buffer;
	size_t i, j, pos = 0;
	double value = 0.0;
	//int size = ( rows.size + 1 ) * rows.size / 2;
	//
	int size = ( rows.size + 1 ) * rows.size / 2;
	size -= rows.size;
	//
	buffer = buff(size, mpi_rank);
	//for ( i=0; i< rows.size; i++ ) {
	for ( i=1; i< rows.size; i++ ) {
		for ( j = rows.begin+1; j <= rows.begin+i; j++ ) {
			value = mat[ i*columns+j ];
			buffer.buffer[ pos ] = value;
			pos ++;
		}
	}
	return buffer;
}


/***
 * @private set_triangular_buffer_data
 * fits triangular part of mat rows block 
 * into a buffer
 * **/
Buff set_triangular_buffer_data(
	Range rows,
	int columns,
	double *mat,
	int mpi_rank
){
	Buff buffer;
	double value = 0.0;
	unsigned long int i, j, pos = 0;
	unsigned long int size = ( rows.size + 1 ) * rows.size / 2;
	buffer = buff(size, mpi_rank);
	#pragma omp parallel for private(i, j) firstprivate(pos)
	for ( i=0; i< rows.size; i++ ) {
		for ( j = rows.begin; j <= rows.begin+i; j++ ) {
			value = mat[ i*columns+j ];
			buffer.buffer[ pos ] = value;
			pos ++;
		}
	}
	return buffer;
}

/************************************
 * @section MATRIX RANGE FUNCTIONS	*
 * **********************************/
/***
 * @private printf_matrix_ranges
 * Printing function
 * **/
 void printf_matrix_ranges(
	double* mat, 
	int row_begin, 
	int row_end, 
	int columns_begin, 
	int columns_end, 
	const char *name
){
	#if PRINT
	int i, j;
	int rows = row_end - row_begin;
	int columns = columns_end-columns_begin;
	for (i = 0; i < rows; i++){
		for (j = 0; j < columns; j++){
			printf(
				"%s[%d][%d] = %.3lf\n", 
				name,
				i+row_begin, 
				j+columns_begin,
				mat[i*columns+j]
			);
		}
	}
	#endif 
}


/***
 * @private printf_matrix_ranges_rank_stage
 * **/
void printf_matrix_ranges_rank_stage(
	double* mat, 
	int row_begin, 
	int row_end, 
	int column_begin, 
	int column_end, 
	const char *name,
	int mpi_rank,
	int stage,
	int block_before,
	_Bool check_triangular
){
	int i, j;
	int columns = column_end-column_begin;
	int rows = row_end-row_begin;
	#if PRINT
	printf(
		"[%d | %d] %s(%d,%d)x(%d,%d)\n", 
		mpi_rank, stage, name,
		row_begin, row_end-1,
		column_begin, column_end-1
	);
	fflush(stdout);
	#endif 
	for (i = 0; i < rows; i++){
		for (j = 0; j < columns; j++){
			#if PRINT
			printf(
				"[%d | %d] %s[%d][%d] = %.3lf\n", 
				mpi_rank, stage, name,
				i+row_begin, 
				j+column_begin,
				mat[i*columns+j+block_before]
			);
			fflush(stdout);
			#endif
			if(check_triangular){
				if (
					mat[i*columns+j+block_before]==0 
					&& 
					(i+row_begin >= j+column_begin)
				){
				printf(
					"[%d | %d] :( %s[%d][%d] = 0 :(\n", 
					mpi_rank, stage, name, i+row_begin, j+column_begin
				);
				printf(
					"[%d | %d] %s(%d,%d)x(%d,%d)\n", 
					mpi_rank, stage, name,
					row_begin, row_end-1,
					column_begin, column_end-1
				);
				error_mpi_abort( MPI_COMM_WORLD, EXIT_FAILURE, mpi_rank, "Not a triangular matrix | Cuidado. numerror = EXIT FAILURE");			
				}
			}
		}
	}
}

/***
 * @private matrix_rangePrintf_rank_stage
 * Printing function
 * **/
void matrix_rangePrintf_rank_stage(MatrixRange mat,const char *name,int mpi_rank,int stage,int block_before){
	#if PRINT
	int i, j;
	#endif 
	printf("[%d | %d] ---> %s --- (%ld,%ld)x(%ld,%ld) -- \n", 
		mpi_rank, stage, 
		name,
		mat.rows.begin, 
		mat.rows.size+mat.rows.begin,
		mat.columns.begin, 
		mat.columns.begin+mat.columns.size
	);
	
	#if PRINT
	fflush(stdout);
	for (i = 0; i < mat.rows.size; i++){
		for (j = 0; j < mat.columns.size; j++){
			printf("[%d | %d] %s[%ld][%ld] = %.3lf\n", mpi_rank, stage, 
			name,i+mat.rows.begin, 	
			j+mat.columns.begin,mat.mat[i*mat.columns.size+j+block_before]);
		fflush(stdout);
		}
	}
	printf("[%d | %d] --- %s --->\n", mpi_rank, stage, name);
	fflush(stdout);
	#endif 
}

/*******************************
 * @section REQUESTS FUNCTIONS *
 * *****************************/
/***
 * @private requestsMakeNULL
 * Makes all values of requests' array MPI_REQUEST_NULL
 * **/
void requestsMakeNull(Requests *requests){
	int req;
	for (req = 0; req < requests -> size; req++){
		requests -> arr[req] = MPI_REQUEST_NULL;
	}
}

/***
 * @private requestsCalloc_abort 
 * Constructor
 * **/
Requests requestsCalloc_abort(size_t size, MPI_Comm comm, int mpi_rank){
	Requests requests;
	requests.size = size;
	requests.arr = calloc_abort(size,sizeof(MPI_Request), "requests", comm, mpi_rank);
	requestsMakeNull(&requests);
	return requests;
}


/***
 * @private requestsMalloc_abort 
 * Constructor
 * **/
Requests requestsMalloc_abort(size_t size, MPI_Comm comm, int mpi_rank){
	Requests requests;
	requests.size = size;
	requests.arr = malloc_abort(size,sizeof(MPI_Request), "requests", comm, mpi_rank);
	return requests;
}

/***
 * @private requestsCalloc_abort 
 * Constructor
 * **/
Statuses statusesCalloc_abort(size_t size, MPI_Comm comm, int mpi_rank){
	Statuses statuses;
	statuses.size = size;
	statuses.arr = calloc_abort(size,sizeof(MPI_Status), "statuses", comm, mpi_rank);
	return statuses;
}


/****************************************************************************************************
 * @section DATATYPE CREATION																		*
 * @brief Functions to create derived datatypes for rectangular and triangular parts of matrices.	*
 * @private Remember | types are done supossing that our data buffer starts at the beginning row.	*
 * **************************************************************************************************/
/***
 * @private create_type_bbox
 * * Type name
 * * * boxes
 * * Description
 * * * Bounding box of rectangle + triangle
 * **/
MPI_Datatype create_type_bbox( Range range, unsigned long int columns ) {
	MPI_Datatype dtype;
	MPI_CHECK( 
		MPI_Type_vector( 
			range.size, 
			range.begin + range.size, 
			columns, 
			MPI_DOUBLE, 
			&dtype 
		) 
	);
	return dtype;
}

MPI_Datatypes mpi_datatypesCalloc(size_t total, char *name, MPI_Comm comm, int mpi_rank){
	MPI_Datatypes type_array = {total, NULL};
	type_array.arr	=calloc_abort(total, sizeof(MPI_Datatype), name, comm, mpi_rank);
	type_array.sizes=calloc_abort(total, sizeof(size_t), name, comm, mpi_rank);
	return type_array;
}

void mpi_datatypesFree(MPI_Datatypes *types){
	size_t t;
	for (t = 0; t < types->total; t++){
		if (types->sizes[t] > 0 && types->arr[t] != MPI_DOUBLE){
			MPI_Type_free(&types->arr[t]);
		}
	}
	
	if  (types -> sizes != NULL){
		free(types->sizes);
	}
	
	if  (types -> arr != NULL){
		free(types-> arr);
	}

	types->total=0;
}

void reset_int_arr(int **arr, int total, char *name, MPI_Comm comm, int mpi_rank){
	if (*arr != NULL){
			free(*arr);
			*arr = NULL;
	}
	*arr = (int*) calloc_abort(total, sizeof(int), "sub block len", comm, mpi_rank);
}

/***
 * @private create_type_bboxARRAY
 * * Type name
 * * * triangle
 * * Description
 * * * Triangular part of panel
 * **/
MPI_Datatypes create_type_bboxArray( Range range, unsigned long int columns, size_t block_limit_quotient_denominator, int mpi_rank, MPI_Comm comm) {
	#if DEBUG
	printf("Create type bbox array\n"); 
	fflush(stdout);
	#endif 
	MPI_Datatypes dtypes = {0, NULL, NULL};
	//-- Auxiliar
	int i, block;
	unsigned long int num_intermediate_blocks;
	unsigned int block_height = 0;
	unsigned int last_block_height = 0;
	unsigned int block_max_size;
	//-- Indexed type variables
	int count; 		//Number of rows in block
	int *sub_block_displ=NULL; 	//Displacements from starting point
	int *sub_block_len=NULL; 	//Lengths
	int width, height;
	//-- Check sizes and num blocks
	#if DEBUG
	printf("block_limit_quotient_denominator is %ld || BLOCK_LIMIT %d\n", block_limit_quotient_denominator, BLOCK_LIMIT);
	fflush(stdout);
	#endif 
	if (block_limit_quotient_denominator == 0 ){
		#if DEBUG
		printf("block_limit_quotient_denominator is 0\n");
		#endif 
		num_intermediate_blocks=0;
		last_block_height = range.size;
		//error_mpi_abort(comm, TYPE_ARRAY_CREATION_ERROR, mpi_rank, "invalid denominator for datatype array generation (bbox)");
	} else {
		block_max_size=MAX(1,BLOCK_LIMIT/block_limit_quotient_denominator);
	}
	height = range.size;
	width = range.begin+range.size;
	#if DEBUG
	printf("[%d] Height %d Width %d\n", mpi_rank, height, width);
	fflush(stdout);
	#endif
	count = height*width; 

	if (block_limit_quotient_denominator > 0){
		block_height=block_max_size/width;
		if (block_height == 0){
			num_intermediate_blocks=0;
			last_block_height = range.size;
		} else {
			num_intermediate_blocks=MAX(0, range.size/block_height);
			last_block_height=height-(block_height*num_intermediate_blocks);
		}
	} else {
		num_intermediate_blocks=0;
		block_height=range.size;
	}
	if (last_block_height > 0){
		dtypes = mpi_datatypesCalloc(num_intermediate_blocks+1, "datatypes array bbox with last block", comm, mpi_rank);
	} else {
		dtypes = mpi_datatypesCalloc(num_intermediate_blocks, "datatype array bbox without last block", comm, mpi_rank);
	}
	count = block_height;
	#if DEBUG
		printf("[%d] Range (%ld,%ld) %ld intermediate blocks || intersize %d || lastsize %d || total %ld\n",mpi_rank, range.begin, range.begin+range.size-1, num_intermediate_blocks, block_height*width, last_block_height*width, dtypes.total);
		fflush(stdout);
	#endif
	//#pragma omp parallel for private(block)
	#if DEBUG
	printf("[%d | (%d, %d) ] Count %d Num Intermediate Blocks %d Last height %d Intermediate height %d Block max size %d || %intermediate Blocks\n", 
		mpi_rank, range.begin, range.size, count, num_intermediate_blocks, last_block_height, block_height, block_max_size
	);
	#endif 
	if (count > 0){
		#if DEBUG
		printf("Count > 0 %ld Width %d\n", count, width);
		fflush(stdout);
		#endif 
		for (block=0; block < num_intermediate_blocks; block++){
			reset_int_arr(&sub_block_len, count, "sub_block_len", comm, mpi_rank);
			reset_int_arr(&sub_block_displ, count, "sub_block_displ", comm, mpi_rank);
			//-- Setup lenghts and displacement for each row
			//#pragma omp parallel for private(i) v- cuidado con .sizes y demás. Daba problemas
			for ( i=0; i < count; i++ ) {
				sub_block_len[ i ] 	 = width;
				sub_block_displ[ i ] = ((block_height*block)+i)*columns;
				dtypes.sizes[block]+=sub_block_len[i];
				#if DEBUG
				printf("[%d | %d | %d] len %d displ %d\n", mpi_rank, block, i, sub_block_len[i], sub_block_displ[i]);
				fflush(stdout);
				#endif
			}
			MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block) );
			MPI_CHECK( MPI_Type_commit( dtypes.arr + block) );
		}
	}
		#if DEBUG
		printf("%d width|| height || Last block height%d\n", width, block_height, last_block_height);
		fflush(stdout);
		#endif 
	if (last_block_height > 0){
		count = last_block_height;
		reset_int_arr(&sub_block_len, count, "sub_block_len", comm, mpi_rank);
		reset_int_arr(&sub_block_displ, count, "sub_block_displ", comm, mpi_rank);
		//-- Setup lenghts and displacement for each row
		//#pragma omp parallel for private(i) -> Daba problemas
		for ( i=0; i < count; i++ ) {
			sub_block_len[ i ] = width;
			sub_block_displ[ i ] = ((block_height*num_intermediate_blocks)+i)*columns;
			dtypes.sizes[num_intermediate_blocks]+=sub_block_len[i];
		}
		MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr +  num_intermediate_blocks) );
		MPI_CHECK( MPI_Type_commit( dtypes.arr + num_intermediate_blocks) );

	}
	return dtypes;
}







/***
 * @private create_type_rectangle
 * * Type name
 * * * rectangle
 * * Description
 * * * Rectangular part of panel
 * **/
MPI_Datatype create_type_rectangle( Range range, unsigned long int columns ) {
	MPI_Datatype dtype;
	MPI_CHECK( 
		MPI_Type_vector( 
			range.size,		//starting count
			range.begin, 	//Length
			columns, 		//Space from init to init (stride)
			MPI_DOUBLE, 	
			&dtype 
		) 
	);
	return dtype;
}


/***
 * @private get_rows_before_trapezoid
 * @brief Gets rows before of panel with block_index id
 * */
unsigned long int get_rows_before_trapezoid(double ideal_block_size, int block_index, int rows_begin){
    double numerator = -1+sqrt(1+8*(block_index*ideal_block_size+((rows_begin-1)*rows_begin)));
    unsigned long int rows_before = (numerator > 0) ? (floor(numerator)/2) : 0;
    return rows_before;
}

void split_trapezoid(Range range, unsigned long int columns, size_t block_limit_quotient_denominator, int mpi_rank, MPI_Comm comm, Sizes *heights, Sizes *disp) { 
	unsigned long int num_intermediate_blocks=0; 
	unsigned long int count=0;
	//int block_height;
	int block_max_size=0;
	int height=0;
	int width_rectangle=0;
	int last_block_height=0; 
	int last_block_rows_begin=0;
	int last_block_rows_end=0;
	int block_rows_begin, block_rows_end=0;
	unsigned int block;
	height = range.size;
	width_rectangle = range.begin;
	count = height*width_rectangle+((height*(height+1)/2)); 
	
	if (block_limit_quotient_denominator ==0 ){
		num_intermediate_blocks=0;
		
	} else {
		block_max_size=BLOCK_LIMIT/block_limit_quotient_denominator;
		num_intermediate_blocks=count/block_max_size;

	}
	#if DEBUG
	printf("[%d | (b%ld, s%ld)] Intermediate Blocks %ld Max Size %d Total Elements %ld\n", mpi_rank, range.begin, range.size, num_intermediate_blocks, block_max_size, count);
	#endif
	heights->size=0;
	heights->buffer=NULL;
	disp->size=0;
	disp->buffer=NULL;
	size_t deleted = 0;
	if (num_intermediate_blocks > 0){
		*heights 	= sizes(num_intermediate_blocks, mpi_rank, comm,  "height array split trapezoid with last block");
		*disp 		= sizes(num_intermediate_blocks, mpi_rank, comm, "disp array split trapezoid with last block");
		heights->size = num_intermediate_blocks;
		disp->size=num_intermediate_blocks;
		unsigned int aux;
		//#pragma omp parallel for private(block)
		for (block=0; block < num_intermediate_blocks; block++){
			if (block==0){
				block_rows_begin=range.begin;
			} else {
				if (heights->buffer[block-deleted-1] > 0){
					block_rows_begin=block_rows_end+1;
				} else {
					block_rows_begin=block_rows_end;
				}
			}
			block_rows_end=get_rows_before_trapezoid(block_max_size, block+1, range.begin);
			aux=MAX(0,range.begin+range.size-1);
			block_rows_end=MIN(block_rows_end,aux);			
			heights->buffer[block-deleted] = MAX(0,block_rows_end-block_rows_begin+1);
		
			if (heights->buffer[block-deleted] > 0){
				disp->buffer[block-deleted] = (block_rows_begin > range.begin) ? block_rows_begin + range.begin : 0;
				#if DEBUG
				printf("[%d | (b%ld,s%ld) | %d] Intermediate Block (%d,%d) | %ld height %ld displ\n", 
					mpi_rank, range.begin, range.size,
					block, block_rows_begin, block_rows_end, heights->buffer[block-deleted], disp->buffer[block-deleted] );
				#endif 
			} else {
				deleted ++;
				heights->size-=1;
				disp->size -= 1;
			}
			//-- Setup lenghts and displacement for each row
		}
			if (heights -> buffer[num_intermediate_blocks-1] > 0){
				last_block_rows_begin = block_rows_end + 1;
			} else {
				last_block_rows_begin = block_rows_end;
			}
		
	} else {
		
		last_block_rows_begin=range.begin;
		last_block_rows_end=range.size+range.begin-1;
		last_block_height=range.size;
	}
	if (block_limit_quotient_denominator == 0 || (heights->size == 0 )){
		num_intermediate_blocks=0;
		last_block_rows_begin=range.begin;
		last_block_rows_end=range.size+range.begin-1;
		last_block_height=range.size;
	} else {
		last_block_rows_end=range.size+range.begin-1;
		//last_block_height=MAX(0,last_block_rows_end-last_block_rows_begin);
		last_block_height=MAX(0,last_block_rows_end-last_block_rows_begin+1);
	}
	if (last_block_height > 0){
		
		if (block_limit_quotient_denominator == 0 || (heights->size < 2 )){
			heights->size+=1;
			disp->size+=1;
			*heights 	= sizes(1, mpi_rank, comm,  "height array split trapezoid with last block");
			*disp 		= sizes(1, mpi_rank, comm, "disp array split trapezoid with last block");
			heights->buffer[0] = range.size;
			disp->buffer[0] = 0;
		} else {
			if (last_block_rows_begin > range.begin && last_block_rows_begin < range.begin+range.size-1){
				heights->size+=1;
				disp->size+=1;
				heights -> buffer=(size_t*)realloc(heights->buffer, sizeof(size_t)*heights->size);
				disp -> buffer=(size_t*)realloc(disp->buffer, sizeof(size_t)*disp->size);
				heights->buffer[heights->size-1] = last_block_height;
				disp->buffer[disp->size-1] =  last_block_rows_begin-range.begin;
			} 
		}
		#if DEBUG
		printf("[%d | %ld ] Last Block (%d,%d) Range (%ld,%ld) | %ld height | %ld displ | datatypes %ld \n", 
			mpi_rank, num_intermediate_blocks, 
			last_block_rows_begin, last_block_rows_end, 
			range.begin, range.begin+range.size-1, 
			heights->buffer[num_intermediate_blocks], disp->buffer[num_intermediate_blocks], heights->size);
		fflush(stdout);
		#endif 
	}
	#if DEBUG
	printf("[%d] (%ld, %ld) last block rows end  = %d Begin = %d Height %d | Total heights: %ld\n", 
	mpi_rank, range.begin, range.size, last_block_rows_end, last_block_rows_begin, last_block_height, heights->size);
	for (block = 0; block < heights -> size; block++) {
		printf(
			"[%d] (%ld, %ld) Height[%d/%ld] = %ld Disp %ld\n",
		 		mpi_rank, range.begin, range.size,
			 	block, heights->size,  heights->buffer[block], disp->buffer[block]
			 );
		fflush(stdout); 
	}
	
	fflush(stdout);
	#endif 
}



/***
 * @private create_type_rectangleARRAY
 * * Type name
 * * * triangle
 * * Description
 * * * Triangular part of panel
 * **/
MPI_Datatypes create_type_rectangleArray( 
	Range range, unsigned long int columns, 
	int mpi_rank, MPI_Comm comm, Sizes heights, Sizes displacements
) {
	MPI_Datatypes dtypes = {0, NULL, NULL};
	//-- Auxiliar
	int i, block;
	//unsigned long int num_intermediate_blocks;
	//int block_height, block_max_size;
	//int last_block_height, last_block_rows_begin, last_block_rows_end;
	//-- Indexed type variables
	int count; 		//Number of rows in block
	int *sub_block_displ=NULL; 	//Displacements from starting point
	int *sub_block_len=NULL; 	//Lengths
	//int width_rectangle;
	int block_rows_begin;
	//int block_rows_end, height;
	//-- Check sizes and num blocks

	dtypes=mpi_datatypesCalloc(heights.size, "mpi rectangle array datatypes", comm, mpi_rank);
	size_t up_space, left_space;
	left_space=0;
	//#pragma omp parallel for private(block)
	if (range.begin > 0 && range.size > 0){
		for (block=0; block < heights.size; block++){
			block_rows_begin=displacements.buffer[block];
			count = heights.buffer[block];
			reset_int_arr(&sub_block_len, count, "sub_block_len", comm, mpi_rank);
			reset_int_arr(&sub_block_displ, count, "sub_block_displ", comm, mpi_rank);
			dtypes.sizes[block]=0;
			//-- Setup lenghts and displacement forleft_space each row
			//#pragma omp paralleft_spacelel for private(i)
			if (count > 0){
				for ( i=0; i < count; i++ ) {
					up_space=(block_rows_begin+i) * columns;
					sub_block_len[ i ] 	= range.begin;
					sub_block_displ[ i ] = up_space+left_space;
					dtypes.sizes[block]+=sub_block_len[i];
					#if DEBUG 
					printf("[%d | (b%ld,s%ld)| RecBlock %d | i %d | %d ] rows before %d displ %ld Up %ld LeftSpace %d len %d \n", 
					mpi_rank, range.begin, range.size, block, i, block_rows_begin, sub_block_displ[i], up_space, left_space, sub_block_len[i], block_rows_begin);
					fflush(stdout);
					#endif 
					MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block) );
					MPI_CHECK( MPI_Type_commit( dtypes.arr + block) );
				}
				#if DEBUG
				printf("[%d | (%ld,%ld)] Rectangle %d Count %d Size %ld row begin %d\n", 
				mpi_rank, range.begin, range.begin+range.size, dtypes.arr[block],count, dtypes.sizes[block], block_rows_begin);
				fflush(stdout);
				#endif
			} 
		}
	}
	return dtypes;
}
/***
 * @private create_type_triangle
 * * Type name
 * * * triangle
 * * Description
 * * * Triangular part of panel
 * **/
MPI_Datatype create_type_triangle( Range range, unsigned long int columns ) {
	MPI_Datatype dtype;
	//-- Auxiliar
	size_t i;
	//-- Indexed type variables
	//-- propuesta macu revisando linea comun: size_t count = (range.size > 0) ? range.size -1 : 0; //Num or elements
	size_t count = range.size;
	int displ[ count ]; //Displacements from starting point
	int blocklen[ count ]; //Lengths
	//-- Setup lenghts and displacement for each row
	#pragma omp parallel for private(i)
	for ( i=0; i < count; i++ ) {
		blocklen[ i ] = i+1;
		displ[ i ] = range.begin + i * columns;
	}
	MPI_CHECK( MPI_Type_indexed( count, blocklen, displ, MPI_DOUBLE, &dtype) );
	return dtype;
}



/***
 * @private create_type_triangleARRAY
 * * Type name
 * * * triangle
 * * Description
 * * * Triangular part of panel
 * **/

MPI_Datatypes create_type_triangleArray( 
	Range range, unsigned long int columns, int mpi_rank, MPI_Comm comm, Sizes heights, Sizes displacements,
	int stage
) {
	MPI_Datatypes dtypes = {0, NULL};
	//-- Auxiliar
	int i, block;
	size_t left_space;
	size_t up_space;
	//unsigned long int num_intermediate_blocks;
	//int block_height;
	//int block_max_size;
	//int last_block_height;
	//int last_block_rows_begin;
	//int last_block_rows_end;
	//-- Indexed type variables
	int count; 		//Number of rows in block
	int *sub_block_displ=NULL; 	//Displacements from starting point
	int *sub_block_len=NULL; 	//Lengths
	//int width_rectangle;
	int block_rows_begin;
	//int block_rows_end;
	//int height;
	//-- Check sizes and num blocks
	dtypes=mpi_datatypesCalloc(heights.size, "mpi rectangle array datatypes", comm, mpi_rank);
	#if DEBUG
	rangePrintf_rank_stage("triangArray", range, mpi_rank, stage);
	fflush(stdout);
	#endif
	#if DEBUG
	printf("[%d] Dtypes %ld\n", 
		mpi_rank, dtypes.total
	);
	fflush(stdout);
	#endif 
	//#pragma omp parallel for private(block)
	for (block=0; block < dtypes.total; block++){
		block_rows_begin=displacements.buffer[block];
		count = heights.buffer[block];
		if (count > 0){
			reset_int_arr(&sub_block_len, count, "sub_block_len", comm, mpi_rank);
			reset_int_arr(&sub_block_displ, count, "sub_block_displ", comm, mpi_rank);
			#if DEBUG
			printf("[%d | %d] TriangBlock %d | begin %d count %d || size %ld\n", mpi_rank, stage, block, block_rows_begin, count, block_rows_begin*columns+(count*(count+1))/2);
			fflush(stdout);
			#endif 
			//-- Setup lenghts and displacement for each row
			//#pragma omp parallel for private(i)
			left_space = (range.begin == 0) ? 0 : range.begin;//range.begin+1;
			//up_space;
			dtypes.sizes[block]=0;
			for ( i=0; i < count; i++ ) {
				up_space=block_rows_begin+i*columns;// + block_rows_begin;
				sub_block_len[ i ] = i+1;
				sub_block_displ[ i ] = up_space + left_space;
				#if DEBUG
				printf("[%d | %d | %d] up %ld len %d disp %d\n",
					mpi_rank, block, i, up_space, sub_block_len[i], sub_block_displ[i]
				);
				fflush(stdout);
				#endif 
				dtypes.sizes[block]+=sub_block_len[i];
				#if DEBUG > 1
				printf("[%d | %d | TriangBlock Range (%ld,%ld) %d | i %d ] rows before %d displ %d Up %ld LeftSpace %ld len %d\n", 
				mpi_rank, stage, range.begin, range.begin+range.size-1, block, i, block_rows_begin, sub_block_displ[i], up_space, left_space, sub_block_len[i]);
				fflush(stdout);
				#endif 
			}
			MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block) );
			MPI_CHECK( MPI_Type_commit( dtypes.arr + block) );
			#if DEBUG > 1
				printf("[%d | %d | TriangBlock Range (%ld,%ld) %d ] Size %ld\n", 
				mpi_rank, stage, range.begin, range.begin+range.size-1, block, sub_block_len[i]);
				fflush(stdout);
			#endif 
		} else {
			#if DEBUG > 1
				printf("[%d | %d | TriangBlock Range (%ld,%ld) %d ] Size %ld\n", 
				mpi_rank, stage, range.begin, range.begin+range.size-1, block, 0);
				fflush(stdout);
			#endif 
			dtypes.sizes[block]=0;
		}
	}

	return dtypes;
}
/*** 
 * @private create_type_combined
 * * Type name: triang2type [trapezoid rectangle + triangle]
 * * Description
 * Combines rectangle and triangle datatypes
 * */
MPI_Datatype create_type_combined( Range range, unsigned long int columns ) {
	//-- Datatypes
	MPI_Datatype dtype;
	MPI_Datatype rectangle 	= create_type_rectangle( range, columns );
	MPI_Datatype triangle 	= create_type_triangle( range, columns );
	
	//-- Setting up dtype (struct)
	int lengths[2] = { 1, 1 }; //num of elements of each subtype
	MPI_Aint displ[2] = { 0, 0 }; //Displacement for each subtype
	MPI_Datatype types[2] = { rectangle, triangle };
	MPI_CHECK( MPI_Type_create_struct( 2, lengths, displ, types, &dtype) );
	
	//-- Free auxiliar datatypes
	MPI_Type_free( &rectangle );
	MPI_Type_free( &triangle );
	return dtype;
}

/***
 * @private create_type_rectangleARRAY
 * * Type name
 * * * triangle
 * * Description
 * * * Triangular part of panel
 * **/
MPI_Datatypes create_type_rectangleArrayCombined( 
	Range range, unsigned long int columns, 
	int mpi_rank, MPI_Comm comm, Sizes heights, Sizes displacements
) {
	#if DEBUG
	printf("[%d] Inside get rectangle vectoret combined\n", mpi_rank);
	fflush(stdout); 
	#endif 
	MPI_Datatypes dtypes = {0, NULL, NULL};
	//-- Auxiliar
	int i, block;
	//unsigned long int num_intermediate_blocks;
	//int block_height, block_max_size;
	//int last_block_height, last_block_rows_begin, last_block_rows_end;
	//-- Indexed type variables
	int count; 		//Number of rows in block
	int *sub_block_displ= NULL; 	//Displacements from starting point
	int *sub_block_len	= NULL; 	//Lengths
	//int width_rectangle;
	int block_rows_begin;
	//int block_rows_displ;
	//int block_rows_end, height;
	//-- Check sizes and num blocks
	#if DEBUG >1
	printf("[%d] Before dtypes  calloc\n", mpi_rank);
	fflush(stdout);
	#endif 
	dtypes=mpi_datatypesCalloc(heights.size, "mpi rectangle array datatypes", comm, mpi_rank);
	#if DEBUG > 1
	printf("[%d] After dtypes  calloc\n", mpi_rank);
	fflush(stdout);
	#endif 
	size_t up_space, left_space;
	size_t total = dtypes.total;
	size_t deleted = 0;
	size_t count_add = 0;
	left_space=0;
	//#pragma omp parallel for private(block)
	#if DEBUG > 1
	printf("[%d] Dtypes.size %ld\n", mpi_rank, dtypes.total);
	fflush(stdout);
	size_t height, width;
	#endif
	//if (range.begin > 0 && range.size > 0){
		#if DEBUG
		printf("[%d] Inside rectangle if\n", mpi_rank);
		fflush(stdout);
		#endif 
		size_t height, width;
		for (block=0; block < total; block++){
			//block_rows_displ=displacements.buffer[block];
			block_rows_begin=count_add;//block_rows_displ/columns;
			//count = heights.buffer[block];
			//height = range.size;
			height = heights.buffer[block-deleted];
			width  = range.begin;
			count  = height;
				if (count > 0 && (block_rows_begin > 0 || range.begin > 0)){
				reset_int_arr(&sub_block_len, count, "sub_block_len", comm, mpi_rank);
				reset_int_arr(&sub_block_displ, count, "sub_block_displ", comm, mpi_rank);
				dtypes.sizes[block-deleted]=0;
				//dtypes.sizes[block]=0;
				//-- Setup lenghts and displacement for each row
				//#pragma omp parallel for private(i)
				#if DEBUG
				printf("[%d] Block Rows Begin %d | Count %d\n", 
				mpi_rank, block_rows_begin, count); fflush(stdout) ;
				#endif 
				for ( i=0; i < count; i++ ) {
					up_space=(block_rows_begin+i) * columns;
					sub_block_len[ i ] 	= width;
					sub_block_displ[ i ] = up_space+left_space;
					dtypes.sizes[block-deleted]+=sub_block_len[i];
					//dtypes.sizes[block]+=sub_block_len[i];
					#if DEBUG 
					printf("[%d | (b%ld,s%ld)| RecBlock %d | i %d ] rows before %d displ %ld Up %ld LeftSpace %d len %ld \n", 
					mpi_rank, range.begin, range.size, block, i,
					 block_rows_begin, 
					 sub_block_displ[i],
					  up_space, left_space, sub_block_len[i]);
					fflush(stdout);
					#endif 
					MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block-deleted) );
					MPI_CHECK( MPI_Type_commit( dtypes.arr + block-deleted) );
					//MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block) );
					//MPI_CHECK( MPI_Type_commit( dtypes.arr + block) );
				}
				#if DEBUG
				printf("[%d | (%ld,%ld)] Rectangle %d Count %d Size %ld row begin %d\n", 
				mpi_rank, range.begin, range.begin+range.size-1, block,count, dtypes.sizes[block], block_rows_begin);
				fflush(stdout);
				#endif
				count_add += count;
			} else {
				dtypes.total -= 1;
				deleted += 1;
				dtypes.sizes[block]=0;
			} 
		}
	/*} else {
		//Aqui no hay problema... asi que estoy segura de que los rectanuglos estarian bien hechos. Pruebo la otra opcion de particionado.
		printf("[%d] Ensuring size 0\n", mpi_rank);
		fflush(stdout);
		for (block=0; block < dtypes.total;block++){
			dtypes.sizes[block]=0;
		}
		printf("[%d] Ensured size 0\n", mpi_rank);
		fflush(stdout);
	}
	printf("[%d] Despues del if \n", mpi_rank);
	fflush(stdout)*/;
	#if DEBUG
	printf("[%d | (%ld,%ld)] Rectangle Dtypes.total %d \n", 
				mpi_rank, range.begin, range.begin+range.size-1, 
				dtypes.total);
				fflush(stdout);
	#endif 
	return dtypes;
}

MPI_Datatypes create_type_triangleArrayCombined__( 
	Range range, unsigned long int columns, int mpi_rank, MPI_Comm comm, Sizes heights, Sizes displacements,
	int stage
) {
	MPI_Datatypes dtypes = {0, NULL};
	//-- Auxiliar
	int i, block;
	size_t left_space;
	size_t up_space;
	//-- Indexed type variables
	int count; 		//Number of rows in block
	int *sub_block_displ=NULL; 	//Displacements from starting point
	int *sub_block_len=NULL; 	//Lengths
	//int width_rectangle;
	int block_rows_begin;
	//int block_rows_end;
	//int height;
	//-- Check sizes and num blocks
	dtypes=mpi_datatypesCalloc(heights.size, "mpi rectangle array datatypes", comm, mpi_rank);
	#if DEBUG
	rangePrintf_rank_stage("triangArray", range, mpi_rank, stage);
	fflush(stdout);
	#endif

	//#pragma omp parallel for private(block)
	for (block=0; block < dtypes.total; block++){
		block_rows_begin=displacements.buffer[block];
		count = heights.buffer[block];
		if (count > 0){
			reset_int_arr(&sub_block_len, count, "sub_block_len", comm, mpi_rank);
			reset_int_arr(&sub_block_displ, count, "sub_block_displ", comm, mpi_rank);
			//-- Setup lenghts and displacement for each row
			//#pragma omp parallel for private(i)
			left_space = block_rows_begin;
			//up_space;
			dtypes.sizes[block]=0;
			#if DEBUG
			printf(
				"[%d | %d  Range (b%ld,s%ld) ] TriangBlock %d | begin %d count %d || Left space %ld || size %d\n", 
				mpi_rank, stage, range.begin, range.size,
				block, block_rows_begin, count, left_space, (count*(count+1))/2
			);
			fflush(stdout);
			#endif 
			for ( i=0; i < count; i++ ) {
				up_space=(block_rows_begin+i)*columns;// + block_rows_begin;
				sub_block_len[ i ] = i+1;
				sub_block_displ[ i ] = up_space + left_space;
				dtypes.sizes[block]+=sub_block_len[i];
				#if DEBUG > 1
				printf("[%d | %d | TriangBlock Range (%ld,%ld) %d | i %d ] rows before %d displ %d Up %ld LeftSpace %ld len %d\n", 
				mpi_rank, stage, range.begin, range.begin+range.size-1, block, i, block_rows_begin, sub_block_displ[i], up_space, left_space, sub_block_len[i]);
				fflush(stdout);
				#endif 
			}
			#if DEBUG
			printf(
				"[%d | %d  Range (b%ld,s%ld) ] TriangBlock %d | First Len %d Last Len %d | First disp %d Last disp %d\n", 
				mpi_rank, stage, range.begin, range.size, block,
				sub_block_len[0],
				sub_block_len[count-1],
				sub_block_displ[0],
				sub_block_displ[count-1]
			);
			fflush(stdout);
			#endif 
			MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block) );
			MPI_CHECK( MPI_Type_commit( dtypes.arr + block) );
			#if DEBUG > 1
				printf("[%d | %d | TriangBlock Range (%ld,%ld) %d ] Size %ld\n", 
				mpi_rank, stage, range.begin, range.begin+range.size-1, block, sub_block_len[i]);
				fflush(stdout);
			#endif 
		} else {
			
			dtypes.sizes[block]=0;
		}
	}

	return dtypes;
}


MPI_Datatypes create_type_triangleArrayCombined( 
	Range range, unsigned long int columns, int mpi_rank, MPI_Comm comm, Sizes heights, Sizes displacements,
	int stage
) {
	MPI_Datatypes dtypes = {0, NULL};
	//-- Auxiliar
	int i, block;
	size_t left_space;
	size_t up_space;
	//-- Indexed type variables
	int count; 		//Number of rows in block
	int *sub_block_displ=NULL; 	//Displacements from starting point
	int *sub_block_len=NULL; 	//Lengths
	//int width_rectangle;
	int block_rows_begin;
	//int block_rows_displ;
	//int block_rows_end;
	//int height;
	//-- Check sizes and num blocks

	
	dtypes=mpi_datatypesCalloc(heights.size, "mpi rectangle array datatypes", comm, mpi_rank);
	#if DEBUG
	rangePrintf_rank_stage("triangArray", range, mpi_rank, stage);
	fflush(stdout);
	#endif
	
	#if DEBUG
	printf("[%d | %d] Heights %ld Dtypes %ld\n", 
		mpi_rank, stage, heights.size, dtypes.total
	);
	fflush(stdout);
	#endif 
	size_t total = dtypes.total;
	size_t deleted = 0;
	size_t count_add = 0;
	left_space=range.begin;
	//#pragma omp parallel for private(block)
	size_t sub_block_len_ = 0;
	for (block=0; block < total; block++){
		//block_rows_displ=displacements.buffer[block];
		block_rows_begin=count_add;//block_rows_displ / columns;
		count = heights.buffer[block-deleted];
		//count = heights.buffer[block];
		#if DEBUG
		printf("[%d | %d] Begin %d Count %d\n", mpi_rank, block, block_rows_begin, count);
		#endif 
		
		if (count > 0){
			reset_int_arr(&sub_block_len, count, "sub_block_len", comm, mpi_rank);
			reset_int_arr(&sub_block_displ, count, "sub_block_displ", comm, mpi_rank);
			//for (i = 0; i < count; i++){
			//	printf("len [%d] = %d\n", i, sub_block_len[i]);
			//	printf("disp [%d] = %d\n", i, sub_block_displ[i]);
			//}
			//-- Setup lenghts and displacement for each row
			//#pragma omp parallel for private(i)
			
			//up_space;
			dtypes.sizes[block-deleted]=0;
			//dtypes.sizes[block]=0;
			#if DEBUG
			printf(
				"[%d | %d  Range (b%ld,s%ld) ] TriangBlock %d | begin %d count %d || Left space %ld || size %d\n", 
				mpi_rank, stage, range.begin, range.size,
				block, block_rows_begin, count, left_space, (count*(count+1))/2
			);
			fflush(stdout);
			#endif 
			for ( i=0; i < count; i++ ) {
				up_space = (block_rows_begin+i)*columns;// + block_rows_begin;
				sub_block_len_=i+1+count_add;
				sub_block_len[ i ] = sub_block_len_;
				sub_block_displ[ i ] = up_space + left_space;
				dtypes.sizes[block-deleted]+=sub_block_len[i];
				//dtypes.sizes[block]+=sub_block_len[i];
				#if DEBUG > 1
				printf("[%d | %d | TriangBlock Range (%ld,%ld) %d | i %d ] rows before %d displ %d Up %ld LeftSpace %ld len %d\n", 
				mpi_rank, stage, range.begin, range.begin+range.size-1, block, i, block_rows_begin, sub_block_displ[i], up_space, left_space, sub_block_len[i]);
				fflush(stdout);
				#endif 
			}
			//left_space += count;
			#if DEBUG
			printf(
				"[%d | %d  Range (b%ld,s%ld) ] TriangBlock %d | count %d| First Len %d Last Len %d | First disp %d Last disp %d\n", 
				mpi_rank, stage, range.begin, range.size, block,count,
				sub_block_len[0],
				sub_block_len[count-1],
				sub_block_displ[0],
				sub_block_displ[count-1]
			);
			fflush(stdout);
			#endif 
			MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block-deleted) );
			MPI_CHECK( MPI_Type_commit( dtypes.arr + block - deleted) );
			//MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block) );
			//MPI_CHECK( MPI_Type_commit( dtypes.arr + block) );
			#if DEBUG > 1
				printf("[%d | %d | TriangBlock Range (%ld,%ld) %d ] Size %ld\n", 
					mpi_rank, stage, range.begin, range.begin+range.size-1, block, dtypes.sizes[block]);
				fflush(stdout);
			#endif 

			count_add += count;
		} else {
			dtypes.total -= 1;
			deleted += 1;
			dtypes.sizes[block]=0;
			#if DEBUG
			printf(
				"[%d | %d  Range (b%ld,s%ld) ] TriangBlock %d | Count 0\n", 
				mpi_rank, stage, range.begin, range.size, block 
				
				);
			#endif
		}
	}
	#if DEBUG
	for (block = 0; block < total; block++){
		printf("[%d | %d | TriangBlock Range Dtypes Check (%ld,%ld) %d ] Size %d\n", 
				mpi_rank, stage, range.begin, range.begin+range.size-1, block, dtypes.arr[block]);
				fflush(stdout);
	}
	#endif 
	return dtypes;
}


/*** 
 * @private create_type_combinedArray
 * * Type name: triang2type [trapezoid rectangle + triangle]
 * * Description
 * Combines rectangle and triangle datatypes
 * */
MPI_Datatypes create_type_combinedArray( 
	Range range, unsigned long int columns, 
	size_t block_limit_quotient_denominator, int mpi_rank, MPI_Comm comm, int stage,  
	Sizes heights, Sizes displacements
) {
	//-- Auxiliar
	size_t block;
	//Subtypes num of elements
	int lengths[2] = { 1, 1 }; 
	//Subtypes displacements
	MPI_Aint displ[2] = { 0, 0 };
	//Subtypes array
	MPI_Datatype types[2];
	//-- Datatypes
	MPI_Datatypes rectangles;
	MPI_Datatypes triangles;
	//Final array
	MPI_Datatypes dtypes = {0, NULL};
	//Rectangle and triangle array
	#if DEBUG
	printf("[%d | %d] Heights TypeArrayCombined %d\n", 
		mpi_rank, stage, heights.size
	);
	#endif 
	rectangles 	= create_type_rectangleArrayCombined(range, columns, mpi_rank, comm, heights, displacements);
	fflush(stdout);
	triangles 	= create_type_triangleArrayCombined(range, columns, mpi_rank, comm, heights, displacements, stage);
	dtypes.total=MAX(rectangles.total, triangles.total);
	dtypes = mpi_datatypesCalloc(dtypes.total, "datatypes", comm, mpi_rank);
	#if DEBUG
	printf("total: rectangles.size = %ld displacements.size = %ld triangles.size = %ld\n", rectangles.total, displacements.size, triangles.total);
	fflush(stdout);
	#endif
	//#pragma omp parallel for private(block) firstprivate(lengths, displ, types)
	for (block=0; block < dtypes.total; block++){
		#if DEBUG
		printf("[%d | %d] Building block %ld\n", mpi_rank, stage, block);
		fflush(stdout);
		printf(
			"[%d | %d] triangle %ld rectangle %ld\n",
			mpi_rank, stage, triangles.sizes[block], rectangles.sizes[block]
		);
		fflush(stdout);
		#endif 
		dtypes.sizes[block]=rectangles.sizes[block]+triangles.sizes[block];
		if (dtypes.sizes[block] > 0){
			#if DEBUG
			printf("[%d | %d | %ld] Block size > 0 \n", mpi_rank, stage, block);
			fflush(stdout);
			#endif 
			if ( rectangles.arr[block] == 0 ){
				
				#if DEBUG
				printf("[%d | %d | %ld] I am a triangle \n", mpi_rank, stage, block);
				fflush(stdout);
				#endif 
				dtypes.arr[block]	=	triangles.arr[block];
				#if DEBUG
				printf(
					"[%d | %d | %ld] triangle_size %ld dtype_size %ld Limit %ld\n",
					mpi_rank, stage, block, triangles.sizes[block], dtypes.sizes[block], BLOCK_LIMIT/block_limit_quotient_denominator
				);
				fflush(stdout);
				#endif 
			} else if (triangles.sizes[block] == 0){
				#if DEBUG
				printf("[%d | %d | %ld] I am a rectangle \n", mpi_rank, stage, block);
				fflush(stdout);
				#endif 
				dtypes.arr[block]=rectangles.arr[block];
				#if DEBUG
				printf(
					"[%d | %d | %ld] rectangle %ld dtype_size %ld Limit %ld\n",
					mpi_rank, stage, block, rectangles.sizes[block], dtypes.sizes[block], BLOCK_LIMIT/block_limit_quotient_denominator
				);
				fflush(stdout);
				#endif 
			} else {
				#if DEBUG
				printf("[%d | %d | %ld] I am a trapezoid \n", mpi_rank, stage, block);
				fflush(stdout);
				#endif 
				lengths[0] = 1;
				lengths[1] = 1;
				types[0]=rectangles.arr[block];
				types[1]=triangles.arr[block];
				//-- Setting up dtype (struct)
				MPI_CHECK( MPI_Type_create_struct( 2, lengths, displ, types, &dtypes.arr[block]) );				
				MPI_CHECK( MPI_Type_commit(&dtypes.arr[block]) );
				MPI_CHECK( MPI_Type_free(&types[0]));
				MPI_CHECK( MPI_Type_free(&types[1]));
				#if DEBUG
				printf(
					"[%d | %d | %ld] triangle %ld rectangle %ld dtype_size %ld Limit %ld\n",
					mpi_rank, stage, block,  triangles.sizes[block], rectangles.sizes[block], dtypes.sizes[block], BLOCK_LIMIT/block_limit_quotient_denominator
				);
				fflush(stdout);
				#endif 
			} 
		} else {
			#if DEBUG
			printf("WTF Empty block\n");
			fflush(stdout);
			#endif 
		}
	}
	freeNotNull(rectangles.arr);
	freeNotNull(triangles.arr);
	freeNotNull(heights.buffer);
	freeNotNull(displacements.buffer);
	return dtypes;
}

/***
 * @private  create_type_trapezoid
 * * Type name
 * Trapezoid
 * * Description 
 * All elements of the block until diagonal
 * **/
MPI_Datatype create_type_trapezoid( Range range, unsigned long int columns ) {
	//-- Datatype
	MPI_Datatype dtype;
	//-- Setting up indexed datatype
	int displ[ range.size ]; 
	int blocklen[ range.size ];
	int i;
	#pragma omp parallel for private(i)
	for ( i=0; i<range.size; i++ ) {
		blocklen[ i ] = range.begin + i + 1; 
		displ[ i ] = i * columns;
	}
	MPI_CHECK( 
		MPI_Type_indexed( 
			range.size, 
			blocklen, 
			displ, 
			MPI_DOUBLE, 
			&dtype
		) 
	);
	return dtype;
}
MPI_Datatypes create_type_trapezoidArray( 
	Range range, unsigned long int columns, int mpi_rank, MPI_Comm comm, Sizes heights, Sizes displacements
) {
	MPI_Datatypes dtypes = {0, NULL, NULL};
	//-- Auxiliar
	int i, block;
	//unsigned long int num_intermediate_blocks;
	//int block_height, block_max_size;
	//int last_block_height, last_block_rows_begin, last_block_rows_end;
	//-- Indexed type variables
	int count; 		//Number of rows in block
	int *sub_block_displ=NULL; 	//Displacements from starting point
	int *sub_block_len=NULL; 	//Lengths
	//int width_rectangle;
	int block_rows_begin;
	//int block_rows_displ;
	//int block_rows_end;
	//int height;
	//-- Check sizes and num blocks
	dtypes=mpi_datatypesCalloc(heights.size, "mpi rectangle array datatypes", comm, mpi_rank);
	size_t count_add = 0;
	size_t deleted=0;
	size_t width = range.begin;
	size_t total=dtypes.total;
	//#pragma omp parallel for private(block) -- Estos pragmas estan mal ... provocan errores
	for (block=0; block < total; block++){
		//block_rows_displ=displacements.buffer[block];
		//block_rows_begin=block_rows_displ/columns;
		block_rows_begin=count_add;
		//width = block_rows_begin+range.begin+(range.begin != 0);//block_rows_displ % columns + range.begin;
		count = heights.buffer[block];
		#if DEBUG
		printf("[%d | (%ld,%ld)]  TrapezoidBlock %d width: %ld | count %d\n", mpi_rank, range.begin, range.begin+range.size-1, block, width, count); 
		fflush(stdout);
		#endif 
		if (count > 0){
			reset_int_arr(&sub_block_len, count, "sub_block_len", comm, mpi_rank);
			reset_int_arr(&sub_block_displ, count, "sub_block_displ", comm, mpi_rank);
			//-- Setup lenghts and displacement for each row
			//#pragma omp parallel for private(i)
			for ( i=0; i < count; i++ ) {
				//sub_block_len[ i ] = range.begin + i+1;
				//sub_block_len[ i ] = count_add + range.begin + i+1;
				sub_block_len[ i ] = width + i + 1;
				sub_block_displ[ i ] = (block_rows_begin+i)*columns;//(count_add+i)*columns;
				dtypes.sizes[block]+=sub_block_len[i];
				#if DEBUG
				if (i == 0 || i == count-1){
				printf(
					"[%d (%ld,%ld)] Trapezoid %d | %d /%d | len %d disp %d | countadd %ld \n", 
					mpi_rank, 
					range.begin, range.begin+range.size-1, 
					block, i, count,
					sub_block_len[i], sub_block_displ[i],
					count_add
				);
				fflush(stdout);

				}
				#endif 
			}
			width+=count;
			MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block-deleted) );
			MPI_CHECK( MPI_Type_commit( dtypes.arr + block - deleted) );
			#if DEBUG
			printf(
				"[%d (%ld,%ld)] TrapezoidBlock %d | size %ld \n", 
				mpi_rank, range.begin, range.begin+range.size-1, 
				block,  dtypes.sizes[block]);
			fflush(stdout);
			#endif 
			count_add += count;
		} else {
			#if DEBUG
			printf("[%d (%ld,%ld)] %d dtypes.totalPrev=%ld\n", 
			mpi_rank, range.begin, range.begin+range.size-1, block,
			dtypes.total);
			#endif 
			//dtypes.sizes[block] = 0;
			dtypes.total -= 1;
			deleted += 1;
			#if DEBUG
			printf("[%d (%ld,%ld)] %d dtypes.total=%ld\n", 
			mpi_rank, range.begin, range.begin+range.size-1, block,
			dtypes.total);
			#endif 
		}
	}
	#if DEBUG
	printf("[%d (%ld,%ld)] dtypes.total=%ld\n", 
	mpi_rank, range.begin, range.begin+range.size-1,
	dtypes.total);
	fflush(stdout);
	#endif 
	return dtypes;
}

MPI_Datatypes create_type_trapezoidArray_( 
	Range range, unsigned long int columns, int mpi_rank, MPI_Comm comm, Sizes heights, Sizes displacements
) {
	MPI_Datatypes dtypes = {0, NULL, NULL};
	//-- Auxiliar
	int i, block;
	//unsigned long int num_intermediate_blocks;
	//int block_height, block_max_size;
	//int last_block_height, last_block_rows_begin, last_block_rows_end;
	//-- Indexed type variables
	int count; 		//Number of rows in block
	int *sub_block_displ=NULL; 	//Displacements from starting point
	int *sub_block_len=NULL; 	//Lengths
	//int width_rectangle;
	int block_rows_begin;
	int block_rows_displ;
	//int block_rows_end;
	//int height;
	//-- Check sizes and num blocks
	dtypes=mpi_datatypesCalloc(heights.size, "mpi rectangle array datatypes", comm, mpi_rank);
	long int count_add = 0;
	size_t width;
	//#pragma omp parallel for private(block) -- Estos pragmas estan mal ... provocan errores
	for (block=0; block < dtypes.total; block++){
		block_rows_displ=displacements.buffer[block];
		block_rows_begin=block_rows_displ/columns;
		width=block_rows_begin=block_rows_displ % columns;
		count = heights.buffer[block];
		
		if (count > 0){
			reset_int_arr(&sub_block_len, count, "sub_block_len", comm, mpi_rank);
			reset_int_arr(&sub_block_displ, count, "sub_block_displ", comm, mpi_rank);
			//-- Setup lenghts and displacement for each row
			//#pragma omp parallel for private(i)
			for ( i=0; i < count; i++ ) {
				//sub_block_len[ i ] = range.begin + i+1;
				//sub_block_len[ i ] = count_add + range.begin + i+1;
				sub_block_len[ i ] = width + i;
				sub_block_displ[ i ] = (count_add+i)*columns;
				dtypes.sizes[block]+=sub_block_len[i];
				#if DEBUG
				if (i == 0 || i == count-1){
				printf(
					"[%d (%ld,%ld)] Trapezoid %d | %d /%d | len %ld disp %d | countadd %ld \n", 
					mpi_rank, 
					range.begin, range.begin+range.size-1, 
					block, i, count,
					sub_block_len[i], sub_block_displ[i],
					count_add
				);
				fflush(stdout);

				}
				#endif 
			}
			MPI_CHECK( MPI_Type_indexed( count, sub_block_len, sub_block_displ, MPI_DOUBLE, dtypes.arr+block) );
			MPI_CHECK( MPI_Type_commit( dtypes.arr + block) );
			#if DEBUG
			printf(
				"[%d (%ld,%ld)] TrapezoidBlock %d | size %ld \n", 
				mpi_rank, range.begin, range.begin+range.size-1, 
				block,  dtypes.sizes[block]);
			fflush(stdout);
			#endif 
			count_add += count;
		} else {
			dtypes.sizes[block] = 0;
		}
	}

	return dtypes;
}
/***
 * @private  create_type_trapezoid_bilateral
 * * Type name
 * Trapezoid bilateral 
 * * Description 
 * All elements of the block until diagonal for multi-diagonal matrices
 * **/
MPI_Datatype create_type_trapezoid_bilateral( Range range, unsigned long int columns ) {
	MPI_Datatype result;
	int displ[ range.size ];
	int blocklen[ range.size ];
	int i;
	for ( i=0; i<range.size; i++ ) {
		blocklen[ i ] = range.begin + 1;
		displ[ i ] = i * (columns + 1);
	}
	MPI_CHECK( MPI_Type_indexed( range.size, blocklen, displ, MPI_DOUBLE, &result ) );
	return result;
}

/************************************
 * @section MAIN AUXILIAR FUNCTIONS	*
 * **********************************/
/***
 * @private show_usage
 * Util function: Show usage 
 * **/
void show_usage( 
	char *argv[], 
	int mpi_rank, 
	char *msg
) {
	if ( mpi_rank == 0 ) {
		fprintf( stderr, "\nError: %s\n", msg );
		fprintf( stderr, "\nUsage: %s <num_rows> <init_scheme> <output> <map_func> <comm_scheme> <comm_shape>\n\n", argv[0] );
		fprintf( stderr, "Init schemes:\n" );
		fprintf( stderr, "\tvalues\tA and B with consecutive fixed values\n" );
		fprintf( stderr, "\tb_id\tB identity, A consecutive values\n" );
		fprintf( stderr, "\ta_id\tA identity, B consecutive values\n" );
		fprintf( stderr, "\ta_id_check\tA identity, B consecutive values, automatically check that result is B\n" );
		fprintf( stderr, "\n" );
		fprintf( stderr, "Output C results:\n" );
		fprintf( stderr, "\toutput\tC is stored in output files\n" );
		fprintf( stderr, "\tno_output\tC is not stored in output files\n" );
		fprintf( stderr, "\n" );
		fprintf( stderr, "Map functions:\n" );
		fprintf( stderr, "\treg\tRegular with remainder in the last proc.\n" );
		fprintf( stderr, "\treg_double\tRegular, but first takes double, remainder in the last proc.\n" );
		fprintf( stderr, "\tbalanced\tBalanced according to the total number of non-zero elements\n" );
		fprintf( stderr, "\n" );
		fprintf( stderr, "Skewed communication schemes:\n" );
		fprintf( stderr, "\towner_skew\tData of A is sent always from owner in skewed mode\n" );
		fprintf( stderr, "\towner_broad\tData of A is sent always from owner using broadcast\n" );
		fprintf( stderr, "\tpipeline\tData of A is moved across consecutive procs in a pipeline\n" );
		fprintf( stderr, "\n" );
		fprintf( stderr, "Communication shapes:\n" );
		fprintf( stderr, "\tfull\tFull rectangular shapes with all the columns\n" );
		fprintf( stderr, "\tboxes\tMinimal rectangular shapes, bounding boxes including the triangle\n" );
		fprintf( stderr, "\ttriang\tExact comm. One comm. with a single derived data type\n" );
		fprintf( stderr, "\ttriang_2buff\tExact comm. Two comms. Rectangular shape with a derived data type plus triangular shape with manual marshalling in a buffer\n" );
		fprintf( stderr, "\ttriang_2type\tExact comm. Two comms. Rectangular shape with a derived data type plus triangular shape with a derived data type\n" );
		fprintf( stderr, "\ttrapezoid\tExact comm. One indexed type with the trapezoid (rectangle + triangle\n" );
		fprintf( stderr, "\n" );
	}
	//MPI_Abort( MPI_COMM_WORLD, INVALID_PROCS_NUM );
	error_mpi_abort( MPI_COMM_WORLD, INVALID_ARGS, mpi_rank, msg);
}

/***
 * @private set_ending_times
 * Util function: set times by the end of execution.
 * **/
void set_ending_times(int mpi_rank,double total_t,double precomp_t,double comm1_t,double comm2_t,double mult_t,double for_t){
	double total_t_red 	= 0.0;
	double prec_t_red 	= 0.0;
	double comm1_t_red 	= 0.0; 
	double mult_t_red	= 0.0;
	double comm2_t_red	= 0.0;
	double comm_acumm_manticore=0.0;
	double comm_acumm_gorgon=0.0;
	double comm_max_manticore=0.0;
	double comm_max_gorgon=0.0;
	double comm_manticore=0.0;
	double comm_gorgon=0.0;
	double total_max_manticore=0.0;
	double total_max_gorgon=0.0;
	double total_acumm_manticore=0.0;
	double total_acumm_gorgon=0.0;
	double total_manticore = 0.0;
	double total_gorgon = 0.0;
	
	char *processor;
	processor=calloc(MPI_MAX_PROCESSOR_NAME, sizeof(char));
	int proc_len=0;
	MPI_Get_processor_name(processor, &proc_len);
	if ( ! strcmp( processor, "manticore" ) ){
		comm_manticore=comm1_t+comm2_t;
		total_manticore=total_t;
	} else {
		comm_gorgon=comm1_t+comm2_t;
		total_gorgon=total_t;
	}
	printf("Server %s\n", processor);
	//printf("[%d | %s ] Mult Time Local %.3lf Comm Time Local %.3lf\n", mpi_rank, processor, mult_t, comm1_t);
	#ifdef __clang_major__
	//printf("[%d | %s] Compilado con CLANG\n", mpi_rank, processor);
	fflush(stdout);
	#else 
	//printf("[%d | %s] Compilado con otra cosa\n", mpi_rank, processor);
	fflush(stdout);
	#endif 
	fflush(stdout);
	double for_t_red	= 0.0;
	MPI_CHECK( 
		MPI_Reduce( 
			&total_gorgon,
			&total_acumm_gorgon,
			1, 
			MPI_DOUBLE, 
			MPI_SUM,
			0, MPI_COMM_WORLD 
		) 
	);
	MPI_CHECK( 
		MPI_Reduce( 
			&comm_gorgon,
			&comm_acumm_gorgon,
			1, 
			MPI_DOUBLE, 
			MPI_SUM,
			0, MPI_COMM_WORLD 
		) 
	);
	MPI_CHECK( 
		MPI_Reduce( 
			&total_manticore,
			&total_acumm_manticore,
			1, 
			MPI_DOUBLE, 
			MPI_SUM,
			0, MPI_COMM_WORLD 
		) 
	);
	MPI_CHECK( 
		MPI_Reduce( 
			&comm_manticore,
			&comm_acumm_manticore,
			1, 
			MPI_DOUBLE, 
			MPI_SUM,
			0, MPI_COMM_WORLD 
		) 
	);
	MPI_CHECK( 
		MPI_Reduce( 
			&total_gorgon,
			&total_max_gorgon,
			1, 
			MPI_DOUBLE, 
			MPI_MAX,
			0, MPI_COMM_WORLD 
		) 
	);
	MPI_CHECK( 
		MPI_Reduce( 
			&total_manticore,
			&total_max_manticore,
			1, 
			MPI_DOUBLE, 
			MPI_MAX,
			0, MPI_COMM_WORLD 
		) 
	);
	MPI_CHECK( 
		MPI_Reduce( 
			&comm_gorgon,
			&comm_max_gorgon,
			1, 
			MPI_DOUBLE, 
			MPI_MAX,
			0, MPI_COMM_WORLD 
		) 
	);
	MPI_CHECK( 
		MPI_Reduce( 
			&comm_manticore,
			&comm_max_manticore,
			1, 
			MPI_DOUBLE, 
			MPI_MAX,
			0, MPI_COMM_WORLD 
		) 
	);

	MPI_CHECK( 
		MPI_Reduce( 
			&total_t,
			&total_t_red,
			1, 
			MPI_DOUBLE, 
			MPI_MAX,
			0, MPI_COMM_WORLD 
		) 
	);
	MPI_CHECK( MPI_Reduce( 
		&precomp_t, 
		&prec_t_red,      
		1, 
		MPI_DOUBLE, 
		MPI_MAX, 
		0, 
		MPI_COMM_WORLD 
		) 
	);
	MPI_CHECK( MPI_Reduce( &comm1_t,	&comm1_t_red,	1,	MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD	) );
	MPI_CHECK( MPI_Reduce( &mult_t,		&mult_t_red,    1,	MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD	) );
	MPI_CHECK( MPI_Reduce( &comm2_t,	&comm2_t_red,   1,	MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD	) );
	MPI_CHECK( MPI_Reduce( &for_t,		&for_t_red,     1,	MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD	) );
	//printf( "\n[%d]LTotal time (red_max): %lf\n", mpi_rank, total_t );
	if ( mpi_rank == 0 ) {
		int sizeeeee;
		MPI_Comm_size(MPI_COMM_WORLD, &sizeeeee);
		//printf( "\nTotal time (red_max): %lf - total mpis %d\n", total_t_red, sizeeeee);
		printf( "\nMult t: %lf\nTotal t: %lf\ntotal mpis %d\n", mult_t_red, total_t_red, sizeeeee);
/*
		printf( "Precomp time (red_max): %lf\n", prec_t_red );
		printf( "Comm1 time (red_max): %lf\n", comm1_t_red );
		printf( "Mult time (red_max): %lf\n", mult_t_red );
		printf( "Comm2 time (red_max): %lf\n", comm2_t_red );
		printf( "For loop time (red_max): %lf\n", for_t_red );
		printf( "SUM_Partials (red_max): %lf\n", prec_t_red+comm1_t_red+mult_t_red+comm2_t_red );
		printf( "Total time - SUM_Partials (red_max): %lf\n", total_t_red - (prec_t_red+comm1_t_red+mult_t_red+comm2_t_red) );
		printf( "Manticore comm time (red_max): %lf\n", comm_max_manticore);
		printf( "Gorgon comm time (red_max): %lf\n", comm_max_gorgon);
		printf( "Manticore total time (red_max): %lf\n", total_max_manticore);
		printf( "Gorgon total time (red_max): %lf\n", total_max_gorgon);
		printf( "Total - comunicaciones (red_max): %lf", total_t_red - comm_max_gorgon - comm_max_manticore);
*/
	}
	MPI_CHECK( MPI_Reduce( &total_t,       &total_t_red,     1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
	MPI_CHECK( MPI_Reduce( &precomp_t,     &prec_t_red,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
	MPI_CHECK( MPI_Reduce( &comm1_t,       &comm1_t_red,     1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
	MPI_CHECK( MPI_Reduce( &mult_t,        &mult_t_red,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
	MPI_CHECK( MPI_Reduce( &comm2_t,       &comm2_t_red,     1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
	MPI_CHECK( MPI_Reduce( &for_t,         &for_t_red,       1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) );
	if ( mpi_rank == 0 ) {
/*
		printf( "\nTotal time (red_sum): %lf\n", total_t_red );
		printf( "Precomp time (red_sum): %lf\n", prec_t_red );
		printf( "Comm1 time (red_sum): %lf\n", comm1_t_red );
		printf( "Mult time (red_sum): %lf\n", mult_t_red );
		printf( "Comm2 time (red_sum): %lf\n", comm2_t_red );
		printf( "For loop time (red_sum): %lf\n", for_t_red );
		printf( "SUM_Partials (red_sum): %lf\n", prec_t_red+comm1_t_red+mult_t_red+comm2_t_red );
		printf( "Total time - SUM_Partials (red_sum): %lf\n", total_t_red - (prec_t_red+comm1_t_red+mult_t_red+comm2_t_red) );
		printf( "Manticore comm time (red_sum): %lf\n", comm_acumm_manticore);
		printf( "Gorgon comm time (red_sum): %lf\n", comm_acumm_gorgon);
		printf( "Manticore total time (red_sum): %lf\n", total_acumm_manticore);
		printf( "Gorgon total time (red_sum): %lf\n", total_acumm_gorgon);
		printf( "Total - comunicaciones (red_sum): %lf ", total_t_red - comm_acumm_gorgon - comm_acumm_manticore);
*/
	}
}

/***
 * @private set_initial_scheme
 * Util function: sets up Params_init_scheme
 * **/
void set_initial_scheme(
	int argc, char *argv[],
	int mpi_rank,
	Params_init *init_scheme,
	int *do_check
){
	if ( ! strcmp( argv[2], "values" ) ){
		*init_scheme = INIT_VALUES;
	} else if (! strcmp( argv[2], "values_check" )){
		*init_scheme = INIT_VALUES;
		*do_check = 1;
	} else if ( ! strcmp( argv[2], "a_id" ) || !strcmp(argv[2], "a_id_triag_triag") ){
		*init_scheme = INIT_A_ID;
	} else if ( ! strcmp( argv[2], "a_id_check") || !strcmp(argv[2], "a_id_triag_triag_check") )  {
		*init_scheme = INIT_A_ID;
		*do_check = 1;
	} else if ( ! strcmp( argv[2], "b_id") || !strcmp(argv[2], "b_id_triag_triag")){
		*init_scheme = INIT_B_ID;

	} else show_usage( argv, mpi_rank, "Unknown initialization scheme");
}

/***
 * @private set_output_result
 * Util function: sets up params_output
 * **/
void set_output_result(
	int argc, char *argv[],
	int mpi_rank,
	Params_output *output
){
	if ( ! strcmp( argv[3], "output" ) ){
		*output = OUTPUT_YES;
		#if DEBUG
		printf("Output will be generated\n");
		#endif 
	}
	else if ( ! strcmp( argv[3], "no_output" ) )
		*output = OUTPUT_NO;
	else show_usage( argv, mpi_rank, "Unknown output selector");
	
}

/***
 * @private set_mapping_function
 * Util function: sets up Map_function
 * **/
void set_mapping_function(
	int argc, 
	char *argv[],
	int mpi_rank,
	Map_function *f_mapping
){
	if ( ! strcmp( argv[4], "reg" ) )
		*f_mapping = map_regular;
	else if ( ! strcmp( argv[4], "reg_double" ) )
		*f_mapping = map_regular_first_double;
	else if ( ! strcmp( argv[4], "balanced" ) )
		*f_mapping = map_balanced;
	else show_usage( argv, mpi_rank, "Unknown name of mapping function" );
}

/***
 * @private set_communication_scheme
 * Util function: sets up Params_comm_scheme and Params_comm_func
 * **/
void set_communication_scheme(
	int arg, char *argv[],
	int mpi_rank,
	Params_comm_scheme *comm_scheme,
	Params_comm_func *comm_func
){
	*comm_func = FUNC_P2P;
	if ( ! strcmp( argv[5], "owner_skew" ) ){
		*comm_scheme = COMM_OWNER;
		fflush(stdout);
		//printf("Owner skew not yet implemented\n");
		//MPI_Abort(EXIT_FAILURE, MPI_COMM_WORLD);
	} else if ( ! strcmp( argv[5], "owner_broad" ) ) {
		*comm_scheme = COMM_OWNER;
		*comm_func = FUNC_BROADCAST;
	} 
	else if ( ! strcmp( argv[5], "owner_broad_sarteco" ) ) {
		*comm_scheme = COMM_OWNER;
		*comm_func = FUNC_BROADCAST_SARTECO;
	} 
	else if ( ! strcmp( argv[5], "pipeline" ) )
		*comm_scheme = COMM_PREVIOUS;
	else show_usage( argv, mpi_rank, "Unknown communication scheme");
}

/***
 * @private set_communication_scheme
 * Util function: sets up Params_comm_shape & Params_comm_mechanism
 * **/
void set_communication_shape(
	int arg, char *argv[],
	int mpi_rank,
	Params_comm_shape *comm_shape,
	Params_comm_mechanism *comm_mechanism,
	Type_function *create_type
){
	if ( ! strcmp( argv[6], "full" ) ){
		#if DEBUG
		printf("COMM SHAPE FULL\n");
		#endif 
		*comm_shape = SHAPE_FULL;
	}
	else if ( ! strcmp( argv[6], "boxes" ) ){
		*comm_shape = SHAPE_BOXES;
		#if DEBUG
		printf("COMM SHAPE BOXES\n");
		#endif 
		*create_type = create_type_bbox;
	}
	else if ( ! strcmp( argv[6], "triang" ) ){
		*comm_shape = SHAPE_TRIANG_1;
		#if DEBUG
		printf("COMM SHAPE TRIANG1\n");
		#endif 
		*create_type = create_type_combined;
	}
	else if ( ! strcmp( argv[6], "triang_2buff" ) ) {
		*comm_shape = SHAPE_TRIANG_2;
		*comm_mechanism = MECH_BUFFER;
		*create_type = create_type_rectangle;
		#if DEBUG
		printf("COMM SHAPE TRIANG2BUFF\n");
		#endif 
	}
	else if ( ! strcmp( argv[6], "triang_2type" ) ) {
		*comm_shape = SHAPE_TRIANG_2;
		*comm_mechanism = MECH_TYPE;
		*create_type = create_type_bbox;
		#if DEBUG
		printf("COMM SHAPE TRIANG2TYPE\n");
		#endif 
	}
	else if ( ! strcmp( argv[6], "trapezoid" ) ){
		*comm_shape = SHAPE_TRAPEZOID;
		*create_type = create_type_trapezoid;
		#if DEBUG
		printf("COMM SHAPE TRAPEZOID\n");
		#endif 
	} 
	else if (! strcmp( argv[6], "trapezoid_bilateral" ) ){
		*comm_shape = SHAPE_TRAPEZOID_BILATERAL;
	}
	else show_usage( argv, mpi_rank, "Unknown name of communication shape");
	#if DEBUG
	printf("Comm shape %d\n", *comm_shape);
	#endif 
}
/***
 * @private set_type_send
 * Sets up type_send1 and type_send2 depending on selected comm_shape and comm_mechanism
 * **/
Buff set_type_send(
	int mpi_rank, 
	MPI_Datatype *type_send1, 
	MPI_Datatype *type_send2,
	Buff *triang_buff_send, 
	MatrixRange *mat_A,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism,
	Params_comm_func comm_fun
){
	//unsigned long int block_size;
	switch (comm_shape){
		case SHAPE_FULL:
		//block_size = mat_A->rows.size * mat_A->columns.size;
		#if DEBUG
			printf("[%d] No sending derived datatype is needed | rank \n", mpi_rank);
			#endif 
		break;
		case SHAPE_BOXES: // Rectangle
			*type_send1 = create_type_bbox( mat_A->rows, mat_A->columns.size );
			MPI_CHECK( MPI_Type_commit( type_send1 ) );
		break;
		case SHAPE_TRIANG_1: // Combined Rectangle + Triangle
			
			//block_size = mat_A->rows.size *(mat_A->rows.begin + mat_A->rows.size);
			*type_send1 = create_type_combined( mat_A->rows, mat_A -> columns.size);
				#if DEBUG
			rangePrintf_rank_stage("rectRg", mat_A -> rows, mpi_rank, -1);
			printf("[%d] count = %ld\n", mpi_rank, mat_A -> rows.size);
			printf("[%d] len = %ld\n",  mpi_rank, mat_A -> rows.begin+1);
			printf("[%d] stride = %ld\n",mpi_rank, mat_A -> columns.size);
			#endif 
			MPI_CHECK( MPI_Type_commit( type_send1 ) );
		break;
		case SHAPE_TRIANG_2: // Rectangle + Triangle
			error_mpi_abort(MPI_COMM_WORLD, INVALID_COMMUNICATION_SHAPE, mpi_rank, "Invalid communication shape");
			//block_size = mat_A->rows.size *(mat_A->rows.begin + mat_A->rows.size);
			*type_send1 = create_type_rectangle(mat_A -> rows, mat_A -> columns.size );
			*type_send1 = create_type_bbox(mat_A -> rows, mat_A -> columns.size );
			MPI_CHECK( MPI_Type_commit( type_send1 ) );
			if ( comm_mechanism == MECH_TYPE ) { // Type
				*type_send2 = create_type_triangle(mat_A -> rows, mat_A -> columns.size );
				MPI_CHECK( MPI_Type_commit( type_send2 ) );
			} else {// Precompute buffer with data
				*triang_buff_send = set_triangular_buffer_data(mat_A->rows, mat_A->columns.size, mat_A -> mat, mpi_rank);
				MPI_CHECK( MPI_Type_commit( type_send1 ) );
			}
		break;
		case SHAPE_TRAPEZOID: // Trapezoid
			error_mpi_abort(MPI_COMM_WORLD, INVALID_COMMUNICATION_SHAPE, mpi_rank, "Invalid communication shape");
			//block_size = mat_A->rows.size *(mat_A->rows.begin + mat_A->rows.size);
			*type_send1 = create_type_trapezoid( mat_A->rows, mat_A -> columns.size);
			MPI_CHECK( MPI_Type_commit( type_send1 ) );
		break;
		case SHAPE_TRAPEZOID_BILATERAL:
			error_mpi_abort(MPI_COMM_WORLD, INVALID_COMMUNICATION_SHAPE, mpi_rank, "Invalid communication shape");
			//block_size = mat_A->rows.size *(mat_A->rows.begin + mat_A->rows.size)/2;
			*type_send1 = create_type_trapezoid_bilateral( mat_A->rows, mat_A -> columns.size);
			MPI_CHECK( MPI_Type_commit( type_send1 ) );
		break;
		default:
			printf("Invalid communication shape\n");
			printf("Communication shape: %d\n", comm_shape);
			printf("Communication shape = %d | %d | %d | %d | %d", SHAPE_FULL, SHAPE_BOXES, SHAPE_TRIANG_1, SHAPE_TRIANG_2, SHAPE_TRAPEZOID);
			error_mpi_abort(MPI_COMM_WORLD, INVALID_COMMUNICATION_SHAPE, mpi_rank, "Invalid communication shape");
		break;
	}
	#if DEBUG
	printf("[%d] Sending derived types setted up \n", mpi_rank);
	#endif 
	fflush(stdout);
	return *triang_buff_send;
}


/***
 * @private set_type_send
 * Sets up type_send1 and type_send2 depending on selected comm_shape and comm_mechanism
 * **/
Buff set_type_sendArray(
	int mpi_rank, 
	MPI_Datatypes *type_send1, 
	MPI_Datatypes *type_send2,
	Buff *triang_buff_send, 
	MatrixRange *mat_A,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism,
	Params_comm_func comm_func,
	size_t block_limit_quotient_denominator, 
	MPI_Comm comm, int stage
){
	#if DEBUG
	printf("[%d | (%ld, %ld) set type send Array\n", mpi_rank, mat_A->rows.begin, mat_A->rows.size)	;
	fflush(stdout);
	#endif 
	//fflush(stdout);
	//unsigned long int block_size;
	Sizes heights={ NULL, 0 }, displacements = { NULL, 0 };
	#if DEBUG
	printf("[%d] set type send \n",mpi_rank);
	#endif
	switch (comm_shape){
		case SHAPE_FULL:
		//block_size = mat_A->rows.size * mat_A->columns.size;
		#if DEBUG
			printf("[%d] No sending derived datatype is needed | rank \n", mpi_rank);
			#endif 
		break;
		case SHAPE_BOXES: // Rectangle
			#if DEBUG
			printf("You have selected BOXES\n");
			fflush(stdout);
			#endif 
			*type_send1 = create_type_bboxArray(mat_A -> rows, mat_A -> columns.size, block_limit_quotient_denominator, mpi_rank, comm);
			//*type_send1 = create_type_bbox( mat_A->rows, mat_A->columns.size );
			//MPI_CHECK( MPI_Type_commit( type_send1 ) );
		break;
		case SHAPE_TRIANG_1: // Combined Rectangle + Triangle
			#if DEBUG
			printf("Estoy aqui shape triang1 send\n");
			printf("mpi rank %d stage %d block limit quotient denominator %ld\n", mpi_rank, stage, block_limit_quotient_denominator);
			#endif
			split_trapezoid(mat_A->rows, mat_A->columns.size, block_limit_quotient_denominator, mpi_rank, comm, &heights, &displacements);
			#if DEBUG
			printf("[%d | %d] Heights After split trapezoid %ld \n", 
				mpi_rank, stage, heights.size
			);
			fflush(stdout);
			#endif 
			//printf("After split trapezoid\n");
			//fflush(stdout);
			//printf("Create type combined Array\n");
			//fflush(stdout);
			*type_send1 = create_type_combinedArray(mat_A -> rows, mat_A -> columns.size, block_limit_quotient_denominator, mpi_rank, comm, stage, heights, displacements);
			//printf("After combined array\n");
			//fflush(stdout);
			//block_size = mat_A->rows.size *(mat_A->rows.begin + mat_A->rows.size);
			//*type_send1 = create_type_combined( mat_A->rows, mat_A -> columns.size);
			//MPI_CHECK( MPI_Type_commit( type_send1 ) );
		break;
		case SHAPE_TRIANG_2: // Rectangle + Triangle
			#if DEBUG
			printf("SHAPE TRIANG 2\n");
			fflush(stdout);
			#endif 
			//split_trapezoid(mat_A->rows, mat_A->columns.size, block_limit_quotient_denominator, mpi_rank, comm, &heights, &displacements); --> COMPROBAR QUE HAY QUE AÑADIRLO
			*type_send1 = create_type_rectangleArray(mat_A -> rows, mat_A -> columns.size, mpi_rank, comm, heights, displacements);
			//*type_send1 = create_type_bbox(mat_A -> rows, mat_A -> columns.size );
			//MPI_CHECK( MPI_Type_commit( type_send1 ) );
			if ( comm_mechanism == MECH_TYPE ) { // Type
				*type_send2 = create_type_triangleArray(mat_A -> rows, mat_A -> columns.size, mpi_rank, comm, heights, displacements, stage);
				#if DEBUG
				printf("[%d] datatypes total %ld  || first %ld\n", mpi_rank, type_send2->total, type_send2->total);
				fflush(stdout);
				#endif 
			} else {// Precompute buffer with data
				*triang_buff_send = set_triangular_buffer_data(mat_A->rows, mat_A->columns.size, mat_A -> mat, mpi_rank);
			}
		break;
		case SHAPE_TRAPEZOID: // Trapezoid
			//error_mpi_abort(MPI_COMM_WORLD, INVALID_COMMUNICATION_SHAPE, mpi_rank, "Invalid communication shape");
			//block_size = mat_A->rows.size *(mat_A->rows.begin + mat_A->rows.size);
			split_trapezoid(mat_A->rows, mat_A->columns.size, block_limit_quotient_denominator, mpi_rank, comm, &heights, &displacements);
			#if DEBUG
			printf("[%d | %d (b%ld, s%ld)] Heights After split trapezoid %ld \n", 
				mpi_rank, stage, mat_A->rows.begin, mat_A->rows.size, heights.size
			);
			fflush(stdout);
			#endif 
			*type_send1 = create_type_trapezoidArray( mat_A->rows, mat_A -> columns.size, mpi_rank, comm, heights, displacements);
			//MPI_CHECK( MPI_Type_commit( type_send1 ) );
		break;
		case SHAPE_TRAPEZOID_BILATERAL:
			error_mpi_abort(MPI_COMM_WORLD, INVALID_COMMUNICATION_SHAPE, mpi_rank, "Invalid communication shape");
			//block_size = mat_A->rows.size *(mat_A->rows.begin + mat_A->rows.size)/2;
			//*type_send1 = create_type_trapezoid_bilateral( mat_A->rows, mat_A -> columns.size);
			//MPI_CHECK( MPI_Type_commit( type_send1 ) );
		break;
		default:
			printf("Invalid communication shape\n");
			printf("Communication shape: %d\n", comm_shape);
			printf("Communication shape = %d | %d | %d | %d | %d", SHAPE_FULL, SHAPE_BOXES, SHAPE_TRIANG_1, SHAPE_TRIANG_2, SHAPE_TRAPEZOID);
			error_mpi_abort(MPI_COMM_WORLD, INVALID_COMMUNICATION_SHAPE, mpi_rank, "Invalid communication shape");
		break;
	}
	#if DEBUG
	printf("[%d] Sending derived types setted up \n", mpi_rank);
	#endif 
	#if DEBUG
	if (comm_shape != SHAPE_FULL) printf("[%d] set type send : %ld\n",mpi_rank, type_send1 -> total);
	#endif
	fflush(stdout);
	return *triang_buff_send;
}



/***
 * @private set_type_recv
 * Sets up type_recv1 and type_recv2 depending on selected comm_shape and comm_mechanism
 * */
void set_type_recv(
	int mpi_rank, 
	Range remote_rows, 
	int columns,
	MPI_Datatype *type_recv1,
	MPI_Datatype *type_recv2,
	Buff *triang_buff_recv,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism
){
	size_t size;
	switch (comm_shape){
		case SHAPE_FULL:
			#if DEBUG
			printf("No receiving derived datatype is needed | rank %d\n", mpi_rank);
			#endif 
		break;
		case SHAPE_BOXES: // Rectangle
			#if DEBUG
			printf("[%d] Bbox ", mpi_rank);
			rangePrintf("remote_rows", remote_rows);
			printf("[%d] columns = %d\n", mpi_rank, columns);
			#endif 
			*type_recv1 = create_type_bbox( remote_rows, columns);
			MPI_CHECK( MPI_Type_commit( type_recv1) );
		break;
		case SHAPE_TRIANG_1: // Combined Rectangle + Triangle
			*type_recv1 = create_type_combined( remote_rows, columns );
			MPI_CHECK( MPI_Type_commit( type_recv1) );
		break;
		case SHAPE_TRIANG_2: // Rectangle + Triangle
			*type_recv1 = create_type_rectangle( remote_rows, columns );
			MPI_CHECK( MPI_Type_commit( type_recv1) );
			if ( comm_mechanism == MECH_TYPE ) { // Type
				*type_recv2 = create_type_triangle( remote_rows, columns );
				MPI_CHECK( MPI_Type_commit( type_recv2) );
			} else {//Buffer
				size = ( remote_rows.size + 1 ) * remote_rows.size / 2;
				*triang_buff_recv = buff(size, mpi_rank);
			}
		break;
		case SHAPE_TRAPEZOID: // Trapezoid
			*type_recv1 = create_type_trapezoid(remote_rows, columns );
			MPI_CHECK( MPI_Type_commit( type_recv1) );
		break;
		case SHAPE_TRAPEZOID_BILATERAL:
			*type_recv1 = create_type_trapezoid_bilateral( remote_rows, columns);
			MPI_CHECK( MPI_Type_commit( type_recv1) );
		break;
		default:
			printf("Invalid communication shape\n");
			printf("Communication shape: %d\n", comm_shape);
			printf("Communication shape = %d | %d | %d | %d | %d", SHAPE_FULL, SHAPE_BOXES, SHAPE_TRIANG_1, SHAPE_TRIANG_2, SHAPE_TRAPEZOID);
			//MPI_Abort(EXIT_FAILURE, INVALID_COMMUNICATION_SHAPE);
			error_mpi_abort(MPI_COMM_WORLD, INVALID_COMMUNICATION_SHAPE, mpi_rank, "Invalid communication shape");
		break;
	}
	#if DEBUG
	printf("Receiving derived types setted up | %d\n", mpi_rank);
	#endif 
	fflush(stdout);
}



/***
 * @private set_type_recv
 * Sets up type_recv1 and type_recv2 depending on selected comm_shape and comm_mechanism
 * */
void set_type_recvArray(
	int mpi_rank, 
	Range remote_rows, 
	int columns,
	MPI_Datatypes *type_recv1,
	MPI_Datatypes *type_recv2,
	Buff *triang_buff_recv,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism,
	size_t block_limit_quotient_denominator, 
	MPI_Comm comm, int stage
){
	#if DEBUG
	printf("[%d] set type recv : \n",mpi_rank);
	#endif
	size_t size;
	Sizes heights = {NULL, 0};
	Sizes displacements = {NULL, 0};
	switch (comm_shape){
		case SHAPE_FULL:
			#if DEBUG
			printf("No receiving derived datatype is needed | rank %d\n", mpi_rank);
			#endif 
		break;
		case SHAPE_BOXES: // Rectangle
			#if DEBUG
			printf("[%d] Bbox ", mpi_rank);
			rangePrintf("remote_rows", remote_rows);
			printf("[%d] columns = %d\n", mpi_rank, columns);
			#endif 
			*type_recv1 = create_type_bboxArray( remote_rows, columns, block_limit_quotient_denominator, mpi_rank, MPI_COMM_WORLD);
			//MPI_CHECK( MPI_Type_commit( type_recv1) );
		break;
		case SHAPE_TRIANG_1: // Combined Rectangle + Triangle
			#if DEBUG
			printf("Estoy aqui shape triang1 recv\n");
			#endif 
			split_trapezoid(remote_rows, columns, block_limit_quotient_denominator, mpi_rank, comm, &heights, &displacements);
			*type_recv1 = create_type_combinedArray( remote_rows, columns, block_limit_quotient_denominator, mpi_rank, comm, stage, heights, displacements);
		break;
		case SHAPE_TRIANG_2: // Rectangle + Triangle
			split_trapezoid(remote_rows, columns, block_limit_quotient_denominator, mpi_rank, comm, &heights, &displacements);
			*type_recv1 = create_type_rectangleArray(remote_rows, columns, mpi_rank, comm, heights, displacements);
			//*type_send1 = create_type_bbox(mat_A -> rows, mat_A -> columns.size );
			//MPI_CHECK( MPI_Type_commit( type_send1 ) );
			if ( comm_mechanism == MECH_TYPE ) { // Type
				*type_recv2 = create_type_triangleArray(remote_rows, columns, mpi_rank, comm, heights, displacements, stage);
				#if DEBUG
				printf("[%d] datatypes total %ld  || first %ld\n", mpi_rank, type_recv2->total, type_recv1->total);
				fflush(stdout);
				#endif 
			} else {// Precompute buffer with data
				size = ( remote_rows.size + 1 ) * remote_rows.size / 2;
				*triang_buff_recv = buff(size, mpi_rank);
			}
		break;
		case SHAPE_TRAPEZOID: // Trapezoid
		//	*type_recv1 = create_type_trapezoid(remote_rows, columns );
		//	MPI_CHECK( MPI_Type_commit( type_recv1) );
			split_trapezoid(remote_rows, columns, block_limit_quotient_denominator, mpi_rank, comm, &heights, &displacements);
			#if DEBUG
			printf("Recv mpi rank %d stage %d block limit %ld\n", mpi_rank, stage, block_limit_quotient_denominator);
			fflush(stdout);
			#endif 
			*type_recv1 = create_type_trapezoidArray( remote_rows, columns, mpi_rank, comm, heights, displacements);
		break;
		//case SHAPE_TRAPEZOID_BILATERAL:
		//	*type_recv1 = create_type_trapezoid_bilateral( remote_rows, columns);
		//	MPI_CHECK( MPI_Type_commit( type_recv1) );
		//break;
		default:
			printf("Invalid communication shape\n");
			printf("Communication shape: %d\n", comm_shape);
			printf("Communication shape = %d | %d | %d | %d | %d", SHAPE_FULL, SHAPE_BOXES, SHAPE_TRIANG_1, SHAPE_TRIANG_2, SHAPE_TRAPEZOID);
			//MPI_Abort(EXIT_FAILURE, INVALID_COMMUNICATION_SHAPE);
			error_mpi_abort(MPI_COMM_WORLD, INVALID_COMMUNICATION_SHAPE, mpi_rank, "Invalid communication shape");
		break;
	}
	#if DEBUG
	printf("Receiving derived types setted up | %d\n", mpi_rank);
	if (comm_shape != SHAPE_FULL) printf("[%d] set type recv : %ld\n",mpi_rank, type_recv1 -> total);
	fflush(stdout);
	#endif
}

/********************************************
 * @section MATRICES GENERATION				*
 * Functions to generate data matrices		*
 * ******************************************/
/**
 * @private mat_identity			
 * Generates identity matrix
 * **/
void mat_identity(Range rows, Range columns, double *mat){
	int i, j;
	for ( i=0; i < rows.size; i++ ) {
		for ( j=0; j < columns.size; j++ ) {
			if ( i + rows.begin == j ) 
				mat[ i * columns.size + j ] = 1;
			else
				mat[ i * columns.size + j ] = 0;
		}
	}
}

/***
 * @private mat_consecutive_values
 * Generates matrices by consecutive values
 * **/
void mat_consecutive_values(Range rows, Range columns, double *mat, int triangular, int global_columns){
	int i, j;
	int elem;
	int block_before = rows.begin*columns.size; //Añadido por Macu por problemas observados en la generacion de A
	for ( i=0; i< rows.size; i++ ) {
		for ( j=0; j< columns.size; j++ ) {
			elem = i*columns.size+j;
			// Triangular matrix
			if (triangular && (j > rows.begin + i) ) {
				mat[elem] = 0;
			} else {
				mat[elem] = ( 
					( i*global_columns+j+columns.begin ) + block_before
				) % VALUES_LIMIT + 1;
				#if MININUMS
				mat[elem] = 1/mat[elem];
				#endif 
				/*mat[elem] = ( 
					( i*global_columns+j+columns.begin )
				) % VALUES_LIMIT + 1;*/
			}
		}
	}
}
/**
 * @private matrix_rangeConsecutive_values
 * Fills up mat by using consecutive values.
 * **/
void matrix_rangeConsecutive_values(MatrixRange *mat, int global_columns){
	int triangular = (mat -> type == TRIANGULAR_D);
	mat_consecutive_values(mat -> rows, mat -> columns, mat -> mat, triangular, global_columns);
}

/***
 * @private matrix_rangePoint_to
 * Changes sizes of Matrix to Remote sizes and changes array pointer to Remote's pointer
 * */
void matrix_rangePoint_to(MatrixRange *mat, MatrixRange remote){
	mat -> rows = remote.rows;
	mat -> columns = remote.columns;
	mat -> mat = remote.mat;
	mat -> type = remote.type;
}

/***
 * @private 
 * Used for setting up mat_A data
 * **/
void matrix_rangeSetA(
	MatrixRange *mat, 
	Params_init init_scheme, 
	int global_columns
){
	#if DEBUG
	printf("Setting mat A\n");
	#endif 
	if ( init_scheme == INIT_A_ID ) {
		#if DEBUG
		printf("Set A ID\n");
		#endif 
		mat_identity(mat -> rows, mat -> columns, mat -> mat);
	} else {
		#if DEBUG
		printf("Set A consecutive\n"); 
		#endif 
		matrix_rangeConsecutive_values(mat, global_columns);
	}
}


/***
 * @private 
 * Used for setting up mat_B data
 * **/
void matrix_rangeSetB(
	MatrixRange *mat, 
	Params_init init_scheme, 
	int global_columns
){
	if ( init_scheme == INIT_B_ID ) {
		mat_identity(mat -> rows, mat -> columns, mat -> mat);
	} else {
		matrix_rangeConsecutive_values(mat, global_columns);
	}
}

/****************************************************************
 * @section CHECK RESULTS										*
 * This part contains functions for checking obtained results.	*
 * **************************************************************/
/***
 * @private check_results_aid
 * Function for checking results when A is identity
 * **/
void check_results_aid(
	Range local_cols, int rows, 
	double *mat_B, double *mat_C,
	int *check_result
){
	//printf("[%d|%d] Check result insi0 %d\n", -1, -1, *check_result);
	int i, j;
	int row, column;
	int flatten;
	for ( i=0; i<rows; i++ ) {
		row = i;
		for ( j=0; j<local_cols.size; j++ ) {
			column = local_cols.begin + j;
			flatten = i*local_cols.size+j;
			if ( mat_C[ flatten ] != mat_B[ flatten ] ) {
				*check_result = 0;
				printf(
					"C[%d][%d] != B[%d][%d] || %.3lf || %.3lf\n",
					row, column, row, column,
					mat_C[ flatten ],
					mat_B[ flatten ]
				);
			}
		}
	}
	//printf("[%d|%d] Check result insi1 %d\n", -1, -1, *check_result);
}
/***
 * @private matrix_product
 * Returns A*B
 * //-- Corregir: no utiliza mat_C para nada
 * **/
double *matrix_product(
	Range local_cols, 
	Range current_range,
	int columns, 
	double *mat_A, 
	double *mat_B, 
	double *mat_C,
	int *check_result,
	int rows,
	int mpi_rank
){
	#if DEBUG
	printf("--> matrix product\n");
	printf("current_range =(%ld,%ld)", current_range.begin, current_range.size);
	printf("local cols = (%ld, %ld)\n", local_cols.begin, local_cols.size);
	fflush(stdout);
	#endif 
	int aux_size = current_range.size * local_cols.size;
	#if DEBUG
	printf("aux size = %d\n", aux_size);
	fflush(stdout);
	#endif
	double *aux = calloc_abort( aux_size, sizeof(double), "aux", MPI_COMM_WORLD, mpi_rank);
	double A = 0.0, B = 0.0;
	int ra, ca, cb;
	for ( ra = 0; ra < current_range.size; ra++ ) {
		for ( ca=0; ca <= columns && ca < rows; ca++) {//i+current_range.begin; k++ ) {
			for ( cb=0; cb < local_cols.size; cb++ ) {
				if ( isnan(mat_A[ ra * columns + ca ]) ) {
					mat_A[ ra * columns + ca ] = 0.0;
				}
				A = mat_A[ ra * columns + ca ];
				B = mat_B[ ca * local_cols.size + cb];
				aux[ra * local_cols.size + cb] +=  A*B;
			}
		}
	}
	#if DEBUG
	printf("--> matrix product\n");
	fflush(stdout);
	#endif 
	return aux;
}


/*** --> En stand by porque igual no es necesario
 * @private matrix_product
 * Returns A*B
 * //-- Corregir: no utiliza mat_C para nada
 * **/
/*double *matrix_rangeProduct(
	MatrixRange mat_A, 
	MatrixRange mat_B, 
	matrixRange *mat_C,
	int *check_result,
){
	unsigned long int aux_size = mat_A.rows.size.size * mat_B.rows.size;
	unsigned long int ra, ca, cb;
	unsigned long int ca_limit;
	double *aux = calloc_abort( aux_size, sizeof(double), "aux", MPI_COMM_WORLD);
	double A = 0.0, B = 0.0;
	if (mat_A.type == RECTANGULAR){
		ca_limit = rows
	}
	for ( ra = 0; ra < mat_A.rows.size; ra++ ) {
		for ( ca=0; ca <= columns && ca < rows; ca++) {//i+current_range.begin; k++ ) {
			for ( cb=0; cb < local_cols.size; cb++ ) {
				if ( isnan(mat_A[ ra * columns + ca ]) ) {
					mat_A[ ra * columns + ca ] = 0.0;
				}
				A = mat_A[ ra * columns + ca ];
				B = mat_B[ ca * local_cols.size + cb];
				aux[ra * local_cols.size + cb] +=  A*B;
			}
		}
	}
	#if DEBUG
	printf("--> matrix product\n");
	fflush(stdout);
	#endif 
	return aux;
}*/


/***
 * @private matrix_equal
 * Checks if two matrices are equal.
 * If true, it returns 1.
 * Otherwise, it returns 0.
 * **/
_Bool matrix_equal(
	Range columns, 
	Range rows,
	double *mat_A, 
	double *mat_B,
	int block_before_A,
	int block_before_B
){
	// Check
	double diff = 0.0;
	double diff_ = 0.0;
	double A, B;
	_Bool eq = 1;
	int row, column, flatten; //Auxiliar variables
	int ra, ca;
	for (ra=0; ra < rows.size; ra++ ) {
		row = ra + rows.begin;
		for (ca=0; ca < columns.size; ca++ ) {
			column = ca + columns.begin;
			flatten = ra*columns.size+ca;
			A = mat_A[flatten+block_before_A];
			B = mat_B[flatten+block_before_B];
			diff_ = fabs( A-B );
			diff += diff_;
			if (diff_ > 0){
				printf(
					"FAILED c[%d][%d] = %.3lf != %.3lf\n", 
					row, column, B, A
				);
			}
		}
	}
	if(diff>0)	{
		printf("PRODUCT ERROR %f\n", diff);
		eq = 0;
	}
	return eq;
}

/***
 * @private check_results_init_values
 * Function for checking results when A is 
 * generated by using consecutive values
 * //-- Mejorar: se puede simplificar utilizando el struct MatrixRange
 * **/
void check_results_init_values(
	Range local_cols, 
	Range current_range,
	int columns, 
	double *mat_A,  
	double *mat_B, 
	double *mat_C,
	int *check_result,
	int rows,
	int mpi_rank
){
	#if DEBUG
	printf("--> Check results init values\n");
	#endif 
	int block_before = current_range.begin*local_cols.size; //Initial position C
	#if DEBUG
	printf_matrix_ranges_rank_stage(mat_C, current_range.begin, current_range.begin+current_range.size, local_cols.begin, local_cols.size+local_cols.begin, "c_check", -1, -2,block_before, 0);
	#endif 
	//-- Calculate product
	double *aux = matrix_product(local_cols, current_range, columns, mat_A, mat_B, mat_C, check_result, rows, mpi_rank);
	#if DEBUG
	printf_matrix_ranges_rank_stage(aux, current_range.begin, current_range.begin+current_range.size, local_cols.begin, local_cols.size+local_cols.begin, "aux", -1, -2,0,0);
	fflush(stdout);
	#endif 
	//-- Check result
	*check_result = matrix_equal(local_cols, current_range, aux, mat_C, 0, block_before);
	free(aux);
	#if DEBUG
	printf("Check results init values -->\n");
	fflush(stdout);
	#endif 
}

/***
 * @private free_types_recv
 * Frees reciving types and buffers and 
 * restores current_A if needed.
 * **/
void free_types_recv(
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism,
	MPI_Datatype *type_recv1, 
	MPI_Datatype *type_recv2, 
	MatrixRange *current_A,
	Buff *triang_buffer_recv
){
	switch (comm_shape)
	{
	case SHAPE_FULL:
		#if DEBUG
		printf("No type needs freedom\n");
		#endif 
		break;
	case SHAPE_TRIANG_2:
		if ( comm_mechanism == MECH_TYPE ) {
			MPI_Type_free( type_recv1 );
			MPI_Type_free( type_recv2 ); //Added
		}
		// Triangle part: buffer unmarshalling 
		if ( comm_shape == SHAPE_TRIANG_2 && comm_mechanism == MECH_BUFFER ) {
			unsigned long int i,j, pos = 0;
			#pragma omp parallel for private(i, j) firstprivate(pos)
			for ( i=0; i<current_A->rows.size; i++ ) {
				for ( j=current_A->rows.begin; j<=i+current_A->rows.begin; j++ ) {
					current_A->mat[ i*current_A -> columns.size+j ] = triang_buffer_recv->buffer[ pos ];
					pos ++;
				}
			}
			// Free triang buffers
			free( triang_buffer_recv->buffer );
		}
		break;
	default:
		MPI_Type_free( type_recv1);
		break;
	}
}


/***
 * @private free_types_recvArray
 * Frees reciving types and buffers and 
 * restores current_A if needed.
 * **/
void free_types_recvArray(
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism,
	MPI_Datatypes *type_recv1, 
	MPI_Datatypes *type_recv2, 
	MatrixRange *current_A,
	Buff *triang_buffer_recv
){
	size_t typ;
	switch (comm_shape)
	{
	case SHAPE_FULL:
		#if DEBUG
		printf("No type needs freedom\n");
		#endif 
		break;
	case SHAPE_TRIANG_2:
		if ( comm_mechanism == MECH_TYPE ) {
			for  (typ = 0; typ < type_recv2->total; typ++){
				if (type_recv2->sizes[typ] > 0){
					MPI_Type_free( &type_recv2->arr[typ] );
				}
			}
			free (type_recv2->arr);
			free (type_recv2->sizes);
		}
		// Triangle part: buffer unmarshalling 
		if ( comm_shape == SHAPE_TRIANG_2 && comm_mechanism == MECH_BUFFER ) {
			unsigned long int i,j, pos = 0;
			#pragma omp parallel for private(i, j) firstprivate(pos)
			for ( i=0; i<current_A->rows.size; i++ ) {
				for ( j=current_A->rows.begin; j<=i+current_A->rows.begin; j++ ) {
					current_A->mat[ i*current_A -> columns.size+j ] = triang_buffer_recv->buffer[ pos ];
					pos ++;
				}
			}
			// Free triang buffers
			free( triang_buffer_recv->buffer );
		}
		//break; --> Asi tengo lo del type1 solo 1 vez
	default:
		for  (typ = 0; typ < type_recv1->total; typ++){
			if (type_recv1->sizes[typ] > 0){
				MPI_Type_free( &type_recv1->arr[typ] );
			}
		}
		free (type_recv1->arr);
		free (type_recv1->sizes);
		break;
	}
}
/****************************************
 * @section OUTPUT RESULTS FUNCTIONS	*
 * **************************************/
/***
 * @private print_output_results
 * Prints mat_C into output files
 * //-- Mejora: Utilizar struct MatrixRange
 * **/
void print_output_results(
	int rows, int columns,
	Range local_cols, 
	int mpi_rank,
	double *mat_C, 
	int mpi_procs,
	char *func
){
	// Open rank file for writting
	char file_name[51];
	int i;
	sprintf( file_name, "./tests/outputs/output_%s.%d.%d.%d.dat", func,rows, mpi_procs, mpi_rank );
	//sprintf( file_name, "%soutput.%d.dat", file_path, mpi_rank );
	FILE *file = fopen( file_name, "wb" );
	if ( file == NULL ) {
		fprintf( stderr, "[print_output_results] Error: Opening the file %s\n", file_name );
		//MPI_Abort( MPI_COMM_WORLD, ERROR_OPENING_FILE );
		error_mpi_abort( MPI_COMM_WORLD, ERROR_OPENING_FILE, mpi_rank, "Error opening file");
	}
	// Write header: Mat sizes and local ranges of the output
	fwrite( &rows, sizeof(int), 1, file );
	fwrite( &columns, sizeof(int), 1, file );
	int output_range[4] = { 0, rows, local_cols.begin, local_cols.size };
	fwrite( output_range, sizeof(int), 4, file );
	// Write data
	for ( i=0; i<rows; i++ ) {
		fwrite( &mat_C[ i*local_cols.size ], sizeof(double), local_cols.size, file );
	}
	fclose( file );
}

void print_output_results_generic(
	int rows, int columns,
	Range local_rows,
	Range local_cols, 
	int mpi_rank,
	double *mat_C, 
	int mpi_procs,
	char *func
){
	// Open rank file for writting
	char file_name[51];
	int i;
	sprintf( file_name, "./tests/outputs/output_%s.%d.%d.%d.dat", func,rows, mpi_procs, mpi_rank );
	//sprintf( file_name, "%soutput.%d.dat", file_path, mpi_rank );
	FILE *file = fopen( file_name, "wb" );
	if ( file == NULL ) {
		fprintf( stderr, "Error: Opening the file %s\n", file_name );
		error_mpi_abort( MPI_COMM_WORLD, ERROR_OPENING_FILE, mpi_rank, "Error opening file");
	}
	// Write header: Mat sizes and local ranges of the output
	fwrite( &rows, sizeof(int), 1, file );
	fwrite( &columns, sizeof(int), 1, file );
	int output_range[4] = { local_rows.begin, local_rows.size, local_cols.begin, local_cols.size };
	#if DEBUG
	printf(
		"Printing C -> {%d|%d|%d|%d} in file %s\n", 
		output_range[0],
		output_range[1],
		output_range[2],
		output_range[3], file_name
	);
	fflush(stdout);
	#endif 
	fwrite( output_range, sizeof(int), 4, file );
	// Write data
	for ( i=0; i<local_rows.size; i++ ) {
		fwrite( &mat_C[ i*local_cols.size ], sizeof(double), local_cols.size, file );
	}
	fclose( file );
}

/***
 * @private matrix_traspose
 * Trasposes mat_A into mat_A_trasp
 * //-- Traducir a matrixRange 
 * //-- ¿Merecería la pena introducir un .block_before a matrixRange?
 * **/
void matrix_traspose(
	double *mat_A, 
	double *mat_A_trasp,
	int rows_A,
	int columns_A,
	int block_before_A, 
	int block_before_A_tras
){	
	int r, c, flatten, flatten_tras;
	for (r=0; r< rows_A; r++) {
		for(c=0; c< columns_A; c++) {
			flatten = block_before_A+(r*columns_A+c);
			flatten_tras = block_before_A_tras+(c*rows_A+r);
			mat_A_trasp[flatten_tras]=mat_A[flatten];
		}
	}
}

/***
 * @section Range functions
 * Printing function
 * **/
int rangePrintf(const char *name, Range range){
	return printf("%s.(s,b) = (%ld,%ld)\n", name, range.size, range.begin);
}
/***
 * @private rangePrintfRank_Stage
 * Printing function
 * **/
int rangePrintf_rank_stage(const char *name, Range range, int mpi_rank, int stage){
	return printf("[%d | %d] %s.(s,b) = (%ld,%ld)\n", mpi_rank, stage, name, range.size, range.begin);
}

/***
 * @private rangesCalloc_abort
 * Builder
 * **/

Ranges rangesCalloc_abort(
	unsigned int total, 
	char *name, 
	MPI_Comm comm,
	int mpi_rank
){
	Ranges ranges;
	ranges.array = calloc_abort(total, sizeof(Range), name, comm, mpi_rank);
	ranges.total = total;
	return ranges;
}

/***
 * @private rangesPrint_sizes
 * Printing function for sizes attributes
 * **/
/*int rangesCalloc_abort(Ranges ranges, char *name, MPI_Comm comm){
	int total_print_result = 0;
	unsigned int i = 0;
	printf("[");
	for (i = 0; i < ranges.total; i++){
		printf("En proceso\n");
	}
	return total_print_result;
}
*/


/**
 * @private mpi_broadcast_splitted_quotient
 * @brief MPI_Ibcast adapt to bigger count 
 */
void mpi_broadcast_splitted_quotient(
    double *buffer,
    unsigned long int count,
    int root, 
    MPI_Comm comm, 
    Requests *requests,//Esto lo pondre como MPI_Requests en el "oficial"
    int mpi_rank,
	MPI_Datatype datatype,
	unsigned int quotient
){
    unsigned long int blocks_total;
	unsigned int block_size_max=BLOCK_LIMIT;
	if (quotient == 0){
		blocks_total=0;
	} else {
		block_size_max=MAX(0,BLOCK_LIMIT/quotient);
		if (block_size_max == 0){
			blocks_total = 0;
		} else {
			blocks_total = MAX(0,count/block_size_max);
		}
		//if (blocks_total > MAX_BLOCKS_TOTAL){
		//	blocks_total = MAX_BLOCKS_TOTAL;
		//	block_size_max = count/MAX_BLOCKS_TOTAL;
		//}
		#if DEBUG
		printf("[%d] Quotient: %u Block_max_size: %u | Total %lu | count %ld\n", mpi_rank, quotient, block_size_max, blocks_total, count);
		fflush(stdout);
		#endif 
	}
	#if DEBUG
	printf("[%d | %d] BLOCKS: %ld | %d size\n", mpi_rank, root, blocks_total, block_size_max);
	#endif 
    unsigned long int block_before = 0;
    unsigned long int block_id = 0;
    unsigned int block_size_last;
	if(blocks_total == 0) {
		block_size_last = count;
		#if DEBUG
		printf("count %ld sizelast %d\n", count, block_size_last);
		#endif 
	} else {
		if (quotient == 0){
			block_size_last=count;
		} else {
			block_size_last=count-blocks_total*block_size_max;
		}
		#if DEBUG
		printf("Total %ld | count %ld sizelast %d\n", blocks_total, count, block_size_last);
		#endif 
	}
    if (requests -> arr != NULL){
		free(requests -> arr);
		requests->arr=NULL;
	}
	if (block_size_last > 0){
		requests -> size = blocks_total + 1;
	} else {
		requests -> size = blocks_total;
	}
	#if DEBUG
	printf("[%d] Requestss size %d\n", mpi_rank, requests -> size); 
	fflush(stdout);
	#endif 
	*requests = requestsCalloc_abort(requests -> size, comm, mpi_rank);
    //#pragma omp parallel for private(block_id)
	for (block_id = 0; block_id < blocks_total; block_id++){
        block_before = block_id*block_size_max;
        MPI_CHECK(
			MPI_Ibcast( 
				buffer + block_before, 
				block_size_max, 
				datatype, 
				root, 
				comm, 
				&requests -> arr[block_id]
			) 
		) ;
		#if DEBUG
		printf("[%d | %d] Bcast Buffer block_id %ld, size %d\n", mpi_rank, root, block_id, block_size_max);
		fflush(stdout);
		#endif 
    }
    if(block_size_last > 0){
		if (quotient == 0){
			block_before = 0;
		} else  {
			block_before = blocks_total*block_size_max;
		}
		block_id=requests -> size-1;
        MPI_CHECK(
			MPI_Ibcast( 
				buffer + block_before, 
				block_size_last, 
				datatype, 
				root, 
				comm,
				&requests -> arr[block_id] 
			) 
		);
		#if DEBUG
		printf("[%d | %d] Bcast Buffer block_id %ld, size (last) %d\n", mpi_rank, root, block_id, block_size_last);
		fflush(stdout);
		#endif 
    }
}




/**
 * @private mpi_broadcast_splitted_quotient
 * @brief MPI_Ibcast adapt to bigger count 
 */
void mpi_broadcast_splitted_quotientArray(
    double *buffer,
    int root, 
    MPI_Comm comm, 
    Requests *requests,//Esto lo pondre como MPI_Requests en el "oficial"
    int mpi_rank,
	MPI_Datatypes datatypes
){
	#if DEBUG
	printf("[%d | %d] Dtypes %ld Requests %d before\n", 
	mpi_rank, root, datatypes.total, requests->size);
	#endif 
	requests -> size = datatypes.total;
	freeNotNull(requests->arr);
	requests->arr=NULL;
	*requests = requestsCalloc_abort(requests -> size, comm, mpi_rank);
	#if DEBUG
	printf("[%d | %d] Dtypes %ld Requests %d after\n", 
		mpi_rank, root, datatypes.total, requests->size
	);
	printf("[%d | %d] Quotient Array requests_size: %d\n", mpi_rank, root, requests->size);
	fflush(stdout);
	#endif 
	size_t block_id, current_req=0;
    //#pragma omp parallel for 
	for (block_id = 0; block_id < datatypes.total; block_id++){
		if (datatypes.sizes[block_id] > 0){
		#if DEBUG
		printf(
			"[%d | %d | %ld] Broadcast QA || Requests %d types %ld\n",
			mpi_rank, root,  block_id,
			requests->size,
			datatypes.total
		);
		fflush(stdout);
		printf(
			"[%d | %d | %ld] Broadcast QA || Current req %ld Request %d Ad %ls\n",
			mpi_rank, root,  block_id,
			current_req, requests->arr[current_req], &requests->arr[current_req]
		);
		fflush(stdout);
		printf(
			"[%d | %d | %ld] Broadcast QA || Type %d Ad %p\n",
			mpi_rank, root,  block_id,
			datatypes.arr[block_id], (void*) &datatypes.arr[block_id]
		);
		fflush(stdout);
		printf(
			"[%d | %d | %ld] Broadcast QA || Buffer %p Comm %p\n",
			mpi_rank, root,  block_id,
			(void*) buffer, (void*) &comm
		);
		fflush(stdout);
		#endif
        	MPI_CHECK(
				MPI_Ibcast( 
					buffer,
					1, 
					datatypes.arr[block_id], 
					root, 
					comm, 
					&requests -> arr[current_req]
				) 
			) ;
			current_req++;
		} else {
			#if DEBUG > 1
			printf("[%d | %d] No data to send\n", mpi_rank, root);
			fflush(stdout);
			#endif 
			requests->size=MAX(0,requests->size-1);
		}
    }
}


/**
 * @private mpi_broadcast_splitted
 * @brief MPI_Ibcast adapt to bigger count 
 */
void mpi_broadcast_splitted(
    double *buffer,
    unsigned long int count,
    int root, 
    MPI_Comm comm, 
    Requests *requests,//Esto lo pondre como MPI_Requests en el "oficial"
    int mpi_rank,
	MPI_Datatype datatype
){
	int total_coms;
	MPI_Comm_size( comm, &total_coms );//Total num of procs
	unsigned int block_size_max=BLOCK_LIMIT/total_coms;
    unsigned long int blocks_total = MAX(0,count/block_size_max);
	if (block_size_max < 1 || blocks_total > BLOCK_LIMIT){
		error_mpi_abort(MPI_COMM_WORLD, INVALID_PROCS_NUM, mpi_rank, "too many procs for broadcast mode for this size");
	}
	#if DEBUG
	printf("[%d] BLOCKS: %ld\n", mpi_rank, blocks_total);
	#endif 
    unsigned long int block_before = 0;
    unsigned long int block_id = 0;
    unsigned int block_size_last = count - blocks_total*block_size_max;
    if (requests -> arr != NULL){
		free(requests -> arr);
		requests->arr=NULL;
	}
	if (block_size_last > 0){
		requests -> size = blocks_total + 1;
	} else {
		requests -> size = blocks_total;
	}
	*requests = requestsCalloc_abort(requests -> size, comm, mpi_rank);
    //#pragma omp parallel for 
	for (block_id = 0; block_id < blocks_total; block_id++){
        block_before = block_id*block_size_max;
        MPI_CHECK(
			MPI_Ibcast( 
				buffer + block_before, 
				block_size_max, 
				datatype, 
				root, 
				comm, 
				&requests -> arr[block_id]
			) 
		) ;
		#if DEBUG
		printf("[%d | %ld] Send %d block_id %d\n", root, block_id, mpi_rank, root);
		fflush(stdout);
		#endif 
    }
    if(block_size_last > 0){
		block_before = blocks_total*block_size_max;
        MPI_CHECK(
			MPI_Ibcast( 
				buffer + block_before, 
				block_size_last, 
				datatype, 
				root, 
				comm,
				&requests -> arr[requests -> size -1] 
			) 
		);
		#if DEBUG
		printf("[%d | %ld] Send %d block_id %ld\n", root, block_id, mpi_rank, blocks_total+1);
		fflush(stdout);
		#endif 
    }
}

/**
 * @private mpi_wait_splitted
 * @brief MPI wait for splitted comms (adapt to bigger count)
 */
void mpi_wait_splitted(Requests *requests, Statuses *status){
	unsigned long int block_id;
	if (status->size < requests->size){
		status->size= requests->size;
		status->arr	= (MPI_Status*) realloc(status->arr, sizeof(MPI_Status)*status->size);
	}
	for (block_id = 0; block_id < requests -> size; block_id++){
		MPI_CHECK(
			MPI_Wait(
				&requests -> arr[block_id], 
				&status -> arr[block_id]
			)
		);
	} 
}



/**
 * @private MPI_Isend_splitted
 * @brief MPI_Isend adapt to bigger count 
 * 
 */
void mpi_isend_splitted(
    double *buffer,
    unsigned long int count,
    int comm_to, 
    MPI_Comm comm, 
    Requests *requests,//Esto lo pondre como MPI_Requests en el "oficial"
    int mpi_rank,
	MPI_Datatype datatype,
	int request_id,
	int max_requests_id,
	int tag, 
	int init
){
	tag*=mpi_rank+1;
	int block_tag, block_req_id;
	unsigned long int blocks_total = MAX(0,count/BLOCK_LIMIT);
    unsigned long int block_before = 0;
    unsigned long int block_id = 0;
    unsigned int block_size_last = count - blocks_total*BLOCK_LIMIT;
	requests->size=(blocks_total+(block_size_last != 0))*max_requests_id;
	int total_coms;
	MPI_Comm_size( comm, &total_coms );//Total num of procs
    if (requests -> arr != NULL){
		free(requests -> arr);
    	requests-> arr = NULL;
	}
	if (init) {
		*requests = requestsCalloc_abort(requests -> size, comm, mpi_rank);
	} else {
		*requests = requestsMalloc_abort(requests -> size, comm, mpi_rank);
	}
	#if DEBUG
		printf("[%d | %d] Send Intermediate Blocks [%ld | %d] \n", mpi_rank, comm_to, blocks_total, requests->size);
		fflush(stdout);
	#endif 
    for (block_id = 0; block_id < blocks_total; block_id++){
		block_req_id=block_id*max_requests_id+request_id;
		block_tag=tag+block_req_id;
		MPI_CHECK(MPI_Isend(buffer + block_before, BLOCK_LIMIT, datatype, comm_to, block_tag, comm, &requests -> arr[block_req_id]) );
        block_before += BLOCK_LIMIT;
		#if DEBUG
			printf("[%d | %d | %ld | %ld | %d] req_send_block %d tag %d\n", mpi_rank, comm_to, block_id, blocks_total, block_req_id, block_size_last, block_tag);
			fflush(stdout);
		#endif 
    }
    if(block_size_last > 0){
		#if DEBUG
			printf("[%d | %d] Send Last Blocks\n", mpi_rank, comm_to);
			fflush(stdout);
		#endif 
		block_req_id=block_id*max_requests_id+request_id;
		block_tag=tag+block_req_id;
		MPI_CHECK(MPI_Isend(buffer + block_before, block_size_last, datatype, comm_to, block_tag, comm, &requests -> arr[block_req_id]) );
		#if DEBUG
		printf("[%d | %d | %ld | %ld | %d] req_send_block %d tag %d\n", mpi_rank, comm_to, block_id, blocks_total, block_req_id, block_size_last, block_tag);
		fflush(stdout);
		#endif 
    }
	#if DEBUG 
	printf("[%d -> %d] Isend %d last %d\n", mpi_rank, comm_to, requests->size, block_size_last );
	fflush(stdout);
	#endif 
}


/**
 * @private mpi_Irecv_splitted
 * @brief MPI_Irecv adapt to bigger count 
 * 
 */
void mpi_irecv_splitted(
    double *buffer,
    unsigned long int count,
    int comm_from, 
    MPI_Comm comm, 
    Requests *requests,//Esto lo pondre como MPI_Requests en el "oficial"
    int mpi_rank,
	MPI_Datatype datatype,
	int request_id,
	int max_requests_id,
	int tag
){
	#if DEBUG
	printf("[%d | %d] --> IRECV SPLITTED \n", mpi_rank, comm_from); fflush(stdout);
	#endif 
	int block_tag, block_req_id;
	unsigned long int blocks_total = MAX(0,count/ BLOCK_LIMIT);
    unsigned long int block_before = 0;
    unsigned long int block_id = 0;
    unsigned int block_size_last = count - (blocks_total* BLOCK_LIMIT);
	tag*=comm_from+1;
	#if DEBUG
	printf("[%d | %d] requests -> size %d before\n", mpi_rank, comm_from, requests->size);
	printf("[%d | %d] requests -> size %ld after \n", mpi_rank, comm_from, (blocks_total+(block_size_last != 0))*max_requests_id);
	fflush(stdout);
	#endif 
	requests -> size =(blocks_total+(block_size_last != 0))*max_requests_id;
	#if DEBUG
	printf("[%d | %d] requests -> size %d\n", mpi_rank, comm_from, requests->size);
	#endif 
	int total_coms;
	MPI_Comm_size( comm, &total_coms );//Total num of procs
    free(requests -> arr);
    requests-> arr = NULL;
	*requests = requestsCalloc_abort(requests -> size, comm, mpi_rank);
	fflush(stdout);
    for (block_id = 0; block_id < blocks_total; block_id++){
		block_req_id=block_id*max_requests_id+request_id;
		block_tag=tag+block_req_id;
		MPI_CHECK(
			MPI_Irecv(
				buffer + block_before, BLOCK_LIMIT, 
				datatype, comm_from, block_tag, comm, 
				&requests -> arr[block_req_id]
			) 
		);
        block_before += BLOCK_LIMIT;
		#if DEBUG
		printf("[%d | %d | %ld | %ld | %d] req_recvd_block %d tag %d\n", mpi_rank, comm_from, block_id, blocks_total, block_req_id, block_size_last, block_tag);
		fflush(stdout);
		#endif 
    }
    if(block_size_last > 0){
		block_req_id=block_id*max_requests_id+request_id;
		block_tag=tag+block_req_id;
		MPI_CHECK(
			MPI_Irecv(
				buffer + block_before, 
				block_size_last, datatype, 
				comm_from, block_tag, comm, 
				&requests -> arr[block_req_id]
			) 
		);
		#if DEBUG
		printf("[%d | %d | %ld | %ld | %ld] req_recvd_block %d tag %d\n", mpi_rank, comm_from, block_id, blocks_total, block_id*max_requests_id+request_id, block_size_last, block_tag);
		fflush(stdout);
		#endif 
    }
	#if DEBUG 
	printf("[%d <- %d] IRecv %d last %d | bt %ld\n", mpi_rank, comm_from, requests->size, block_size_last, blocks_total);
	fflush(stdout);
	printf("[%d | %d] IRECV SPLITTED --> \n", mpi_rank, comm_from); fflush(stdout);
	#endif 
}

void mpi_isend_array(
	double *buffer,
    int comm_to, 
    MPI_Comm comm, 
    Requests *requests,//Esto lo pondre como MPI_Requests en el "oficial"
    int mpi_rank,
	MPI_Datatypes datatypes,
	int request_id,
	int max_requests_id,
	int tag, 
	int init
){
	tag*=mpi_rank+1;
	int block_tag, block_req_id;
	unsigned long int block_id = 0;
    requests->size=datatypes.total;
	int total_coms;
	MPI_Comm_size( comm, &total_coms );//Total num of procs
    if (requests -> arr != NULL){
		free(requests -> arr);
    	requests-> arr = NULL;
	}
	if (init) {
		*requests = requestsCalloc_abort(requests -> size, comm, mpi_rank);
	} else {
		*requests = requestsMalloc_abort(requests -> size, comm, mpi_rank);
	}
    for (block_id = 0; block_id < datatypes.total; block_id++){
		if (datatypes.sizes[block_id] > 0){
			block_req_id=block_id*max_requests_id+request_id;
			block_tag=tag+block_req_id;
			MPI_CHECK(MPI_Isend(buffer, 1, datatypes.arr[block_id], comm_to, block_tag, comm, &requests -> arr[block_req_id]) );
		}
    }
    	
}


/**
 * @private mpi_Irecv_array
 * @brief MPI_Irecv adapt to bigger count 
 * 
 */
void mpi_irecv_array(
    double *buffer,
    int comm_from, 
    MPI_Comm comm, 
    Requests *requests,//Esto lo pondre como MPI_Requests en el "oficial"
    int mpi_rank,
	MPI_Datatypes datatypes,
	int request_id,
	int max_requests_id,
	int tag
){
	#if DEBUG
	printf("[%d | %d] --> IRECV SPLITTED \n", mpi_rank, comm_from); fflush(stdout);
	#endif 
	int block_tag, block_req_id;
	unsigned long int block_id = 0;
    tag*=comm_from+1;
	#if DEBUG
	printf("[%d | %d] requests -> size %d before\n", mpi_rank, comm_from, requests->size);
	fflush(stdout);
	#endif 
	requests -> size = datatypes.total;
	#if DEBUG
	printf("[%d | %d] requests -> size %d\n", mpi_rank, comm_from, requests->size);
	#endif 
	int total_coms;
	MPI_Comm_size( comm, &total_coms );//Total num of procs
    free(requests -> arr);
    requests-> arr = NULL;
	*requests = requestsCalloc_abort(requests -> size, comm, mpi_rank);
	fflush(stdout);
    for (block_id = 0; block_id < datatypes.total; block_id++){
		if (datatypes.sizes[block_id] > 0){
			block_req_id=block_id*max_requests_id+request_id;
			block_tag=tag+block_req_id;
			MPI_CHECK(
				MPI_Irecv(
					buffer, 1, 
					datatypes.arr[block_id], comm_from, block_tag, comm, 
					&requests -> arr[block_req_id]
				) 
			);
		}
    }
	#if DEBUG 
	printf("[%d | %d] IRECV SPLITTED --> \n", mpi_rank, comm_from); fflush(stdout);
	fflush(stdout);
	#endif 
}




/**
 * @private buffIbcast
 * Broadcast for buffer
 * **/
void buffIbcast(Buff *buffer, int proc, Requests *requests, MPI_Comm comm, int mpi_rank){
	#if SPLITTED
	mpi_broadcast_splitted(
		buffer -> buffer,
		buffer -> size,
		proc, 
		comm, 
		requests, 
		mpi_rank,
		MPI_DOUBLE
	);
	#else 
	MPI_CHECK(
		MPI_Ibcast( 
			buffer -> buffer, 
			buffer -> size, 
			MPI_DOUBLE, 
			proc, comm, 
			&requests -> arr[proc] 
		) 
	);
	#endif 
}

/**
 * @private buffIbcast_quotient
 * Broadcast for buffer
 * **/
void buffIbcast_quotient(Buff *buffer, int proc, Requests *requests, MPI_Comm comm, int mpi_rank, unsigned int quotient){
	if (buffer -> size > 0){
		#if DEBUG
		printf("[%d | %d] Bcast buffer %d\n", mpi_rank, proc, buffer->size);
		#endif 
		requests -> size =1;
		*requests = requestsCalloc_abort(1, comm, mpi_rank);
		mpi_broadcast_splitted_quotient(
			buffer -> buffer,
			buffer -> size,
			proc, 
			comm, 
			requests, 
			mpi_rank,
			MPI_DOUBLE, quotient 
		);
	} else {
		#if DEBUG
		printf("[%d | %d] No data to bcast in buffer\n", mpi_rank, proc);
		#endif 
		requests -> size =0;
	}
	#if DEBUG
	printf("[%d | %d] Bcast buffer %d --> \n", mpi_rank, proc, buffer->size);
	fflush(stdout);
	#endif 
}
 


/***
 * @section MATRIXRANGE FUNCTIONS
 * **/
/**
 * @private matrix_range
 * mat declaration
 * **/
MatrixRange matrix_range(Matrix_type type){
	MatrixRange mat = {NULL, {0,0}, {0,0}, type};
	return mat;
}

/**
 * @private matrix_rangeNull
 * mat declaration
 * **/
MatrixRanges matrix_ranges(Matrix_type type, unsigned int total){
	unsigned int i;
	MatrixRanges mats = {total, NULL};
	mats.array=(MatrixRange *) calloc_abort(total, sizeof(MatrixRange), "matrix ranges", MPI_COMM_WORLD, -1);
	for (i = 0; i < total; i++){
		mats.array[i] = matrix_range(type);
	}
	return mats;
}

/**
 * @private matrix_rangeMat_calloc_abort
 * mat -> mat constructor
 * **/
void matrix_rangeMat_calloc_abort(
	MatrixRange *mat, const char *name, 
	MPI_Comm comm, int mpi_rank 
){
	size_t size = mat -> rows.size*mat -> columns.size;
	#if DEBUG > 1
	printf("size = %ld\n", size);
	#endif 
	fflush(stdout);
	mat -> mat = (double *) calloc_abort(size, sizeof(double), name, comm, mpi_rank);
}

/**
 * @private matrix_rangeCopy_null
 * Copies orig structure into new matrix and sets it's vlaues to 0
 * **/
MatrixRange matrix_rangeCopy_null(
	MatrixRange orig, 
	const char *name, 
	MPI_Comm comm, 
	int mpi_rank
){
	MatrixRange dest;
	size_t size = orig.rows.size*orig.columns.size;
	dest.columns =  orig.columns;
	dest.rows = orig.rows;
	dest.mat = (double *) calloc_abort(size, sizeof(double), name, comm, mpi_rank);
	return dest;
}
double matrix_rangeCheckSizeLimit(
	MatrixRange mat, int mpi_rank
){
	unsigned long int current_size = mat.rows.size*mat.columns.size;
	//char *info= (char*)calloc_abort(50 + current_size+BLOCK_LIMIT, sizeof(char), "info", MPI_COMM_WORLD, mpi_rank);
	//sprintf(info, "Block size %ld > Max bcast block size %d", current_size, BLOCK_LIMIT);
	char *info =  "Too few procs for this matrix range";
	if (current_size > BLOCK_LIMIT){
		error_mpi_abort(MPI_COMM_WORLD, INVALID_PROCS_NUM, mpi_rank, info);
	}
	return current_size;
}

/***
 * @private matrix_rangeIbcast
 * Broadcast for MatrixRange
 * **/
void matrix_rangeIbcast(
	MatrixRange *mat, 
	int proc, 
	Requests *requests, 
	MPI_Comm comm, 
	unsigned long int block_before, 
	int mpi_rank
){
	unsigned long int current_size=matrix_rangeCheckSizeLimit(*mat, mpi_rank);
	MPI_CHECK(
		MPI_Ibcast(
			mat -> mat + block_before, 
			current_size, 
			MPI_DOUBLE, 
			proc, 
			comm, 
			&requests -> arr[proc]
		) 
	);
}
void matrix_rangeISend(
	MatrixRange *mat, int comm_from, int comm_to, 
	Requests *requests, MPI_Comm comm, unsigned long int block_before, 
	unsigned int req_id
){
	//matrix_rangeCheckSizeLimit(*mat);
	MPI_CHECK(
		MPI_Isend(
			mat -> mat + block_before, 
			mat -> rows.size*mat -> columns.size, 
			MPI_DOUBLE, 
			comm_to, comm_from, comm, &requests -> arr[req_id]
		) 
	);
}

void matrix_rangeIRecv(
	MatrixRange *mat, int comm_from, int comm_to, 
	Requests *requests, MPI_Comm comm, unsigned long int block_before,
	unsigned int req_id
){
	MPI_CHECK(
		MPI_Irecv(
			mat -> mat + block_before, 
			mat -> rows.size*mat -> columns.size, 
			MPI_DOUBLE, 
			comm_from, comm_to, comm, &requests -> arr[req_id]
		) 
	);
}

/***
 * @private set_time
 * Function for setting current time for a time variable
 * */
double set_time(struct timeval time){
	double total_time, diff_sec, diff_usec;
	struct timeval aux_time;
	gettimeofday(&aux_time, NULL);
	diff_sec 	= (double) aux_time.tv_sec - time.tv_sec;
	diff_usec 	= (double) aux_time.tv_usec - time.tv_usec;
    total_time 	= (double) diff_sec + diff_usec/1000000.0;	
	return total_time;
}





/***
 * @private error_num_name
 * Gets name from error
 * **/

char *error_num_name(MM_TRIANGULAR_ERROR error_value, char *more_info){
	printf("ERRORNUMNAME\n");
	fflush(stdout);
	char *error = (char*) calloc_abort(strlen(more_info)+50,sizeof(char), "error name", MPI_COMM_WORLD, -1);
	switch (error_value)
	{
	case INVALID_PROCS_NUM:
		sprintf(error, "INVALID_PROCS_NUM %s", more_info);
		break;
	case INVALID_ARGS:
		sprintf(error, "INVALID_ARGS %s", more_info);
		break;
	case INVALID_COMMUNICATION_SHAPE:
		sprintf(error, "INVALID_COMMUNICATION_SHAPE %s", more_info);
		break;
	case ERROR_OPENING_FILE:
		sprintf(error, "ERROR_OPENING_FILE %s", more_info);
		break;
	case MEMORY_ERROR:
		sprintf(error, "MEMORY_ERROR %s", more_info);
		break;
	case INVALID_BLOCKS_NUM:
		sprintf(error, "INVALID BLOCKS NUM %s", more_info);
		break;
	case INVALID_ROWS_BEFORE:
		sprintf(error, "INVALID_ARGS %s", more_info);
		break;
    case EMPTY_BLOCK:
		sprintf(error, "EMPTY BLOCK %s", more_info);
		break;
    case HORIZONTAL_MATRIX:
		sprintf(error, "HORIZONTAL MATRIX %s", more_info);
		break;
    case INVALID_RANGE:
		sprintf(error, "INVALID RANGE %s", more_info);
		break;
    case INVALID_ITERATION_RANGE:
    	sprintf(error, "INVALID ITERATION RANGE %s", more_info);
		break;
	case INVALID_RECTANGLE_NUM:
		sprintf(error, "INVALID RECTANGLE NUM %s", more_info);
		break;
    case INVALID_SIZE:
		sprintf(error, "INVALID SIZE %s", more_info);
		break;
    case INVALID_INVERSION_INDEXES:
		sprintf(error, "INVALID INVERSION INDEXES %s", more_info);
		break;
    case INVALID_INITIAL_SHAPE:
		sprintf(error, "INVALID INITIAL SHAPE %s", more_info);
		break;
	default:
		return "UNKNOWN_ERROR";
		break;
	}
	return error;
}

FILE *fopen_abort(char *filename, char *mode, _Bool abort){
	FILE *f=NULL;
	size_t intentos = 0;
	do {
		f=fopen(filename, mode);
		intentos++;
	}while(f==NULL && intentos < 10);
	
	if (f==NULL){
		fprintf(stderr, "Couldn't open %s in mode %s\n", filename, mode);
		fflush(stdout);
			int mpi_rank;
			MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
			char *file=(char*) calloc(25+300+mpi_rank, sizeof(char*));
			sprintf(file, "Error_Writting_%d_%s", mpi_rank, filename);
		if (abort){
			MPI_CHECK(MPI_Abort(MPI_COMM_WORLD, ERROR_OPENING_FILE+2););
		} else {
			int written=0;
			FILE *f = fopen(file, "w");
				do {
					written = fprintf(f, "Error:%s: %d\n", error_num_name(ERROR_OPENING_FILE, filename), ERROR_OPENING_FILE+2);
				} while (written <= 0);
			fclose(f);
			 f = fopen(file, mode);
			} 
		
	}
	return f;
}

void ferror_mpi_abort(
	MPI_Comm comm, 
	MM_TRIANGULAR_ERROR error_value,
	int mpi_rank, 
	char *more_info
){
	char *error_filename = calloc_abort(ERROR_CHAR_SIZE, sizeof(char), "error_filename",  MPI_COMM_WORLD, mpi_rank);
	char *errorname = calloc_abort(30, sizeof(char), "errorname", MPI_COMM_WORLD, mpi_rank);
	int written=0;
	FILE *error_file_name, *errorf;
	
	sprintf(errorname, "error_filename_%d", mpi_rank );
	printf("[%d] --EE--> ERROR FILENAME Error: %d: %s <--EE--\n", mpi_rank, error_value+2, error_num_name(error_value, more_info));
	fflush(stdout);
	error_file_name = fopen_abort(errorname, "r", 0);
	if (error_file_name != NULL || strcmp( error_filename, "")) {
		printf("[%d] Getting error filename\n", mpi_rank); 
		fflush(stdout);
		fgets(error_filename, ERROR_CHAR_SIZE, error_file_name );
		printf("[%d] file: %s\n", mpi_rank, error_filename);
		fflush(stdout);
		fclose(error_file_name);
		errorf = fopen_abort(error_filename, "w", 0);
		do {
			fprintf(stderr, "--EE--> ERROR FILENAME Error: %d: %s <--EE--\n", error_value+2, error_num_name(error_value, more_info));
			written = fprintf(errorf, "Error:%s_%d\n", error_num_name(error_value, more_info), error_value+2);
		} while (written <= 0);	
		fclose(errorf);
	} else {
		printf("Empty filename %s \n", error_filename);
		fprintf(stderr, "--EE--> ERROR FILENAME Error: %d: %s <--EE--\n", error_value+2, error_num_name(error_value, more_info));
		fflush(stdout);
		MPI_CHECK(MPI_Abort(comm, error_value+2););
	} 
}



/***
 * @private error_mpi_abort
 * Aborts and saves error
 * **/
void error_mpi_abort(
	MPI_Comm comm, 
	MM_TRIANGULAR_ERROR error_value,
	int mpi_rank, 
	char *more_info
){
	printf("--> Error mpi abort\n");
	fflush(stdout);
	printf("--> ferror MPI_Abort");
	ferror_mpi_abort(comm, error_value, mpi_rank, more_info);
	printf("--> MPI_Abort");
	MPI_CHECK(MPI_Abort(comm, error_value+2););
	printf("--> MPI_Finalize\n");
	MPI_Finalize();
	printf("Error_mpi_abort --> \n");
	exit(error_value+2);
}

void write_error_filename(char *error_filename, int mpi_rank){
	char *errorname = calloc_abort(30, sizeof(char), "errorname", MPI_COMM_WORLD, mpi_rank);
	int num_chars=0;
	sprintf(errorname, "error_filename_%d", mpi_rank );
	if(strcmp( error_filename, "")){
		#if DEBUG 
		printf("wrtFilename %s\n", error_filename);
		#endif 
		FILE *error_file_name = fopen_abort(errorname,"w", 0);
		do {
			num_chars=fprintf(error_file_name, "%s", error_filename);
		} while (num_chars <= 0);
		fclose(error_file_name);
	} else {
		#if DEBUG
		printf("Empty filename %s\n", error_filename);
		#endif 
	}
}

double matrix_rangeInfinityNorm_(MatrixRange myMat, MatrixRange mat, MPI_Comm comm, int mpi_rank){
	MPI_Barrier(comm);
	Buff myRowsSumm;
	Buff rowsSumm;
	double norm = 0.0;
	int myRow, myColumn, row, last_column;
	myRowsSumm 	= buff(mat.rows.size, mpi_rank);
	rowsSumm 	= buff(mat.rows.size, mpi_rank);
	for (
		myRow = 0; 
		myRow < myMat.rows.size; 
		myRow++
	){
		row = myRow + myMat.rows.begin;
		myRowsSumm.buffer[row] = 0.0;
		if(myMat.type == TRIANGULAR_D) {
			//printf("Soy triangular\n");
			fflush(stdout);
			last_column = MAX(0,row-myMat.columns.begin+1);
		} else {
			//printf("Soy rectangular\n");
			fflush(stdout);
			last_column = myMat.columns.size;
		}
			
		for (
			myColumn = 0; 
			myColumn < last_column; 
			myColumn++
		){
			myRowsSumm.buffer[row] += fabs(myMat.mat[myRow * myMat.columns.size + myColumn]);
		}
	}
	MPI_Barrier(comm);
	#if DEBUG
	for (row = 0; row < mat.rows.size; row++){
		printf("[%d] myAdd[%d] = %.3lf \n", mpi_rank, row, myRowsSumm.buffer[row]);
		fflush(stdout);
	}
	#endif 
	MPI_CHECK( 
		MPI_Reduce( 
			myRowsSumm.buffer, 
			rowsSumm.buffer, 
			(int) rowsSumm.size, 
			MPI_DOUBLE, 
			MPI_SUM, 
			0, 
			MPI_COMM_WORLD 
		) 
	);
	if (mpi_rank == 0){
		norm = 0.0;
		for (row = 0; row < mat.rows.size; row++){
			//printf("FinalRowSum[%d] = %.3lf\n", row, rowsSumm.buffer[row]);
			norm = MAX(rowsSumm.buffer[row],norm);
		}
		FILE *normf = fopen("norm_ff.txt", "w+");
		printf ("\nNorm: %.3lf\n", norm); //For parser
		fprintf(normf, "%.3lf", norm); //For small check

	}
	return norm;

}

double matrix_rangeInfinityNorm(MatrixRange myMat, MatrixRange mat, MPI_Comm comm, int mpi_rank){
	#if DEBUG
	matrix_rangePrintf_rank_stage(myMat, "myMat", mpi_rank, -1, 0);
	#endif
	MPI_Barrier(comm);
	Buff myRowsSumm;
	Buff rowsSumm;
	double norm = 0.0;
	double norm_aux;
	int myRow=0, myColumn=0, row=0, last_column=0;
	//int grow, gcolumn;
	myRowsSumm 	= buff(mat.rows.size, mpi_rank);
	rowsSumm 	= buff(mat.rows.size, mpi_rank);
	#if DEBUG
	for (row = 0; row < mat.rows.size; row++){
		printf("[%d | %d] %.3lf\n", mpi_rank, row, myRowsSumm.buffer[row]) ;
		fflush(stdout);
	}
	#endif 
	
	
	row = myMat.rows.begin;
	for (
		myRow = 0; 
		myRow < myMat.rows.size; 
		myRow++
	){
		#if DEBUG
		printf(
			"Row: %d row | matRows %ld | myRowsSumm.size %d\n", 
			row, mat.rows.size, myRowsSumm.size
		);
		fflush(stdout);
		#endif 
		myRowsSumm.buffer[row] = 0.0;
		if(myMat.type == TRIANGULAR_D && mat.type == TRIANGULAR_D) {
			//printf("Soy triangular\n");
			fflush(stdout);
			last_column = MAX(0,row-myMat.columns.begin+1);
		} else {
			//printf("Soy rectangular\n");
			fflush(stdout);
			last_column = myMat.columns.size;
		}
			
		for (
			myColumn = 0; 
			myColumn < last_column; 
			myColumn++
		){
			#if DEBUG
			printf("[%d | %d] myRow %d myColumn %d Csize %ld Rsize %ld\n",
				mpi_rank, row, myRow, myColumn, myMat.columns.size, myMat.rows.size);
				fflush(stdout);
				#endif 
 			norm_aux=myMat.mat[myRow * myMat.columns.size + myColumn];
			#if DEBUG
			printf("[%d | %d] norm_auxPreFabs = %.3lf\n", 
				mpi_rank, myRow, norm_aux); 
			fflush(stdout);
			#endif 
			norm_aux=fabs(norm_aux);
			myRowsSumm.buffer[row] += norm_aux;
			#if DEBUG
			printf("[%d | %d] norm_aux = %.3lf\n", mpi_rank, myRow, norm_aux); 
			fflush(stdout);
			printf("[%d | %d] myRowAdd = %.3lf\n", mpi_rank, myRow, myRowsSumm.buffer[row]);
			fflush(stdout);
			#endif 
		}
		row++;
	}
	MPI_Barrier(comm);
	#if DEBUG
	printf("[%d] myRowSumm.size = %d\n", mpi_rank, myRowsSumm.size);
	fflush(stdout);
	
	for (row = 0; row < myRowsSumm.size; row++){
		printf("[%d] row: %d\n", mpi_rank, row); fflush(stdout);
		printf("[%d] myAdd: %.3lf\n", mpi_rank, myRowsSumm.buffer[row]); fflush(stdout);
		printf("[%d] myAdd[%d] = %.3lf \n", 
		mpi_rank, row, myRowsSumm.buffer[row]
		);
		fflush(stdout);
	}
	#endif 
	MPI_CHECK( 
		MPI_Reduce( 
			myRowsSumm.buffer, 
			rowsSumm.buffer, 
			(int) myRowsSumm.size, 
			MPI_DOUBLE, 
			MPI_SUM, 
			0, 
			MPI_COMM_WORLD 
		) 
	);
	if (mpi_rank == 0){
		norm = 0.0;
		for (row = 0; row < mat.rows.size; row++){
			//printf("FinalRowSum[%d] = %.3lf\n", row, rowsSumm.buffer[row]);
			norm = MAX(rowsSumm.buffer[row],norm);
		}
		FILE *normf = fopen("norm_ff.txt", "w+");
		printf ("\nNorm: %.3lf\n", norm); //For parser
		fprintf(normf, "%.3lf", norm); //For small check

	}
	return norm;

}
