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

/********************************************************
 * @section MAPPING FUNCTIONS 							*
 * @brief Mapping functions divides Range into blocks.	*
 * They are used for tilling data matrices.				*
 * ******************************************************/
/******
 * @subsection Usual mappings
 * */
/***
 * @private map_regular
 * @brief Mapping function. 
 * Blocks: 
 * * Regular blocks 
 * * Remainder in the last block.
 * **/
Range map_regular( 
	unsigned long int num_rows, 
	unsigned int num_procs, 
	unsigned int ind_proc,
	int mpi_rank
) {
	unsigned long int equitative_size = num_rows / num_procs;
	unsigned long int remainder = num_rows % num_procs;
	Range panel_ind_proc;
	//-- Sizes
	if (equitative_size < 1) equitative_size = 1;
	if (
			num_procs != 1 
		&& 	equitative_size*num_procs > num_rows
	) {
		printf(
			"Invalid num of processors %d for size %ld using map_regular\n",
			num_procs, num_rows
		);
		printf(
			"Eq_size = %ld, Eq_size*procs=%ld\n", 
			equitative_size, 
			equitative_size*num_procs
		);

		//fflush(stdout);
		error_mpi_abort(MPI_COMM_WORLD, INVALID_PROCS_NUM, mpi_rank, "Too few procs for map regular");
	}
	//-- Panel
	panel_ind_proc.begin = ind_proc * equitative_size;
	panel_ind_proc.size = equitative_size;
	//-- Add remainder to last proc's panel
	if ( ind_proc == num_procs - 1 ) {
		panel_ind_proc.size += remainder;
	}
	#if DEBUG
	rangePrintf_rank_stage("mapreg", panel_ind_proc, ind_proc, -1);
	//fflush(stdout);
	#endif 
	return panel_ind_proc;
}

/***
 * @private map_regular_first_double
 * @brief Mapping function: 
 * 	Blocks: 
 * 	* First block with double size, 
 * 	* the rest are regular 
 * 	* Remainder in the last block. 
 * **/
Range map_regular_first_double( 
	unsigned long int num_rows, 
	unsigned int num_procs, 
	unsigned int ind_proc,
	int mpi_rank
) {
	//-- Sizes
	unsigned long int equitative_size = num_rows / (num_procs+1);
	unsigned long int remainder = num_rows % (num_procs+1);
	//-- Panel
	Range panel_ind_procs;
	if (
			num_procs != 1 
		&& (equitative_size*(num_procs+1) > num_rows)
	){
		error_mpi_abort( MPI_COMM_WORLD, INVALID_PROCS_NUM, mpi_rank, "ERROR: invalid num of procs (total sizes excedes total rows");
	}
	//-- 
	if ( ind_proc == 0 ) {//First panel double-sized
		panel_ind_procs.begin = 0;
		panel_ind_procs.size = (num_procs == 1) ? num_rows : equitative_size*2;
	} else {//Regular panels: eq size
		panel_ind_procs.begin = (ind_proc+1) * equitative_size;
		panel_ind_procs.size = equitative_size;
	}
	//-- Add remainder to last proc's panel
	if ( 
			ind_proc == num_procs - 1 
		&& 	num_procs!= 1
	) {
		panel_ind_procs.size += remainder; //Errores de redondeo?
		//panel_ind_procs.size = num_rows-panel_ind_procs.begin;
	}
	//rangePrintf_rank_stage("regular_range", panel_ind_procs, mpi_rank, ind_proc);
	return panel_ind_procs;
}
/******
 * @subsection BALANCED MAPPING
 * */
/***
 * @private get_triangle_num_elements
 * Computes elements of a regular triangle from its rows number
 * **/
unsigned long int get_triangle_num_elements(
	unsigned long int triangle_rows
){
    unsigned int num_elements;
    num_elements = 0.5*(triangle_rows*(triangle_rows+1));
    return num_elements;
}

/***
 * @private get_matrix_num_elements_under_diagonal
 * Computes elements of a regular triangle from its rows number
 * **/
unsigned long int get_matrix_num_elements_under_diagonal(
    unsigned long int matrix_rows, 
    unsigned long int matrix_cols
){
    unsigned long int num_elements = 0;
    unsigned long int num_elements_triangle = 0; 
    unsigned long int num_elements_rectangle = 0;
    if (matrix_rows <= matrix_cols){
        num_elements_triangle = get_triangle_num_elements(matrix_rows);
    } else {
        num_elements_triangle = get_triangle_num_elements(matrix_cols);
        num_elements_rectangle = matrix_cols*(matrix_rows-matrix_cols);
    }
    num_elements = num_elements_triangle + num_elements_rectangle;
    return num_elements;
}

/***
 * @private get_ideal_block_size
 * Computes ideal blocks sizes
 * **/
double get_ideal_block_size(
    unsigned long int matrix_rows, 
    unsigned long int matrix_cols, 
    unsigned int num_blocks, 
	int mpi_rank
){
    double ideal_block_size;
    unsigned long int num_elements;
    num_elements = get_matrix_num_elements_under_diagonal(matrix_rows, matrix_cols);
    ideal_block_size = (double) num_elements/ (double)num_blocks;
    if (ideal_block_size < 1) {
        error_mpi_abort(MPI_COMM_WORLD, INVALID_BLOCKS_NUM, mpi_rank, "Ideal block size smaller than one: too much procs for this row number");
    }
    ideal_block_size = (ideal_block_size < 1) ? 1 : ideal_block_size;
    return ideal_block_size;
}

/***
 * @private get_rows_before
 * @brief Gets rows before of panel with block_index id
 * rows*(rows+1)/2 = (id-1)*size
 * rows^2+rows-2(id-1)*size=0
 * rows=-1+sqrt(1+8(id-1)*size)
 * */
unsigned long int get_rows_before(double ideal_block_size, int block_index){
    double numerator = -1+sqrt(1+8*(block_index*ideal_block_size));
    unsigned long int rows_before = (numerator > 0) ? (floor(numerator)/2) : 0;
    return rows_before;
}


int check_range(Range rg){
    int exit_val = EXIT_SUCCESS;
    if (rg.size <= 0){
        exit_val = INVALID_RANGE;
		error_mpi_abort(MPI_COMM_WORLD, INVALID_RANGE, 0, "invalid range: too much procs for this mapping mode");
    }
    return exit_val;
}

int check_ranges(Ranges rgs){
    int exit_val = EXIT_SUCCESS;
    Range rg1, rg2;
    unsigned int i = 0;
    if (rgs.total == 1){
        exit_val = check_range(rgs.array[0]);
    } else {
        for (i = 0; i < rgs.total-1 && (exit_val == EXIT_SUCCESS); i++){
            rg1 = rgs.array[i];
            rg2 = rgs.array[i+1];
            exit_val = check_range(rg1);
            if (exit_val == EXIT_SUCCESS){
                if (
                    rg2.begin != rg1.begin + rg1.size 
                ){
                    printf("Invalid ranges %d - %d\n", i, i+1);
                    printf(
                        "rg2.begin %ld != rg1.begin %ld + rg1.size %ld = %ld\n", 
                        rg2.begin, 
                        rg1.begin, 
                        rg1.size,
                        rg1.begin + rg1.size
                    );
        
                    //fflush(stdout);
                    exit_val = INVALID_ITERATION_RANGE;
                }
            } else {
                exit_val = INVALID_ITERATION_RANGE;
            }
        }
        if (exit_val == EXIT_SUCCESS){
            if ((!check_range(rgs.array[rgs.total-1])) == EXIT_SUCCESS ){
                exit_val = INVALID_ITERATION_RANGE;
                //fflush(stdout);
                error_mpi_abort(MPI_COMM_WORLD, exit_val, 0, "Invalid last range: debug specific case");
            }
        }
    }
    return exit_val;
}   


/***
 * @private map_balanced (recursivo)
 * @brief Mapping function
 * Blocks: 
 * * Balanced according to the number of non-zero elements 
 * * Each block contains the most symilar possible amount of 
 * * non-zero elements
 * **/
Range map_balanced_( 
	unsigned long int num_rows, 
	unsigned int num_procs, 
	unsigned int ind_proc 
) {
	//-- Panel
   	Range result;
	//-- Sizes
   	unsigned long int triangle_m = num_rows;
	unsigned long int rectangle_m = 0;
   	unsigned long int equitative_size;
	#if DEBUG
   		unsigned long int size_triangle, size_rectangle;
	#endif 
   	unsigned long int total_elements  = 0.5 * num_rows * ( num_rows + 1 );
	//-- Indexes
   	unsigned long int i;
	//-- Auxiliar
	double rows_formula ;
   	//-- Previous rows loop
	for ( i = 0; i <= fmin(ind_proc, num_procs-1); i++ ){
		equitative_size = total_elements/(num_procs - i); //Division entera -- por eso quito ceil
		rows_formula = -2.0 * rectangle_m -1 +
                        sqrt( 4 * rectangle_m * rectangle_m + 1
                        + 4* rectangle_m + 8 * equitative_size) ;
		triangle_m      = ceil(rows_formula) / 2;
      	// REMARK: It is necessary to check that we are not consuming more rows than possible
      	// If we do not perform next if, then it could happen (due to rounding) that a proc got 0 rows
      	if ( num_rows - rectangle_m - triangle_m < num_procs - i - 1){
			triangle_m -= (num_procs - i - 1) - (num_rows - rectangle_m - triangle_m);
      	}
      	if ( triangle_m > num_rows ){
        	triangle_m = num_rows;
	  	}
	  	#if DEBUG //Check sizes inside
	   		if(ind_proc == num_procs-1){
				size_rectangle = rectangle_m*triangle_m;
				size_triangle = (int) ceil(0.5 * triangle_m * ( triangle_m + 1));
				printf("[%ld] total_elements = %ld\n", i, total_elements);
				printf("[%ld] equitative_size = %ld\n", i, equitative_size);
				printf("[%ld] rectangle_m = %ld\n", i, rectangle_m);
				printf("[%ld] triangle_m = %ld\n", i, triangle_m);
				printf("[%ld] n_rectangle = %ld\n", i, size_rectangle);
				printf("[%ld] n_triangle = %ld\n", i, size_triangle );
				//printf("[%d] total_elements_next = %d - (%d+%d), %d\n", i, total_elements, size_rectangle, size_triangle, 
				//total_elements-size_rectangle-size_triangle);
				//fflush(stdout);
	   		}
	   	#endif
      	total_elements -= rectangle_m * triangle_m + 0.5 * triangle_m * ( triangle_m + 1 );
      	rectangle_m    += triangle_m;
		#if DEBUG > 1
		printf("[%ld] total_elements_next_calc = %ld\n", i, total_elements);
		//fflush(stdout);
		#endif 
   	}
   	// Check if it is last panel (would take all the remaining elements)
	if ( ind_proc == num_procs - 1 ){
      	result.begin = rectangle_m  - triangle_m;
      	result.size  = num_rows - rectangle_m + triangle_m;
   	}
   	else {
    	result.begin = rectangle_m - triangle_m;
      	result.size  = triangle_m;
   	}
   return result;
}
/***
 * @private map_balanced (directo)
 * @brief Mapping function
 * Blocks: 
 * * Balanced according to the number of non-zero elements 
 * * Each block contains the most symilar possible amount of 
 * * non-zero elements
 * **/
Range map_balanced(
	unsigned long int num_rows, 
	unsigned int num_procs, 
	unsigned int ind_proc,
	int mpi_rank
){
	#if USE_OLD_BALANCED
	return map_balanced_(num_rows, num_procs, ind_proc);
	#else
    Range result;
    unsigned int rows_begin, rows_end;
    unsigned int ideal_equitative_size = ceil(get_ideal_block_size(num_rows, num_rows, num_procs, mpi_rank));
    #if DEBUG
    printf("\t>> MAP BALANCED\n");
    //fflush(stdout);
    printf("\t\t Computing rows begin\n");
    //fflush(stdout);
    #endif 
    rows_begin  = get_rows_before(ideal_equitative_size, ind_proc);
    #if DEBUG
    printf("\t\t Computing rows end\n");
    //fflush(stdout);
    #endif 
    rows_end    = get_rows_before(ideal_equitative_size, ind_proc+1);
    result.begin    = rows_begin;
    result.size     = rows_end-rows_begin;
	//Reajuste ultimo y primero
	if (ind_proc==0){
		result.begin=0;
		if (num_procs > 1){
			result.size=get_rows_before(ideal_equitative_size, 1);
		} else {
			result.size=num_rows;
		}
	}
	if (num_rows > 1 && ind_proc == num_procs - 1){
		result.size=num_rows-result.begin;
	}
	if(check_range(result) != EXIT_SUCCESS){
		error_mpi_abort(MPI_COMM_WORLD, INVALID_PROCS_NUM, mpi_rank, "Invalid range. Debug specific case");
	}
    #if DEBUG
    printf("\t\t[%d] range.begin: %ld range.size: %ld\n", ind_proc, result.begin, result.size);
	printf("\tMAP BALANCED <<\n");
    #endif
    return result;
	#endif 
}




/***
 * @private  @section DIRECTIONED FUNCTIONS
 * */
void rangeFlip(
	Range *rg, 
	unsigned int num_procs
){
	rg -> begin=num_procs-rg -> begin;
}

Range map_balanced_non_regular(
    unsigned long int num_rows, 
	unsigned long int num_cols,
    int num_procs, 
    int ind_proc,
	Direction dir
){
	Range block;
	int mpi_rank = ind_proc;
	unsigned int block_index=ind_proc=0;
    unsigned int rectangular_block_rows, rectangular_block_index=0;
    unsigned int triangular_block_num=0;
	unsigned int rectangular_block_num=0;
    unsigned long int remainder = 0;
	unsigned long int ideal_equitative_size=0;
	if (dir == HORIZONTAL){
		block_index=num_procs-block_index;
	}
	if (num_rows > num_cols){
		if (num_procs == 1){
        	if (dir == VERTICAL){
            	error_mpi_abort(MPI_COMM_WORLD, INVALID_PROCS_NUM, mpi_rank, "For vertical map balanced sizes triangular and rectangular, we need at least 2 processors");
        	} else {
            	block.begin = 0;
            	block.size = num_cols;
        	}
    	} else {
			rectangular_block_num = 0;
			ideal_equitative_size = get_ideal_block_size(num_rows, num_cols, num_procs, mpi_rank);
			rectangular_block_rows = ideal_equitative_size/num_cols;
			double aux= ((double)num_rows-num_cols)/(double)rectangular_block_rows;
			rectangular_block_num=aux > 1 ? floor (aux) : 1;
		}	if (rectangular_block_num == num_procs){
            error_mpi_abort(MPI_COMM_WORLD, INVALID_RECTANGLE_NUM, ind_proc, "Rectangles wants all for them, give more procs, please\n");
        }
		if (rectangular_block_num == 1 && ind_proc) {
            block.begin = num_cols+remainder;
            block.size = rectangular_block_rows;
        } else if (rectangular_block_num >= num_procs  ){
            //printf("No puede salir mas rectangulos que numprocesos ni igual");
            //printf("Recs: %d nprocs: %d\n", rectangular_block_num, num_procs);
            //exit(EXIT_FAILURE);
			error_mpi_abort(MPI_COMM_WORLD, INVALID_RECTANGLE_NUM, ind_proc, "Rectangles wants all (or more) for them, give more procs, please\n");
        } else {
			if (
						(block_index >= num_procs-rectangular_block_num+1)
					&& 	(block_index <= num_procs)
			){
                rectangular_block_index = block_index-num_procs+rectangular_block_num;
                block=map_regular(num_rows-num_cols, rectangular_block_num, rectangular_block_index, mpi_rank);
                block.begin += num_cols;
            }
            
            //Primero y ultimo,    
			if (block_index==num_procs-rectangular_block_num){
            	if (dir == VERTICAL) {
                	block.begin = num_cols;
            	} else {
                	if (num_cols-block.begin > 0){
	                    block.begin = num_cols;
    	                block.size -= num_cols-block.begin;
        	        }
	            }
				if (rectangular_block_num > 1 && rectangular_block_index==0){
					Range block_aux=map_regular(num_rows-num_cols, rectangular_block_num, 1, mpi_rank);
	            	block.size = block_aux.begin - num_cols;

				}
			}
			if (block_index == num_procs-1){
            	block.size = num_rows-block.begin;
			}
        } 

        //Reajuste tamaño
        #if DEBUG
        printf("START TRIANGULAR BLOCKS\n");
        //fflush(stdout);
        #endif 
        triangular_block_num = num_procs-rectangular_block_num;
        //ideal_equitative_size_triangles = get_ideal_block_size(num_cols, num_cols, triangular_block_num);
    
        //Calculo paneles trapezoidales
		if (block_index >= 1 && block_index < triangular_block_num){
		    block = map_balanced(num_cols, num_procs, block_index, mpi_rank);
        }

		if (block_index==0){
        	block.size = map_balanced(num_cols, num_procs,1, mpi_rank).begin;
		}
		if (block_index==triangular_block_num-1){
        	block.size = num_cols - block.begin;
		}
        if (dir == HORIZONTAL && block_index == triangular_block_num-1){//Enlosamiento
			Range block_aux_= map_balanced(num_cols, num_procs,block_index+1, mpi_rank);
            block.size = block_aux_.begin - block.begin;
        }
	} else {
		block=map_balanced(num_rows, num_procs, ind_proc, mpi_rank);
	}
	if (dir == HORIZONTAL){
		rangeFlip(&block, num_procs);
	}
	return block;
}
