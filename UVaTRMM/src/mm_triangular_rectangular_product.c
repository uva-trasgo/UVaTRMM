#include "mm_triangular.h"



void mm_triangular_rectangular_productRECTANGLE(
	MatrixRange mat_A, 
    MatrixRange mat_B,
    MatrixRange *mat_C,
    size_t block_before_C,
    unsigned int stage,
    unsigned int mpi_rank
){
    // Variables for MKL calls
	//-- C = alpha A*B+beta C
    double ALPHA=1.0; 
    double BETA	=0.0;
	//-- Dgemm rectangle * rectangle
	const CBLAS_LAYOUT 		layout = CblasRowMajor;
	const CBLAS_TRANSPOSE 	TransA = CblasNoTrans;
	const CBLAS_TRANSPOSE 	TransB = CblasNoTrans; 
	CBLAS_INDEX M, N, K;
	CBLAS_INDEX LDA, LDB, LDC; 
	//--- Set up sizes and ranges
	//--- Original rectangles size
	size_t mat_A_rectangle_rows 	= 0;
	size_t mat_A_rectangle_columns  = mat_A.rows.begin+1;//mat_A->rows.begin+1; 
	#if DEBUG
	printf("[%d | %d] Mat_A_RectangleColumns %d\n", mpi_rank, stage, mat_A_rectangle_columns); 
	fflush(stdout);
	#endif 
	size_t block_before_B=0;
	//block_before_B=mat_A.rows.begin;
	size_t mat_B_rectangle_rows     = mat_B.rows.size;//-block_before_B;
	//block_before_B=block_before_B*mat_B.columns.size;
	if (mat_A.rows.begin == 0 ){
		mat_A_rectangle_rows= mat_A.rows.size;
		//mat_A_rectangle_rows = 1; --> Cambiado 
	}else {
		mat_A_rectangle_rows=mat_A.rows.size;
	}
	#if DEBUG
	printf("[%d | %d] mat_A_rectangles_rows=%ld\n", mpi_rank, stage, mat_A_rectangle_rows);
	#endif 
	#if DEBUG > 1
    size_t mat_B_rectangle_columns 	= mat_B.columns.size;
    size_t index1;
	#endif 
	
	#if DEBUG
		printf("[%d | %d] cols_recA %ld rows_recB %ld\n", mpi_rank,stage ,mat_A_rectangle_columns, mat_B_rectangle_rows);
		//matrix_rangePrintf_rank_stage(*mat_A, "fullA",mpi_rank, stage, 0);
	#endif 
	//mat_A_rectangle_columns = MIN(mat_A_rectangle_columns, mat_B_rectangle_rows);
	mat_B_rectangle_rows = MIN(mat_A_rectangle_columns, mat_B_rectangle_rows);
	#if DEBUG  > 1
		for (long int f = 0; 
			f < mat_A.rows.size;
			f++
		){
			for (
				long int c = mat_B.rows.begin; 
				c < mat_B.rows.begin+mat_A_rectangle_columns; 
				c++
			){
				printf("[%d | %d] recA[%ld][%ld] = %.0lf\n", 
					mpi_rank, stage, 
					f + mat_A.rows.begin, 
					c + mat_A.columns.begin, 
					mat_A.mat[f* mat_A.columns.size+c]
				);
			}
			printf("\n[%d | %d] -- | --\n", mpi_rank, stage);
		}
		matrix_rangePrintf_rank_stage(mat_A, "fullAhm",mpi_rank, stage, 0);
	#endif 
	// De A tenemos que coger 
	// Desde 0 hasta begin+1
	// De B tenemos que coger [MIN(Filas B, Columnas A)] x [Columnas B]
	#if DEBUG>1
		for (
			long int f = 0; f < mat_A_rectangle_columns; f++){
			for (long int c = 0; c < mat_B_rectangle_columns; c++){
				index1=f*mat_B.columns.size+c;
				printf(
					"[%d | %d] recB[%ld][%ld] = %.0lf\n", 
					mpi_rank, stage, 
					f + mat_B.rows.begin, 
					c + mat_B.columns.begin, 
					mat_B.mat[index1]
				);
			}
		}
		printf("\n[%d | %d] -- | --\n", mpi_rank, stage);
	#endif 
	// C ababa con [Filas A] x [Columnas B]
	M = mat_A_rectangle_rows;	//mat_A -> rows.size; 
	N = mat_B.columns.size;	//mat_B -> columns.size;
	K = mat_A_rectangle_columns;
	#if DEBUG 	
	printf("[%d | %d] Rect M (A rows)= %ld\tN= %ld (B columns) \tK (A columns) = %ld\n", mpi_rank, stage, M, N, K);
	fflush(stdout);
	#endif 
	//Avoid errors of 0-dimension matrices
	if (M > 0 && N > 0 && K > 0){
		//Intercambiados
		LDA = mat_A.columns.size; //Leading dimension A
		LDB = mat_B.columns.size; //Leading dimension B
		LDC = mat_C -> columns.size; //Leading dimension C
		#if DEBUG
			printf("[%d | %d] Rect M = %ld\tN= %ld\tK = %ld\n", mpi_rank, stage, M, N, K);
			fflush(stdout);
			fflush(stdout);
    	matrix_rangePrintf_rank_stage(*mat_C, "Cbdgemm", mpi_rank, stage, 0);
		#endif 
		cblas_dgemm (  
			layout, TransA, TransB,
		    M, N, K, 
		    ALPHA, 
		    mat_A.mat,// + mat_A -> rows.begin, 
		    LDA,
            mat_B.mat+block_before_B, 
            LDB,BETA, 
            mat_C->mat+block_before_C, LDC
	    );
		#if DEBUG
		printf("[%d | %d] M = %ld\tN= %ld\tK = %ld\n", mpi_rank, stage, M, N, K);
		fflush(stdout);
		fflush(stdout);
		matrix_rangePrintf_rank_stage(mat_A, "Aadgemm", mpi_rank, stage, 0);
		matrix_rangePrintf_rank_stage(mat_B, "Badgemm", mpi_rank, stage, 0);
		matrix_rangePrintf_rank_stage(*mat_C, "Cadgemm", mpi_rank, stage, 0);
		#endif 
	}
}

void mm_triangular_rectangular_productTRIANGLE(
	MatrixRange mat_A, 
    MatrixRange mat_B,
    MatrixRange *mat_C,
    size_t block_before_C,
    unsigned int stage,
    unsigned int mpi_rank
){
	#if DEBUG 
		printf("[%d | %d] --> Product triangle\n", mpi_rank, stage);
	#endif 
    // Variables for MKL calls
	//-- C = alpha A*B+beta C
    double ALPHA=1.0; 
	const CBLAS_LAYOUT 		layout = CblasRowMajor;
	const CBLAS_TRANSPOSE 	TransA = CblasNoTrans;
	//const CBLAS_TRANSPOSE 	TransB = CblasNoTrans;
	//-- Dtrmm triangle LU * rectangle
	const CBLAS_SIDE sideA 	= CblasLeft;		//A is left (A*B)
	const CBLAS_UPLO triagA = CblasLower;		//A is triangular L
	const CBLAS_DIAG unit 	= CblasNonUnit;		//A is not Identity (1)
	CBLAS_INDEX M, N;
	CBLAS_INDEX LDA, LDB; 

	//--- Set up sizes and ranges
	//--- Original rectangles size
    ///Negative index ended up in problems 
	unsigned int triangle_first_pos, triangle_last_pos;
	size_t mat_A_triangle_rows;
	size_t mat_A_triangle_columns;
	size_t mat_B_rectangle_columns;
	size_t rows_begin_B;//, space_after_triangle_init_B;
	size_t space_above_triangle, space_left_to_triangle;
	size_t block_before_A, block_before_B;
	size_t index1, index2;
	
	/*if (mat_A.rows.begin==0) {
		triangle_first_pos			= 0;
		triangle_last_pos			= mat_A.rows.size;
		mat_A_triangle_rows			= mat_A.rows.size;
	} else {*/
		triangle_first_pos			= mat_A.rows.begin+1;
		triangle_last_pos			= mat_A.rows.begin+mat_A.rows.size;
		mat_A_triangle_rows			= MAX(0,triangle_last_pos-triangle_first_pos);
	//}
	mat_A_triangle_rows		= 	MIN(mat_A_triangle_rows, mat_A.columns.size);
	mat_A_triangle_rows		=	MIN(mat_A_triangle_rows, mat_B.rows.size);
	mat_A_triangle_columns 	= 	mat_A_triangle_rows; 
	mat_B_rectangle_columns = 	MIN(mat_B.rows.begin+mat_B.rows.size, mat_B.columns.size);
	space_left_to_triangle 	= 	MAX(0,mat_A.rows.begin+mat_A.rows.size - mat_A_triangle_rows);
	rows_begin_B 			= 	space_left_to_triangle;
	space_above_triangle	=	MAX(0,mat_A.rows.size-mat_A_triangle_rows);
	
	block_before_B = mat_B.columns.size*rows_begin_B;//Filas encima del triangulo
	block_before_A = space_above_triangle*mat_A.columns.size+space_left_to_triangle;
	
	M = mat_A_triangle_columns; 	//mat_A_triangle_rows;			//B rows
	N = mat_B_rectangle_columns;	//B columns
	#if DEBUG 
	//printf("[%d | %d] Range (%d,%d) Space left to triangle=%d\n",
	//	mpi_rank, stage, mat_A.rows.begin, mat_A.rows.begin+mat_A.rows.size-1, space_left_to_triangle
	//);
	printf("[%d | %d] Range (%ld,%ld) TrRows %ld Up %ld Left %ld | block before A %ld\n", 
		mpi_rank, stage, mat_A.rows.begin, mat_A.rows.begin+mat_A.rows.size-1, 
		mat_A_triangle_rows, space_above_triangle, space_left_to_triangle, block_before_A
	);
	fflush(stdout);
	rangePrintf_rank_stage("Arows", mat_A.rows, mpi_rank, stage);
	rangePrintf_rank_stage("Acolss", mat_A.columns, mpi_rank, stage);
	fflush(stdout);
	printf("Triange rows: %ld\n", mat_A_triangle_rows);
	fflush(stdout);
		printf("[%d | %d] Triangle M %ld N %ld\n", mpi_rank, stage, M, N);
	#endif
	//Avoid errors of 0-dimension matrices
	if (M > 0 && N > 0){
		#if DEBUG
		printf("mat_A_triangle_rows %ld\n", mat_A_triangle_rows);
		fflush(stdout);
		printf("[%d | %d] block before A %ld\n", mpi_rank, stage, block_before_A);
		printf("[%d | %d] block before B %ld\n", mpi_rank, stage, block_before_B);
		printf("[%d | %d] B rows %ld\n", mpi_rank, stage, M);
		printf("[%d | %d] B columns %ld\n", mpi_rank, stage, N);

		matrix_rangePrintf_rank_stage(mat_A, "matAbeforetr", mpi_rank, stage, 0);
		matrix_rangePrintf_rank_stage(mat_B, "matBbeforetr", mpi_rank, stage, 0);
		fflush(stdout);
		#endif 
		LDA = mat_A.columns.size;		//Leading dimension A
		LDB = mat_B.columns.size;		//Leading dimension B
		cblas_dtrmm (   
			layout, sideA, triagA, TransA, unit,
       		M, N, ALPHA, 
			mat_A.mat + block_before_A, 
			LDA,
			mat_B.mat + block_before_B, 
			LDB 
		);
		#if DEBUG
			matrix_rangePrintf_rank_stage(mat_B, "matBaftertr", mpi_rank, stage, 0);
			printf(
				"[%d | %d]  [%ld,%ld] x [%d,%ld]\n", 
				mpi_rank, stage, 
				mat_A.rows.begin+1, 
				mat_A.rows.begin+1+mat_A_triangle_rows,
				0, mat_B.columns.size
			);
		#endif 
		//-- Setup C from product result
		//#pragma omp parallel for collapse(2)
		for ( 
			long int f = rows_begin_B; 
			f < rows_begin_B+mat_A_triangle_columns; 
			f++ 
		){
			for ( 
				long int c = 0; 
				c < mat_B_rectangle_columns; 
				c++
			){
				//index1 = f*mat_C -> columns.size+c;
				//index1 = (f+mat_B.rows.begin)*mat_C -> columns.size+c;
				index2 = f*mat_B.columns.size+c;
				index1=index2;
				//index2=f+mat_A->rows.begin)*mat_B -> columns.size+c;
				//printf("[%d] [f][c]=[%ld][%ld]\n", mpi_rank, f, c);
				//fflush(stdout);
				#if DEBUG > 1
				printf("[%d | %d] f %ld c %ld id1 %ld id2 %ld Bcols %ld\n", 
				mpi_rank, stage, f, c, index1, index2, mat_B_rectangle_columns);
				printf("[%d | %d] C[%ld] = %.3lf\n", mpi_rank, stage, index1, mat_C -> mat[index1]);
				printf("[%d | %d] B[%ld] = %.3lf\n", mpi_rank, stage, index2, mat_B.mat[index2]);
				fflush(stdout);
				#endif 
				mat_C -> mat[index1] +=  mat_B.mat[index2];

				#if DEBUG > 1
				printf("[%d | %d] C[%ld] = %.3lf\n", mpi_rank, stage, index2, mat_B.mat[index2]);
				fflush(stdout);
					printf("[%d | %d] triang * rec [%ld][%ld]\n", mpi_rank, stage,f, c);
					fflush(stdout);
				#endif 
			}
		}
		#if DEBUG > 1
		printf("[%d | %d] |_\n", mpi_rank, stage);
		#endif 
	}
	#if DEBUG 
		printf("[%d | %d] Product triangle --> \n", mpi_rank, stage);
	#endif
}


void mm_rectangular_rectangular_product(
	MatrixRange mat_A, 
    MatrixRange mat_B,
    MatrixRange *mat_C,
    size_t block_before_C,
    unsigned int stage,
    unsigned int mpi_rank
){
	#if DEBUG
	printf("[ %d | %d ] Just before product\n",mpi_rank, stage);
	fflush(stdout);
	#endif 
    // Variables for MKL calls
	//-- C = alpha A*B+beta C
    double ALPHA=1.0; 
    double BETA	=0.0;
	//-- Dgemm rectangle * rectangle
	const CBLAS_LAYOUT 		layout = CblasRowMajor;
	const CBLAS_TRANSPOSE 	TransA = CblasNoTrans;
	const CBLAS_TRANSPOSE 	TransB = CblasNoTrans; 
	CBLAS_INDEX M = mat_A.rows.size;
	CBLAS_INDEX N = mat_B.columns.size; 
	CBLAS_INDEX K = mat_A.columns.size;
	CBLAS_INDEX LDA = mat_A.columns.size;
	CBLAS_INDEX LDB = mat_B.columns.size;
	CBLAS_INDEX LDC = mat_C -> columns.size;
	//Avoid errors of 0-dimension matrices
	if (M > 0 && N > 0 && K > 0){
		#if DEBUG 
			printf("[%d | %d] M = %ld\tN= %ld\tK = %ld\n", mpi_rank, stage, M, N, K);
			fflush(stdout);
		#endif 
		
		//Intercambiados
		LDA = mat_A.columns.size; //Leading dimension A
		LDB = mat_B.columns.size; //Leading dimension B
		LDC = mat_C -> columns.size; //Leading dimension C
		#if DEBUG > 1
    		matrix_rangePrintf_rank_stage(*mat_C, "Cbdgemm", mpi_rank, stage, 0);
			fflush(stdout);
		#endif 
		cblas_dgemm (  
			layout, TransA, TransB,
		    M, N, K, 
		    ALPHA, 
		    mat_A.mat,// + mat_A -> rows.begin, 
		    LDA,
            mat_B.mat, 
            LDB,BETA, 
            mat_C->mat+block_before_C, LDC
	    );
		#if DEBUG > 1
		matrix_rangePrintf_rank_stage(*mat_C, "Cadgemm", mpi_rank, stage, 0);
		fflush(stdout);
		#endif 
	} else {
		#if DEBUG > 1
		printf("0-dimension matrix\n");
		#endif 
	}
}

void mm_triangular_triangular_productRectangle(
    MatrixRange mat_A, 
    MatrixRange mat_B,
    MatrixRange *mat_C,
    size_t block_before_C,
    unsigned int stage,
    unsigned int mpi_rank
){
    // Variables for MKL calls
	//-- C = alpha A*B+beta C
    double ALPHA=1.0; 
    double BETA	=1.0;
	//-- Dgemm rectangle * rectangle
	const CBLAS_LAYOUT 		layout = CblasRowMajor;
	const CBLAS_TRANSPOSE 	TransA = CblasNoTrans;
	const CBLAS_TRANSPOSE 	TransB = CblasNoTrans; 
	CBLAS_INDEX M, N, K;
	CBLAS_INDEX LDA, LDB, LDC; 
	//--- Set up sizes and ranges
	//--- Original rectangles size
	//size_t mat_A_rectangle_rows		= 0;
	size_t mat_A_rectangle_columns 	= 0;
	if (mat_A.rows.begin!=0){
		mat_A_rectangle_columns=MIN(mat_A.columns.size,mat_A.rows.begin+1); //(mpi_rank == stage) ? 1 : mat_B.rows.size; 
		mat_A_rectangle_columns=MAX(0, mat_A_rectangle_columns - mat_B.rows.begin);
		mat_A_rectangle_columns=MIN(mat_A_rectangle_columns, mat_B.rows.size);
	} 
	#if DEBUG
	printf("[%d | %d] r(b%ld,s%ld) Rectangle mat_A_Columns %ld \n", mpi_rank, stage, mat_A.rows.begin, mat_A.rows.size, mat_A_rectangle_columns);
	#endif 
	/*if (mat_A.rows.begin == 0 ){
		mat_A_rectangle_rows = 1;
	}else {
		mat_A_rectangle_rows=mat_A.rows.size;
	}*/
	unsigned long int block_before_A 	= mat_B.rows.begin;
	// De B tenemos que coger [MIN(Filas B, Columnas A)] x [Columnas B]
    size_t mat_B_rectangle_columns 		= MIN(mat_B.rows.begin+mat_B.rows.size, mat_B.columns.size);
	size_t mat_B_rectangle_rows			= mat_B.rows.size;
	mat_B_rectangle_rows = MIN(mat_A_rectangle_columns, mat_B_rectangle_rows);
		// C ababa con [Filas A] x [Columnas B]
		M = mat_A.rows.size; 
		N = mat_B_rectangle_columns;
		K = mat_A_rectangle_columns;
		//Avoid errors of 0-dimension matrices
		if (M > 0 && N > 0 && K > 0){
			//Intercambiados
			LDA = mat_A.columns.size;	//Leading dimension A
			LDB = mat_B.columns.size; 	//Leading dimension B
			LDC = mat_C -> columns.size;
			#if DEBUG > 1
    		printf("[%d | %d] Rectangle block_before_A %ld block_before_C %ld M %ld N %ld K %ld Al %.3lf Bet %.3lf LDA %ld LDB %ld LDC %ld\n", 
			mpi_rank, stage, 
			block_before_A, block_before_C, M, N, K, ALPHA, BETA, LDA,LDB, LDC);
			matrix_rangePrintf_rank_stage(mat_A, "matA", mpi_rank, stage, 0);
			matrix_rangePrintf_rank_stage(mat_B, "matB", mpi_rank, stage, 0);
			matrix_rangePrintf_rank_stage(*mat_C, "matC", mpi_rank, stage, 0);
			fflush(stdout);
			#endif 
			cblas_dgemm (   
	    		layout, TransA, TransB,
			    M, N, K, 
			    ALPHA, 
			    mat_A.mat + block_before_A, //+ mat_A.rows.begin, 
			    LDA,
                mat_B.mat,
                LDB,BETA, 
                mat_C->mat+block_before_C, LDC
		    );
		}

}

void mm_triangular_triangular_productTriangle_(
    MatrixRange mat_A, 
    MatrixRange mat_B,
    MatrixRange *mat_C,
    size_t block_before_C,
    unsigned int stage,
    unsigned int mpi_rank
){

    // Variables for MKL calls
    double ALPHA=1.0; 
	const CBLAS_LAYOUT 		layout = CblasRowMajor;
	const CBLAS_TRANSPOSE 	TransA = CblasNoTrans; 
	CBLAS_INDEX M, N;
	CBLAS_INDEX LDA, LDB; 
	const CBLAS_SIDE sideA 	= CblasLeft;		//A is left (A*B)
	const CBLAS_UPLO triagA = CblasLower;		//A is triangular L
	const CBLAS_DIAG unit 	= CblasNonUnit;		//A is not Identity (1)

	size_t mat_A_triangle_rows 		= (mat_A.rows.size == 0) ? 0 : mat_A.rows.size - 1;
	size_t mat_B_rectangle_columns 	= mat_A.rows.begin+mat_A.rows.size;
	size_t block_before_A = mat_B.columns.size+mat_A.rows.begin+1;
	size_t block_before_B = mat_B.columns.size;

	M = mat_A_triangle_rows;		//A rows
	N = mat_B_rectangle_columns;	//B columns

	//M = mat_A.rows.size;
	//N = mat_B.columns.size;
	if (M > 0 && N > 0){
		LDA=mat_A.columns.size;
		LDB=mat_B.columns.size;
		cblas_dtrmm ( 
			layout, sideA, triagA, TransA, unit,
			M, N,ALPHA, 
			mat_A.mat+block_before_A, 
			LDA,mat_B.mat+block_before_B, 
			LDB 
		);
	}

	///Negative index ended up in problems 
	/*size_t mat_A_triangle_rows 		= (mat_A -> rows.size == 0) ? 0 : mat_A -> rows.size - 1;
	size_t mat_B_rectangle_columns 	= mat_A -> rows.begin+mat_A -> rows.size;
	size_t block_before_A = mat_B -> columns.size+mat_A -> rows.begin+1;
	M = mat_A_triangle_rows;			//A rows
	N = mat_B_rectangle_columns;	//B columns
		//Avoid errors of 0-dimension matrices
		if (M > 0 && N > 0){
			size_t block_before_B = mat_B -> columns.size;//mat_B -> columns.size*mat_A->rows.begin+mat_B->columns.size;
			printf("[%d | %d] block before B %ld\n", mpi_rank, stage, block_before_B);
			LDA = mat_A -> columns.size;		//Leading dimension A
			LDB = mat_B -> columns.size;		//Leading dimension B
			cblas_dtrmm (   
				layout, sideA, triagA, TransA, unit,
           		M, N, ALPHA, 
				mat_A -> mat+block_before_A, 
				LDA,
				mat_B -> mat + block_before_B, 
				LDB 
			);
			#if DEBUG
			for ( i = 0; i < mat_A_triangle_rows; i++ ){
				for ( j = 0; j < mat_B -> columns.size; j ++ ){
					printf(
						"[%d] Bcopyafter[%ld][%ld]=%.0lf\n",
						mpi_rank, 
						i + 1, //Real coordinate
						j,
						mat_B -> mat[i*mat_B_rectangle_columns+j]
					);
				}
			}
		#endif 
				//-- Setup C from product result
			for ( long int f = 1; f < 1+mat_A_triangle_rows; f++ ){
				for ( long int c = 0; c < mat_B_rectangle_columns; c++){
					index1 = f*mat_C -> columns.size+c;
					index2 = f*mat_B_rectangle_columns+c;//(f-1)*mat_B_rectangle_columns+c;
					mat_C -> mat[index1] +=  mat_B -> mat[index2];
					#if DEBUG
					printf(
						"[%d | %d] triang * rec [%ld][%ld]\n", mpi_rank, stage,
						f, 
					c
					);
					#endif 
					}
				#if DEBUG
				printf("[%d | %d] |_\n", mpi_rank, stage);
				#endif 
			}
		}
			#if DEBUG
				printf("[%d | %d] AFTER TRIANGULAR PART\n", mpi_rank, stage);
				fflush(stdout);
				matrix_rangePrintf_rank_stage(*mat_C, "TLC", mpi_rank, stage, 0);
			#endif */
}


void mm_triangular_triangular_productTriangle(
    MatrixRange mat_A, 
    MatrixRange mat_B,
    unsigned int stage,
    unsigned int mpi_rank
){
    // Variables for MKL calls
    double ALPHA=1.0; 
	const CBLAS_LAYOUT 		layout = CblasRowMajor;
	const CBLAS_TRANSPOSE 	TransA = CblasNoTrans; 
	CBLAS_INDEX M, N;
	CBLAS_INDEX LDA, LDB; 
	const CBLAS_SIDE sideA 	= CblasLeft;		//A is left (A*B)
	const CBLAS_UPLO triagA = CblasLower;		//A is triangular L
	const CBLAS_DIAG unit 	= CblasNonUnit;		//A is not Identity (1)
	unsigned int triangle_first_pos, triangle_last_pos;
	size_t mat_A_triangle_rows;
	size_t mat_A_triangle_columns;
	size_t mat_B_rectangle_columns;
	size_t rows_begin_B, space_after_triangle_init_B;
	size_t space_above_triangle, space_left_to_triangle;
	size_t block_before_A, block_before_B;
	
	//size_t index1, index2;
	
	if (mat_A.rows.begin==0) {
		triangle_first_pos			= 0;
		triangle_last_pos			= mat_A.rows.size;
		mat_A_triangle_rows			= mat_A.rows.size;
	} else {
		triangle_first_pos			= mat_A.rows.begin+1;
		triangle_last_pos			= mat_A.rows.begin+mat_A.rows.size;
		mat_A_triangle_rows			= MAX(0,triangle_last_pos-triangle_first_pos);
	}
	mat_A_triangle_rows			= MIN(mat_A_triangle_rows, mat_A.columns.size);
	mat_A_triangle_rows			= MIN(mat_A_triangle_rows, mat_B.rows.size);
	#if DEBUG
	printf("[%d | %d] r(b%ld,s%ld) | mat_A triangle rows:%ld (F%d,L%d)\n", mpi_rank, stage, mat_A.rows.begin, mat_A.rows.size, mat_A_triangle_rows, triangle_first_pos, triangle_last_pos);
	#endif 
	space_left_to_triangle 		= MAX(0,mat_A.rows.begin+mat_A.rows.size - mat_A_triangle_rows);
	
	rows_begin_B 				= space_left_to_triangle - mat_A.rows.begin;
	if (rows_begin_B==0){
		space_after_triangle_init_B = mat_B.rows.size;
	} else {
		space_after_triangle_init_B = (mat_B.rows.size > rows_begin_B) ? mat_B.rows.size - rows_begin_B : 0;
	}
	mat_A_triangle_columns			= MIN(mat_A_triangle_rows, space_after_triangle_init_B);
	/*printf("[%d | %d] Readj r(b%d,s%d) | mat_A triangle rows:%ld (F%d,L%d) | Space after triangle_init B %ld\n", 
		mpi_rank, stage, 
		mat_A.rows.begin, mat_A.rows.size, 
		space_after_triangle_init_B,
		mat_A_triangle_rows, triangle_first_pos, triangle_last_pos
	);*/
	#if DEBUG
	printf("[%d | %d] r(b%ld,s%ld) | mat_A triangle rows:%ld mat_A_triangle columns %ld\n", 
		mpi_rank, stage, 
		mat_A.rows.begin, mat_A.rows.size, 
		mat_A_triangle_rows,
		mat_A_triangle_columns
	);
	printf("[%d | %d] Readj r(b%ld,s%ld) | Space after triangle_init B %ld\n", 
		mpi_rank, stage, 
		mat_A.rows.begin, mat_A.rows.size, 
		space_after_triangle_init_B
	);
	#endif
	space_above_triangle		= MAX(0,mat_A.rows.size-mat_A_triangle_rows);
	
	rows_begin_B 				= space_left_to_triangle - mat_A.rows.begin;
	mat_B_rectangle_columns 	= MIN(mat_B.rows.begin+mat_B.rows.size, mat_B.columns.size);
	
	block_before_B 	= mat_B.columns.size*rows_begin_B;
	block_before_A 	= space_above_triangle*mat_A.columns.size+space_left_to_triangle;
	M = mat_A_triangle_columns; 	//mat_A_triangle_rows;			//B rows
	N = mat_B_rectangle_columns;	//B columns
	#if DEBUG
	printf("[%d | %d] mat_B rectangle columns:%ld\n", mpi_rank, stage, mat_B_rectangle_columns);
	printf("[%d | %d] Block Before A %ld Space Left %ld Space Above %ld\n", mpi_rank, stage, block_before_A, space_left_to_triangle,space_above_triangle);
	printf("[%d | %d] Block Before B %ld Begin Row %ld\n", mpi_rank, stage, block_before_B, rows_begin_B);
	fflush(stdout);	
	#endif 

	if (M > 0 && N > 0){
		#if DEBUG > 1
		printf("[%d | %d] M %ld N %ld\n",
		mpi_rank, stage, M, N);
		fflush(stdout);	
		matrix_rangePrintf_rank_stage(mat_A, "mataprodtri", mpi_rank, stage, 0);
		matrix_rangePrintf_rank_stage(mat_B, "matbprodtri", mpi_rank, stage, 0);
		printf("[%d | %d] block_before A %ld B %ld\n ", mpi_rank, stage, block_before_A, block_before_B);
		#endif 
		
		LDA=mat_A.columns.size;
		LDB=mat_B.columns.size;
		cblas_dtrmm ( 
			layout, sideA, triagA, TransA, unit,
			M, N,ALPHA, 
			mat_A.mat+block_before_A, LDA,
			mat_B.mat+block_before_B, LDB 
		);
	}
}



void mm_triangular_triangular_productTriangle_no_copy(
    MatrixRange mat_A, 
    MatrixRange mat_B,
    MatrixRange *mat_C,
    size_t block_before_C,
    unsigned int stage,
    unsigned int mpi_rank
){
    // Variables for MKL calls
    double ALPHA=1.0; 
	const CBLAS_LAYOUT 		layout = CblasRowMajor;
	const CBLAS_TRANSPOSE 	TransA = CblasNoTrans; 
	CBLAS_INDEX M, N;
	CBLAS_INDEX LDA, LDB; 
	const CBLAS_SIDE sideA 	= CblasLeft;		//A is left (A*B)
	const CBLAS_UPLO triagA = CblasLower;		//A is triangular L
	const CBLAS_DIAG unit 	= CblasNonUnit;		//A is not Identity (1)

	/*size_t mat_A_triangle_rows 		= (mat_A.rows.size == 0) ? 0 : mat_A.rows.size - 1;
	size_t mat_B_rectangle_columns 	= mat_B.columns.size; //mat_B.rows.begin+mat_B.rows.size;
	size_t block_before_B 			= 0; //mat_B.columns.size; --HERE
	size_t block_before_A 			= 0; //mat_A.columns.size+mat_A.rows.begin+1; --HERE
	block_before_A = mat_A.columns.size+mat_A.rows.begin+1;;//(mat_A.rows.size - mat_A_triangle_rows)*mat_A.columns.size+1; //--HERE
	
	M = mat_A_triangle_rows;		//A rows
	N = mat_B_rectangle_columns;	//B columns*/
	size_t mat_A_triangle_rows 		= (mat_A.rows.size == 0) ? 0 : mat_A.rows.size - 1;
	size_t mat_B_rectangle_columns 	= mat_B.rows.begin+mat_B.rows.size;
	size_t block_before_B 			= mat_B.columns.size;
	size_t block_before_A 			= mat_A.columns.size+mat_A.rows.begin+1; //mat_B.columns.begin+mat_A.columns.size+1; //mat_B.columns.size+mat_A.rows.begin+1;

	M = mat_A_triangle_rows;		//A rows
	N = mat_B_rectangle_columns;	//B columns

	if (M * N > 0){
		LDA=mat_A.columns.size;
		LDB=mat_B.columns.size;
		cblas_dtrmm ( 
			layout, sideA, triagA, TransA, unit,
			M, N,ALPHA, 
			mat_A.mat+block_before_A, LDA,
			mat_B.mat+block_before_B, LDB 
		);
	}
}



void mm_triangular_triangular_product(
    MatrixRange mat_A, 
    MatrixRange mat_B,
    MatrixRange *mat_C,
    size_t block_before_C,
    unsigned int stage,
    unsigned int mpi_rank
){
	/*
    // Triangular part
    // -- Triangular part * full rectangle
	// Only in last part*/

	unsigned long int id, jd, index1, index2;
	#if COPY_B
	MatrixRange mat_B_copy = matrix_range(TRIANGULAR_D);
	mat_B_copy.columns = mat_B.columns;
	mat_B_copy.rows = mat_B.rows;
	matrix_rangeMat_calloc_abort(&mat_B_copy, "Bcopy", MPI_COMM_WORLD, mpi_rank );
	#endif 
	#if DEBUG > 1
	matrix_rangePrintf_rank_stage(mat_A, "matA before rectangle", mpi_rank, stage, 0);
	matrix_rangePrintf_rank_stage(mat_B, "mat B Before Rectangle", mpi_rank, stage, 0);
	matrix_rangePrintf_rank_stage(*mat_C, "mat C Before Rectangle", mpi_rank, stage, block_before_C);
	#endif 
	mm_triangular_triangular_productRectangle(mat_A, mat_B, mat_C, block_before_C, stage, mpi_rank);
	
	#if DEBUG > 1
	matrix_rangePrintf_rank_stage(*mat_C, "mat C After Rectangle", mpi_rank, stage, block_before_C);
	#endif 
	if (stage == mpi_rank){                                       
		#if COPY_B
		for (unsigned long int i = 0; i < mat_B.rows.size; i++)     {
			for (unsigned long int j = 0; j < mat_B.columns.size; j++){
				index1=i*mat_B_copy.columns.size+j;
				mat_B_copy.mat[index1]=mat_B.mat[index1];
			}
		}
		
		mm_triangular_triangular_productTriangle(mat_A, mat_B_copy, stage, mpi_rank);
		#else 
				
		mm_triangular_triangular_productTriangle(mat_A, mat_B, stage, mpi_rank);
		#endif 
		#if DEBUG > 1
		printf("[%d | %d] mat_B(%ld,%ld)\n", mpi_rank, stage, mat_B.rows.size, mat_B.columns.size);
		fflush(stdout);
		#endif 
		unsigned long int rows_limit, columns_limit;
		rows_limit = MIN(mat_B.rows.size, mat_A.rows.size);		
		unsigned int rows_begin = (mat_A.rows.begin == 0) ? 0 : 1;
		//#pragma omp parallel for
		for(
			id=rows_begin; 
			id<rows_limit; 
			id++
		){
			columns_limit=mat_B.rows.begin+id;
			columns_limit = MIN(mat_C -> columns.size, columns_limit);
			for (
				jd=0; 
				jd<=columns_limit;
				jd++
			){
				index2=id*mat_B.columns.size+jd;
				index1=index2;
				#if COPY_B
				mat_C -> mat[index1] += mat_B_copy.mat[index2];
				#else 
				
					mat_C -> mat[index1] += mat_B.mat[index2];
				
				#endif 
			}
		}
		#if DEBUG > 1
		matrix_rangePrintf_rank_stage(*mat_C, "mat C After triangle", mpi_rank, stage, block_before_C);
		#endif
		fflush(stdout);
	} 

	
	
	#if COPY_B
	free(mat_B_copy.mat);
	#endif  

	#if DEBUG > 1
		printf("[%d | %d] AFTER STAGE PROD\n", mpi_rank, stage);
		fflush(stdout);
		matrix_rangePrintf_rank_stage(*mat_C, "Cstg", mpi_rank, stage, 0);
		fflush(stdout);
		matrix_rangePrintf_rank_stage(mat_A, "fullATP",mpi_rank, stage, 0);
	#endif 
}



void mm_triangular_rectangular_product(
    MatrixRange mat_A, 
    MatrixRange mat_B,
    MatrixRange *mat_C,
    size_t block_before_C,
    unsigned int stage,
    unsigned int mpi_rank
){	
	int mpi_procs;
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );//Total num of procs
	MatrixRange mat_B_copy = matrix_range(RECTANGULAR);
	int f, c;
	unsigned long int index1, index2;
	if (mpi_procs > 1){
		mat_B_copy.columns = mat_B.columns;
		mat_B_copy.rows = mat_A.rows;
		matrix_rangeMat_calloc_abort(&mat_B_copy, "Bcopy", MPI_COMM_WORLD, mpi_rank);
		for (f=0; f < mat_A.rows.size; f++){
			for (c=0; c < mat_B.columns.size; c++){
				index1=f*mat_B.columns.size+c;
				index2=(f+mat_A.rows.begin)*mat_B.columns.size+c;
				mat_B_copy.mat[index1]=mat_B.mat[index2];
			}
		}

	}
	// Sizes [horizontal]x[vertical]
	// Rectangle [0:current_range.begin]x[0:current_range.size],
	// Triangle	[current_range.begin+1,current_range.begin+current_range.size]x[1:current_range.size]
	
	mm_triangular_rectangular_productRECTANGLE(mat_A, mat_B, mat_C, block_before_C, stage, mpi_rank);
	
	//mm_triangular_triangular_productRectangle(mat_A, mat_B, mat_C, block_before_C, stage, mpi_rank);
	#if DEBUG
	printf("[%d | %d] AFTER RECTANGULAR PART\n", mpi_rank, stage);
	#endif
	#if DEBUG  > 1
	matrix_rangePrintf_rank_stage(*mat_C, "RLC", mpi_rank, stage, 0);
	fflush(stdout);
	#endif
    // Triangular part
    // -- Triangular part * full rectangle
	// Only in last part
	//-- Probar a sacarlo fuera del for
	
	mm_triangular_rectangular_productTRIANGLE(mat_A, mat_B, mat_C, block_before_C, stage, mpi_rank);
	
	if (mpi_procs > 1){
		for (f=0; f < mat_A.rows.size; f++){
			for (c=0; c < mat_B.columns.size; c++){
				index1=f*mat_B.columns.size+c;
				index2=(f+mat_A.rows.begin)*mat_B.columns.size+c;
				mat_B.mat[index2]=mat_B_copy.mat[index1];
			}
		}
		free(mat_B_copy.mat);
	}

	//mm_triangular_triangular_productTriangle(mat_A, mat_B, stage, mpi_rank);
	#if DEBUG > 1
		printf("[%d | %d] AFTER STAGE PROD\n", mpi_rank, stage);
		fflush(stdout);
		matrix_rangePrintf_rank_stage(*mat_C, "Cstg", mpi_rank, stage, 0);
		fflush(stdout);
	#endif 
}

