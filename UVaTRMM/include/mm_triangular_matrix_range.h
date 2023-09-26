#ifndef MM_TRIANGULAR_MATRIX_RANGE_H
    #define MM_TRIANGULAR_MATRIX_RANGE_H
    #include "matrix_range.h"
    #include "mm_triangular_params.h"
    #include "buff.h"
    /****
     * @private @section Matrix generation functions
     * */
    extern void mat_identity(Range rows, Range columns, double *mat);
    extern void mat_consecutive_values(Range rows, Range columns, double *mat, int triangular, int global_columns);
    
    extern void matrix_rangePoint_to(MatrixRange *Matrix, MatrixRange Remote);


    //--> REVISAR SI SE USAN EN ALGUN LADO O SE PUEDEN BORRAR
    extern void set_mat_A(Range rows, int columns, double *mat_A, Params_init init_scheme, int global_columns);
    extern void set_mat_B(int rows, Range columns, double *mat_B, Params_init init_scheme, int global_columns);
    extern void set_mat_A_range(Range rows, Range columns, double *mat_A, Params_init init_scheme, int global_columns);
    extern void set_mat_B_range(Range rows, Range columns, double *mat_B, Params_init init_scheme, int global_columns);

    void buff__matrix_range(Buff buffer, MatrixRange *mat);

    /*
        Matrix functions
    */
    void matrix_traspose(double *mat_A, double *mat_B,int rows_A,int columns_A,int block_before_A, int block_Before_B);
    _Bool matrix_equal(Range columns, Range rows,double *mat_A, double *mat_B,int block_before_A,int block_before_B);
    MatrixRange matrix_rangeCopy_null(MatrixRange orig, const char *name, MPI_Comm comm, int mpi_rank);
    MatrixRange matrix_rangeCalloc(size_t rows, size_t columns, const char *name, MPI_Comm comm, int mpi_rank);
    void matrix_rangeSetA(MatrixRange *mat, Params_init init_scheme, int global_columns);
    void matrix_rangeSetB(MatrixRange *mat, Params_init init_scheme, int global_columns);
#endif
