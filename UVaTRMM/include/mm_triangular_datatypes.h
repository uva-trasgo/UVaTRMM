#ifndef MM_DATATYPES_H
    #define MM_DATATYPES_H
    #include "mkl_range.h"
    #include "mm_triangular_params.h"
    
    /***
     * @section STRUCTURE DEFINITIONS
     * */
    /**
     * @private struct MPI_Datatypes
     * arr: buffer of MPI_Datatype
     * total: number of types colected
     * */
    typedef struct{
        size_t total;
        MPI_Datatype *arr;
        size_t *sizes;
    } MPI_Datatypes;
        /**
    * @private Structure Buff
    * Simplifies the use of buffers
    * * Buffer is buffer itself
    * * Size is the number of elements in the buffer
    * */
    typedef struct {
        /*@{*/
    	size_t *buffer; /**< Double array */
    	size_t size;       /**< Current buffer size */
        /*@}*/
    } Sizes;
     /**
     * @private Sizes
     * @brief Struct Sizes builder
     * @param[in] size      Size of the buff we want to create
     * @param[in] mpi_rank  Current proccessor identificator
     * @return Buff of size "size" or aborts in case of error
     * */
    Sizes sizes(int size, int mpi_rank, MPI_Comm comm, char *name);
    
    /***
     * @section TYPE CREATION FUNCTIONS
     * */
    /**
     * @private Shapes type generation prototype
     * @param[in] range     is the range of the rows (or columns) we can build the shape at
     * @param[in] integer   is the total number of rows (or columns) in the matirx
     * */
    typedef	MPI_Datatype (*Type_function)(Range, unsigned long int);
    
    /**
     * @private create_type_box
     * Creates a type that represents the bounding box of rectangle + triangle of matrix's range block
     * */
    extern MPI_Datatype create_type_bbox( Range range, unsigned long int columns );
    /**
     * @private create_type_box
     * Creates a type that represents the rectangle of matrix's range block
     * */
    extern MPI_Datatype create_type_rectangle( Range range, unsigned long int columns );
    /**
     * @private create_type_box
     * Creates a type that represents the triangle of matrix's range block
     * */
    extern MPI_Datatype create_type_triangle( Range range, unsigned long int columns );
    /**
     * @private create_type_box
     * Creates a type that combines both rectangle and triangle of matrix's range block
     * */
    extern MPI_Datatype create_type_combined( Range range, unsigned long int columns );
    /**
     * @private create_type_box
     * Creates a type that represent the complete part down-diagonal of matrix's range block
     * */
    extern MPI_Datatype create_type_trapezoid( Range range, unsigned long int columns );
    /**
     * @private Not in use
     * */
    extern MPI_Datatype create_type_trapezoid_bilateral( Range range, unsigned long int columns);
    /***
     * @section MPI_DATATYPES FUNCTIONS
     * */
    extern MPI_Datatypes mpi_datatypesCalloc(size_t total, char *name, MPI_Comm comm, int mpi_rank);
    void mpi_datatypesFree(MPI_Datatypes *types);
    
    /**
     * @private Create type bbox array
     * Creates an array of types equivalents to bbox in order to not to overflow MPI's send-receives
     * */
    MPI_Datatypes create_type_bboxArray( Range range, unsigned long int columns, size_t block_limit_quotient_denominator, int mpi_rank, MPI_Comm comm);
    MPI_Datatypes create_type_rectangleArray( Range range, unsigned long int columns, int mpi_rank, MPI_Comm comm, Sizes heights, Sizes displacements);
    MPI_Datatypes create_type_triangleArray( Range range, unsigned long int columns, int mpi_rank, MPI_Comm comm, Sizes heights, Sizes displacements, int stage);
    MPI_Datatypes create_type_combinedArray( Range range, unsigned long int columns, size_t block_limit_quotient_denominator, int mpi_rank, MPI_Comm comm, int stage, Sizes heights, Sizes displacements);

#endif 

