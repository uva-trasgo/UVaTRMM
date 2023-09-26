#ifndef MAP_FUNCS_H
    #define MAP_FUNCS_H
    /***
    * Matrices maping function headers
    * **/
    #include "mkl_range.h"
    /***
    * @private Direction
    * Direccion de particion (bloque vertical/bloque horizontals)
    * */
    typedef enum {
        VERTICAL,
        HORIZONTAL
    } Direction;
    /*
    Mapping functions
    */
    extern Range map_regular( unsigned long int num_rows, unsigned int num_procs, unsigned int ind_proc, int mpi_rank );
    extern Range map_regular_first_double( unsigned long int num_rows, unsigned int num_procs, unsigned int ind_proc, int mpi_rank );
    extern Range map_balanced( unsigned long int num_rows, unsigned int num_procs, unsigned int ind_proc, int mpi_rank );
#endif 

