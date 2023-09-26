#ifndef MM_TRIANGULAR_PARAMS
	#define MM_TRIANGULAR_PARAMS
	#include "matrix_range.h"
	#include "mm_triangular_datatypes.h"
	#include "buff.h"
	/**************************************
 	* ENUM TYPES FOR SELECTOR PARAMETERS *
 	* ************************************/ 
	/***
 	* @private Params_init
 	* Matrices initialization modes
 	* **/
	typedef enum {
		INIT_VALUES,/// both A and B initialized using consecutive values
		INIT_A_ID, 	/// A identity, B consecutive values
		INIT_B_ID	/// A consecutive values, B identity
	} Params_init;
	/***
 	* @private Params_output
 	* Output generation flag
 	* **/
	typedef enum {
		OUTPUT_YES, OUTPUT_NO
	} Params_output;

	/***
 	* @private Params_map
 	* Matrix mapping options
 	* */
	typedef enum {
		MAP_REGULAR, 		///Divide matrix into regular blocks
		MAP_REGULAR_DOUBLE, ///Divide matrix into regular blocks, except first one (double-sized)
		MAP_BALANCED		///Divides matrix into blocks with as equal as posible number of non-zero elements
	} Params_map;
	/***
 	* @private Params_comm_scheme
 	* **/
	typedef enum {
		COMM_OWNER, 	///Owner sends to all
		COMM_PREVIOUS	///Previous sends to next
	} Params_comm_scheme;
	/***
 	* @private Params_comm_shape
 	* **/
	typedef enum {
		SHAPE_FULL, 	//Full rows
		SHAPE_BOXES, 	//Smallest rectangle with all non-zero elements
		SHAPE_TRIANG_1, 
		SHAPE_TRIANG_2, 
		SHAPE_TRAPEZOID, 
		SHAPE_TRAPEZOID_BILATERAL
	} Params_comm_shape;
	typedef enum {
		MECH_BUFFER, MECH_TYPE
	} Params_comm_mechanism;
	typedef enum {
		FUNC_P2P, FUNC_BROADCAST, FUNC_BROADCAST_SARTECO
	} Params_comm_func;

	typedef	Range (*Map_function)(unsigned long int, unsigned int, unsigned int, int);


	/*
	Initial parameters setting up functions
	*/
	extern void set_initial_scheme(int argc, char *argv[],int mpi_rank,Params_init *init_scheme,int *do_check);
	extern void set_output_result(int argc, char *argv[],int mpi_rank,Params_output *output);
	extern void set_mapping_function(int argc, char *argv[],int mpi_rank,Map_function *f_mapping);
	extern void set_communication_scheme(int arg, char *argv[],int mpi_rank,Params_comm_scheme *comm_scheme,Params_comm_func *comm_func);
	extern void set_communication_shape(int arg, char *argv[],int mpi_rank,Params_comm_shape *comm_shape,Params_comm_mechanism *comm_mechanism, Type_function *create_type);
    /*
        Type set-free functions
    */
   	//-- Originals
    void free_types_recv(Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,MPI_Datatype *type_recv1, MPI_Datatype *type_recv2, MatrixRange *current_A, Buff *triang_buffer_recv);
    Buff set_type_send(int mpi_rank, MPI_Datatype *type_send1,MPI_Datatype *type_send2,Buff *triang_buff_send, MatrixRange *mat_A,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism, Params_comm_func comm_func);
   	void set_type_recv(int mpi_rank, Range remote_rows, int columns, MPI_Datatype *type_recv1,MPI_Datatype *type_recv2,Buff *triang_buff_recv,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism);
	//-- Array
	void free_types_recvArray(Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,MPI_Datatypes *type_recv1, MPI_Datatypes *type_recv2, MatrixRange *current_A,Buff *triang_buffer_recv);
	Buff set_type_sendArray(int mpi_rank, MPI_Datatypes *type_send1, MPI_Datatypes *type_send2,Buff *triang_buff_send, MatrixRange *mat_A,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,Params_comm_func comm_func,size_t block_limit_quotient_denominator, MPI_Comm comm, int stage);
	void set_type_recvArray(int mpi_rank, Range remote_rows, int columns, MPI_Datatypes *type_recv1,MPI_Datatypes *type_recv2,Buff *triang_buff_recv,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism, size_t block_limit_quotient_denominator, MPI_Comm comm, int stage);
#endif //MM_TRIANGULAR_PARAMS
