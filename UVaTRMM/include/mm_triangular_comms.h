#ifndef MM_TRIANG_COMMS_H
    #define MM_TRIANG_COMMS_H
    #include "buff.h"
    #include "matrix_range.h"
    #define SPLITTED 1
    /**
    * @private Structure MPI_Comms
    * Simplifies the use of buffers
    * * Buffer is buffer itself
    * * Size is the number of elements in the buffer
    * */
    typedef struct {
    	MPI_Comm *buffer;
    	int size;
    } MPI_Comms;

    /***
     * @private @section Matrix setting functions
     * */
    
    ///@private Sets up mat_A_current
    void matrix_rangeSet_mat_A_current(MatrixRange *mat_A_current,MatrixRange mat_A_local,MatrixRange mat_A_remote, int stage, int mpi_procs, int mpi_rank, Map_function f_mapping, Ranges *mat_A_current_rows);
    /***
     * @private @section Send-receive functions
     * */
    void comm_func_bcast_sarteco(MatrixRange mat_A_local,MatrixRange *mat_A_global,MatrixRange *mat_A_current,int stage,Requests *requests,Requests *requests2,int mpi_rank,int mpi_procs,Map_function f_mapping,Type_function create_type,int rows,MPI_Datatype *type1,MPI_Datatype *type2,Buff *triang_buff,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism);
    void comm_func_bcast(MatrixRange mat_A, MatrixRange *current_A,int stage, Requests *broad_requests,  Requests *broad_requests2,int mpi_rank, int mpi_procs,Map_function f_mapping,int rows,MPI_Datatype *type_send1,MPI_Datatype *type_send2,Buff *triang_buff_send,MPI_Datatype *type_recv1,MPI_Datatype *type_recv2,Buff *triang_buff_recv,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism);
    void comm_func_pipeline_recv(MatrixRange *mat_A_remote, int stage,Requests *requests, int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatype *type_recv1,MPI_Datatype *type_recv2,Buff *triang_buff_recv,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,int comm_from);
    void comm_func_pipeline_send(MatrixRange mat_A,int stage,Requests *requests,int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatype *type_send1,MPI_Datatype *type_send2,Buff *triang_buff_send,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,int comm_to);
    void comm_func_pipeline_recvArray(MatrixRange *mat_A_remote, int stage,Requests *requests, int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatypes *type_recv1,MPI_Datatypes *type_recv2,Buff *triang_buff_recv,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,int comm_from, unsigned int quotient, MPI_Comm comm);
    void comm_func_pipeline_sendArray(MatrixRange mat_A,int stage,Requests *requests,int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatypes *type_send1,MPI_Datatypes *type_send2,Buff *triang_buff_send,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,int comm_to,unsigned int quotient, MPI_Comm comm);
    void comm_func_pipeline_recv_(MatrixRange *mat_A_remote, int stage,Requests *requests, int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatype *type_recv1,MPI_Datatype *type_recv2,Buff *triang_buff_recv,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,int comm_from);
    void comm_func_pipeline_send_(MatrixRange mat_A,int stage,Requests *requests,int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatype *type_send1,MPI_Datatype *type_send2,Buff *triang_buff_send,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,int comm_to);

    ///p2p = owner skew
    void comm_func_p2p_recv(int comm_from,MatrixRange mat_A, MatrixRange *current_A,int stage, Requests *requests, int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatype *type_recv1,MPI_Datatype *type_recv2,Buff *triang_buff_recv,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism);
    void comm_func_p2p_send(int comm_to,MatrixRange mat_A,int stage, Requests *requests, int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatype *type_send1, MPI_Datatype *type_send2, Buff *triang_buff_send, Params_comm_shape comm_shape, Params_comm_mechanism comm_mechanism);
    void comm_func_p2p_recvArray(int comm_from,MatrixRange mat_A, MatrixRange *current_A,int stage, Requests *requests, int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatypes *type_recv1,MPI_Datatypes *type_recv2,Buff *triang_buff_recv,Params_comm_shape comm_shape,Params_comm_mechanism comm_mechanism,unsigned int quotient, MPI_Comm comm );
    void comm_func_p2p_sendArray(int comm_to,MatrixRange mat_A,int stage, Requests *requests, int mpi_rank,int mpi_procs,Map_function f_mapping,int rows,MPI_Datatypes *type_send1, MPI_Datatypes *type_send2, Buff *triang_buff_send, Params_comm_shape comm_shape, Params_comm_mechanism comm_mechanism,unsigned int quotient, MPI_Comm comm );
    void mpi_broadcast_splitted(double *buffer,unsigned long int count,int root, MPI_Comm comm, Requests *requests, int mpi_rank, MPI_Datatype datatype);
    void mpi_broadcast_splitted_quotient(double *buffer,unsigned long int count,int root, MPI_Comm comm, Requests *requests, int mpi_rank,MPI_Datatype datatype,unsigned int quotient);
    void mpi_isend_splitted(double *buffer, unsigned long int count, int comm_to, MPI_Comm comm, Requests *requests, int mpi_rank,MPI_Datatype datatype,int request_id,int max_requests_id, int tag, int initialize);
    void mpi_irecv_splitted(double *buffer, unsigned long int count, int comm_from, MPI_Comm comm, Requests *requests,int mpi_rank,MPI_Datatype datatype,int request_id,int max_requests_id, int tag);
    void mpi_isend_array(double *buffer,int comm_to, MPI_Comm comm, Requests *requests, int mpi_rank,MPI_Datatypes datatypes,int request_id,int max_requests_id, int tag, int initialize);
    void mpi_irecv_array(double *buffer,int comm_from, MPI_Comm comm, Requests *requests,int mpi_rank,MPI_Datatypes datatypes,int request_id,int max_requests_id, int tag);
    void mpi_wait_splitted(Requests *requests, Statuses *status);

    
    void comm_func_bcast_quotient(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_requests,
	Requests *broad_requests2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatype *type_send1,
	MPI_Datatype *type_send2,
	Buff *triang_buff_send,
	MPI_Datatype *type_recv1,
	MPI_Datatype *type_recv2,
	Buff *triang_buff_recv,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism, unsigned int quotient
);
    
    /***
     * @section Buff communication functions
     * */
    void buffIbcast(Buff *buffer, int proc, Requests *requests, MPI_Comm comm, int mpi_rank);
        /***
     * @section Buff communication functions
     * */
    void buffIbcast_quotient(Buff *buffer, int proc, Requests *requests, MPI_Comm comm, int mpi_rank, unsigned int quotient);
    /***
     * @section MatrixRange communication functions
     * */
    void matrix_rangeIbcast(MatrixRange *mat, int proc, Requests *requests, MPI_Comm comm, unsigned long int block_before, int mpi_rank);
    void matrix_rangeISend(MatrixRange *mat, int comm_from, int comm_to, Requests *requests, MPI_Comm comm, unsigned long int block_before, unsigned int req_id);
    void matrix_rangeIRecv(MatrixRange *mat, int comm_from, int comm_to, Requests *requests, MPI_Comm comm, unsigned long int block_before, unsigned int req_id);
    /***
     * @private
     * @section Communication functions
     * Porrefinar 
     * **/
    void mm_triangular_triangularMpi_commSplit(MPI_Comms *comms, MPI_Comm main_comm, unsigned int mpi_rank, unsigned int mpi_procs);
    void mpi_datatypeBcast_sarteco(
	double *buff, 		//Array 
	MPI_Datatype *type, //Type
	int proc, 			//Process sending/receiving
	Requests *requests,	//MPI_Request buffer
	MPI_Comm comm, 		//Communicator
	size_t block_before,
	int mpi_rank);

    void mpi_datatypeBcast_sarteco_quotient(
	double *buff, 		//Array 
	MPI_Datatype *type, //Type
	int proc, 			//Process sending/receiving
	Requests *requests,	//MPI_Request buffer
	MPI_Comm comm, 		//Communicator
	size_t block_before,
	int mpi_rank, unsigned int quotient);

    void mpi_datatypeBcast_sarteco_quotientArray(
	double *buff, 		//Array 
	MPI_Datatypes *type, //Type
	int proc, 			//Process sending/receiving
	Requests *requests,	//MPI_Request buffer
	MPI_Comm comm, 		//Communicator
	size_t block_before,
	int mpi_rank
);
void mpi_broadcast_splitted_quotientArray(
    double *buffer,
    int root, 
    MPI_Comm comm, 
    Requests *requests,//Esto lo pondre como MPI_Requests en el "oficial"
    int mpi_rank,
	MPI_Datatypes datatypes
);

void comm_func_bcast_quotientArray(
	MatrixRange mat_A, 
	MatrixRange *current_A,
	int stage, 
	Requests *broad_requests,
	Requests *broad_requests2,
	int mpi_rank,
	int mpi_procs,
	Map_function f_mapping,
	int rows,
	MPI_Datatypes *type_send1,
	MPI_Datatypes *type_send2,
	Buff *triang_buff_send,
	MPI_Datatypes *type_recv1,
	MPI_Datatypes *type_recv2,
	Buff *triang_buff_recv,
	Params_comm_shape comm_shape,
	Params_comm_mechanism comm_mechanism, 
    unsigned int quotient, 
	MPI_Comm comm 
);
#endif
