########################
# PATHS AND EXTENSIONS #
########################
MAIN_PATH=./src/
EXEC_PATH=./exec/
EXEC_EXTENSION=.x

#Main versions
T=_triangular
VERSION_LU_FULL=mm$(T)

###########
# GENERAL #
###########
OBJS=./exec/*.x

#####################
# COMPILATION FLAGS #
#####################
#Compilers
GCC=mpicc

#Mandatory
FLAGS_MANDATORY_COMMON=-O3 -fopenmp
FLAGS_MANDATORY_GCC=-std=c99

#Entire lines
FLAGS_GCC= $(FLAGS_MANDATORY_COMMON) $(FLAGS_MANDATORY_GCC)
FLAGS_AOCC=$(FLAGS_MANDATORY_COMMON) $(FLAGS_MANDATORY_GCC)

######################
# EXTERNAL LIBRARIES #
######################
LIBS=-lm
LIBS_MKL=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread # MKL Seq
LIBS_MKL_MTH=-lmkl_core -lmkl_rt -lpthread -lmkl_gnu_thread -lgomp # MKL mth

HEADER_MKL=-I/LUSTRE/SOFT/calendula2/intel/oneapi/mkl/2021.3.0/include
LIBRARY_MKL=-L/LUSTRE/SOFT/calendula2/intel/oneapi/mkl/2021.3.0/lib/intel64

##################
# INTERNAL FILES #
##################

#Internal libraries
HEADER_COMMON=-I./ -I./include -I./src/

#Auxiliar files needed to comppile before
SRC_PATH=./src/
COMMON=$(SRC_PATH)mm_triangular_common.c
BROADCAST=$(SRC_PATH)comm/mm_triangular_broadcast.c
SARTECO=$(SRC_PATH)comm/mm_triangular_broadcast_sarteco.c
P2P=$(SRC_PATH)comm/mm_triangular_p2p.c
PIPELINE=$(SRC_PATH)comm/mm_triangular_pipeline.c
COMMS=$(SRC_PATH)comm/mm_comm_split.c
PRODUCT=$(SRC_PATH)mm_triangular_rectangular_product.c
MAPPING=$(SRC_PATH)mm_mapping.c

MAIN_AUXILIAR=$(COMMON) $(MAPPING) $(BROADCAST) $(SARTECO) $(P2P) $(PIPELINE) $(PRODUCT) $(COMMS) 

todo: todo_aux todo

todo_aux: print_mat checks
todo: todo_sequential todo_multithreaded

todo_sequential: sequential 

todo_multithreaded: multithreaded 

clean: rm -f $(OBJS)

#Build checking executables 
checks: check check_triang check_triang_triang

#Build main executables
sequential: $(VERSION_LU_FULL)
multithreaded: $(VERSION_LU_FULL)_mth

#####################
# COMPILATION LINES #
#####################
comp_line_common=$(HEADER_MKL) $(HEADER_COMMON) $(MAIN_AUXILIAR) $< $(LIBRARY_MKL) -o $(EXEC_PATH)$@$(EXEC_EXTENSION) $(LIBS) 
repo_prefix=$(HEADER_MKL) $(HEADER_COMMON) $(MAIN_AUXILIAR) $(LIBRARY_MKL)

comp_line_gcc	= $(GCC) $(FLAGS_GCC) $(comp_line_common)

# MAIN PROGRAMS

## Sequential
$(VERSION_LU_FULL): $(MAIN_PATH)$(VERSION_LU_FULL).c 
	$(comp_line_gcc) $(LIBS_MKL)

## Multithread
$(VERSION_LU_FULL)_mth: $(MAIN_PATH)$(VERSION_LU_FULL).c 
	$(comp_line_gcc) $(LIBS_MKL_MTH)

######################
# AUXILIAR FUNCTIONS #
######################

#Printing matrices
print_mat: $(SRC_PATH)print_mat.c
	gcc $(FLAGS_GCC) $< -o $(EXEC_PATH)$@$(EXEC_EXTENSION)

#Checking results
check_triang: ./src/values_init_check_triang.c
	$(GCC) $(FLAGS_GCC) $(HEADER_MKL)  $< $(LIBRARY_MKL)  -o $(EXEC_PATH)$@$(EXEC_EXTENSION) $(LIBS) $(LIBS_MKL)