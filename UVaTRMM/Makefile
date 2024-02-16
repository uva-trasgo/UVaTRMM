#####################
# PATHS & NAMES     #
#####################
MAIN_PATH=./src/
EXEC_PATH=./exec/
EXEC_EXTENSION=.x

#####################
# COMPILATION FLAGS #
#####################
GCC=mpicc
FLAGS_GCC=-O3 -fopenmp -std=c99
DNMI=-DMININUMS=1

######################
# EXTERNAL LIBRARIES #
######################
LIBS=-lm
LIBS_MKL=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread # MKL Seq

##################
# INTERNAL FILES #
##################
HEADER_COMMON=-I./ -I./include -I./src/
SRC_PATH=./src/
COMMON=$(SRC_PATH)mm_triangular_common.c
COMMS=$(SRC_PATH)mm_comm_split.c
MAPPING=$(SRC_PATH)mm_mapping.c

all: uvatrmm

clean: 
        rm -f $(EXEC_PATH)*$(EXEC_EXTENSION)

#####################
# COMPILATION LINES #
#####################
uvatrmm: $(MAIN_PATH)mm_triangular.c
        $(GCC) $(FLAGS_GCC) $(HEADER_MKL) $(HEADER_COMMON) $(COMMON) $(MAPPING) $(COMMS) $< $(LIBRARY_MKL) -o $(EXEC_PATH)uvatrmm$(EXEC_EXTENSION) $(LIBS) $(DNMI) $(LIBS_MKL)
