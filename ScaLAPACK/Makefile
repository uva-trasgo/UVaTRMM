C_FLAGS=-std=c99 -O3
LIBS=-lm
LIBS_MKL=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread # MKL Seq
LIBS_MKL_MTH=-lmkl_core -lmkl_rt -lpthread -lmkl_gnu_thread -lgomp # MKL mth
LIBS_SCA=-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

all: sca_trmm_seq sca_trmm_mth sca_gemm_seq sca_gemm_mth

sca_trmm_seq:
	mpicc scalapack_sample_trmm.c ${C_FLAGS} ${LIBS_MKL} ${LIBS_SCA} ${LIBS} -o scalapack_sample_trmm_seq

sca_trmm_mth:
	mpicc scalapack_sample_trmm.c ${C_FLAGS} ${LIBS_MKL_MTH} ${LIBS_SCA} ${LIBS} -o scalapack_sample_trmm_mth

sca_gemm_seq:
	mpicc scalapack_sample_gemm.c ${C_FLAGS} ${LIBS_MKL} ${LIBS_SCA} ${LIBS} -o scalapack_sample_gemm_seq

sca_gemm_mth:
	mpicc scalapack_sample_gemm.c ${C_FLAGS} ${LIBS_MKL_MTH} ${LIBS_SCA} ${LIBS} -o scalapack_sample_gemm_mth
