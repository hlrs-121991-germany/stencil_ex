# matrix multiplication with openmp intended as a Caliper test

CC=icc

INC=
LIB=-qopenmp


X_SIZE=5000
Y_SIZE=2000
# X_SIZE=15
# Y_SIZE=10
STENCIL_TYPE=0

VARS=-DX_SIZE=${X_SIZE} 
VARS+=-DY_SIZE=${Y_SIZE} 
VARS+=-DSTENCIL_TYPE=${STENCIL_TYPE} 

# VARS+=-DDO_IO=1
VARS+=-DFILE_NAME=\"openmp_mesh.out\"

CALI_INC=-I${CALIPER_DIR}/include
CALI_LIB=-L${CALIPER_DIR}/lib64 -lcaliper

all: stencil

stencil: stencil.c
	${CC} -g -o test_stencil ${INC} stencil.c stencil_patterns.c ${LIB} ${VARS}

clean:
	rm -f *.out test_stencil *.o *.cali *.json
	# rm -rf MULTI__*

clena: clean

