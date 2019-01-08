# matrix multiplication with openmp intended as a Caliper test

CC=icc

INC=
LIB=-qopenmp

X_SIZE=10000
Y_SIZE=20000

VARS=-DX_SIZE=${X_SIZE} 
VARS+=-DY_SIZE=${Y_SIZE} 

CALI_INC=-I${CALIPER_DIR}/include
CALI_LIB=-L${CALIPER_DIR}/lib64 -lcaliper

all: stencil

stencil: stencil.c
	${CC} -g -o test_stencil ${INC} stencil.c ${LIB} ${VARS}

clean:
	rm -f *.out test_stencil *.o *.cali *.json
	# rm -rf MULTI__*

