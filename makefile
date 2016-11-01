EXECS=mpi_psrs
MPICC?=mpicc

all: ${EXECS}

mpi_psrs: mpi_psrs.c
	${MPICC} -o mpi_psrs mpi_psrs.c

clean:
	rm -f ${EXECS}