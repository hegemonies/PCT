#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);

	int rank;
	int max_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &max_rank);

	int sizeMessage = 1024;

	char *sbuf = malloc(sizeof(char) * sizeMessage);
	char *rbuf = malloc(sizeof(char) * sizeMessage);

	double t;
	t -= MPI_Wtime();
	
	MPI_Request requests[max_rank];

	for (int i = 0; i < max_rank; i++) {
		MPI_Isend(sbuf, sizeMessage, MPI_CHAR, i, 0, MPI_COMM_WORLD, &requests[i]);
	}

	MPI_status status;
	for (int i = 0; i < max_rank; i++) {
		MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_IRecv(rbuf, sizeMessage, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &requests[status.MPI_SOURCE]);
		printf("[%d] recv from [%d]\n", rank, status.MPI_SOURCE);
	}

	t += MPI_Wtime();
	printf("[%d] Elapsed time is %lf", rank, t);
	
	MPI_Waitall(max_rank * 2, requests, MPI_STATUS_IGNORE);

	MPI_Finalize();
	return 0;
}