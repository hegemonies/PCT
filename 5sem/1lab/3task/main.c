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

	char *rbuf;
	if (rank == 0) {
		rbuf = malloc(sizeof(char) * sizeMessage * max_rank);
	}
	double t;
	for (int j = 0; j < 2; j++) {
		t = 0;
		t -= MPI_Wtime();
		if (rank == 0) {
			for (int i = 1; i < max_rank; i++) {
				MPI_Recv(rbuf + (sizeMessage * i), sizeMessage, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			MPI_Sendrecv(sbuf, sizeMessage, MPI_CHAR, 0, 0, rbuf, sizeMessage, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			t += MPI_Wtime();
			printf("[%d] Task %d is done  Elapsed time: %lf\n", rank, j, t);
		} else {
			MPI_Send(sbuf, sizeMessage, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
			t += MPI_Wtime();
			printf("[%d] Task %d is done  Elapsed time: %lf\n", rank, j, t);
		}

		sizeMessage *= 1024;
		sbuf = (char *)realloc(sbuf, sizeMessage);
		if (rank == 0) {
			rbuf = (char *)realloc(rbuf, sizeMessage * max_rank);
		}
	}

	MPI_Finalize();
	return 0;
}