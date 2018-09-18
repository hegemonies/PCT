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

	char *rbuf;

	for (int j = 0; j < 2; j++) {
		rbuf = calloc(sizeMessage, sizeof(char));
		if (rank == 0) {
			double t = 0;
			t -= MPI_Wtime();
			char *sbuf = calloc(sizeMessage, sizeof(char));
			int st = (sizeMessage / max_rank) * rank;
			int fn = st + (sizeMessage / max_rank);

			for (int i = st; i < fn; i++) {
				sbuf[i] = 128;
			}

			for (int i = 1; i < max_rank; i++) {
				MPI_Send(sbuf, sizeMessage, MPI_CHAR, i, 0, MPI_COMM_WORLD);
			}
			MPI_Status status;
			MPI_Sendrecv(sbuf, sizeMessage, MPI_CHAR, 0, 0, rbuf, sizeMessage, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			int count;
			MPI_Get_count(&status, MPI_CHAR, &count);
			
			printf("[%d] recv %d B\n", rank, count);
			t += MPI_Wtime();
			if (rank == 0) {
				printf("[%d] Task %d is done\nElapsed time: %lf\n", rank, j, t);
			}
		} else {
			MPI_Status status;
			MPI_Recv(rbuf, sizeMessage, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
			int count;
			MPI_Get_count(&status, MPI_CHAR, &count);
			
			printf("[%d] recv %d B\n", rank, count);

			rbuf = (char *) realloc(rbuf, sizeMessage);
		}
		
		sizeMessage *= 1024;
		free(rbuf);
	}

	MPI_Finalize();
	return 0;
}