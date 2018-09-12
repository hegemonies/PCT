#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);

	int rank;
	int max_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &max_rank);

	int size_of_message = 1;
	char *message = malloc(sizeof(char) * size_of_message);

	int count = 0;

	while (count < 3) {
		if (rank == 0) {
			*message = 'f';
			MPI_Send(message, size_of_message, MPI_CHAR, rank + 1, 0, MPI_STATUS_IGNORE);
			char *buf = malloc(sizeof(char) * size_of_message);
			MPI_Recv(buf, size_of_message, MPI_CHAR, max_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		} else if (rank > 0-&& rank < ) {
			MPI_Recv(message, size_of_message, MPI_CHAR, rank - 1)
		}
	}
	
	
	MPI_Finalize();
	return 0;
}