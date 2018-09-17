#include <string.h>
#include <stdio.h>
#include <stdlib.h>
// #include <time.h>
#include <mpi.h>

char *printArrayBytes(unsigned char *arr, int size)
{
	// for (int i = 0; i < size; i++) {
	// 	printf("%d ", arr[i]);
	// }
	// printf(" \b");

	char *buf = calloc(0, sizeof(char) * size * 12);
	for (int i = 0; i < size; i++) {
		if (arr[i] == 0) {
			strcat(buf, "0000 0000 ");
			continue;
		}
		for (int j = 7; j >= 0; j--) {
			if ((arr[i] & (0x1 << j)) != 0) {
				strcat(buf, "1");
				// printf("1");
			} else {
				strcat(buf, "0");
				// printf("0");
			}
			if (j == 4) {
				strcat(buf, " ");
				// printf(" ");
			}
		}
		// strcat(buf, "\n");
		// printf(" ");
		strcat(buf, " ");
	}
	// printf("\n");
	// printf("%s\n", buf);
	// strcat(buf, "\0");
	return buf;
}

int main(int argc, char **argv) {
	// srand(time(NULL));
	MPI_Init(&argc, &argv);

	int rank;
	int max_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &max_rank);

	// FILE *file;
	// char rank_str[2];
	// rank_str[0] = (char)rank;
	// char *name_of_file = strcat("log_rank", rank_str);
	// printf("test = %s\n", name_of_file);
	// remove(name_of_file);
	// if ((file = fopen(strcat("log_rank", rank_str), "ab")) == NULL) {
	// 	printf("Error of open the file.\n");
	// }

	int size_of_message = 1;
	int old_size;
	unsigned char *message = calloc(size_of_message, sizeof(char));
	if (rank < 8) {
		*message |= 0x1 << rank;
	} else {
		*message = 0;
	}
	// printf("[%d] START position = %s\n", rank, printArrayBytes(message, size_of_message));
	
	// unsigned char *buf = calloc(size_of_message, sizeof(char));

	double time = 0;
		
	// printf("I %d of %d\n", rank, max_rank);

	int count = 0;
	for (; count < 3; ) {
		if (rank == 0) {
			time -= MPI_Wtime();
		}
		for (int i = 0; i < max_rank - 1; i++) {
			if (rank == max_rank - 1) {
				// MPI_Send(message, size_of_message, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
				// MPI_Recv(message, size_of_message, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(message, size_of_message, MPI_CHAR, 0, 0, message, size_of_message, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("[%d] recv: %d\n", rank, *message);
				// printf("[%d] recv: %s\n", rank, printArrayBytes(message, size_of_message));
				// printArrayBytes(message, size_of_message);
			} else if (rank == 0) {
				// MPI_Send(message, size_of_message, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD);
				// MPI_Recv(message, size_of_message, MPI_CHAR, max_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(message, size_of_message, MPI_CHAR, rank + 1, 0, message, size_of_message, MPI_CHAR, max_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("[%d] recv: %d\n", rank, *message);
				// printf("[%d] recv: %s\n", rank, printArrayBytes(message, size_of_message));
				// printArrayBytes(message, size_of_message);
			} else {
				// MPI_Send(message, size_of_message, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD);
				// MPI_Recv(message, size_of_message, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(message, size_of_message, MPI_CHAR, rank + 1, 0, message, size_of_message, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// printf("[%d] recv: %d\n", rank, *message);
				// printf("[%d] recv: %s\n", rank, printArrayBytes(message, size_of_message));
				// printArrayBytes(message, size_of_message);
			}
		}
		if (rank == 0) {
			time += MPI_Wtime();
			printf("Time elapsed [count = %d]: %f\n", size_of_message, time);
			printf("Stage %d complete\n", count);
		}

		count++;
		size_of_message = count == 1 ? 1024 : 1024 * 1024;
		message = (char *)realloc(message, size_of_message);

		int st = (size_of_message / max_rank) * rank;
		int fn = st + (size_of_message / max_rank);
		for (int i = st; i < fn; i++) {
			message[i] = 255;
		}
		// printArrayBytes(message, size_of_message);

		/*
		count++;
		old_size = size_of_message;
		size_of_message = 1024 ? count == 1 : 1024 * 1024;
		realloc(message, size_of_message);
		realloc(buf, size_of_message);
		memset(message, 0, old_size);
		memset(buf, 0, old_size);
		for (int j = rank * 64; j < (rank + 1) * 64; j++) {
			message[j] <= 8;
		}
		*/
	}
	
	MPI_Finalize();
	return 0;
}