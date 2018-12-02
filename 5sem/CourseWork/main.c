#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <limits.h>
#include <time.h>

#define INF_PERC 40

static int rank;

int min(int a, int b) {
    return (a < b) ? a : b;
}

void data_random_init(int *arr, int n) {
    srand(time(0));

    int rnd;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                rnd = rand() % 1000;
                if ((rnd % 100) < INF_PERC) {
                    arr[i * n + j] = INT_MAX;
                } else {
                    arr[i * n + j] = rnd + 1;
                }
            } else {
                arr[i * n + j] = 0;
            }
        }
    }
}

void dummy_data_init(int *arr, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (i == j) {
                arr[i * n + j] = 0;
            } else if (i == 0) {
                arr[i * n + j] = j;
            } else {
                arr[i * n + j] = -1;
            }
            arr[j * n + i] = arr[i * n + j];
        }
    }
} 

void serial_Floyd(int *arr, int n) {
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((arr[i * n + k] != INT_MAX) && (arr[k * n + j] != INT_MAX)) {
                    arr[i * n + j] = min(arr[i * n + j], arr[i * n + k] + arr[k * n + j]);
                }
            }
        }
    }
}

void print_matrix(int *arr, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (arr[i * n + j] == INT_MAX) {
                printf("INF\t");
            } else {
                printf("%d\t", arr[i * n + j]);
            }
        }
        printf("\n");
    }
}

void par_print_matrix(int *arr, int n, int count_rows, int rank, int commsize) {
    for (int p = 0; p < commsize; p++) {
        if (p == rank) {
            printf("rank = %d\n", rank);
            for (int i = 0; i < count_rows; i++) {
                for (int j = 0; j < n; j++) {
                    if (arr[i * n + j] == INT_MAX) {
                        printf("INF\t");
                    } else {
                        printf("%d\t", arr[i * n + j]);
                    }
                }
                printf("\n");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void copy(int *start, int n, int *out) {
    if (!start || !out) {
        return;
    }

    for (int i = 0; i < n; i++) {
        out[i] = start[i];
    }
}

void row_distr(int *arr, int n, int count_rows, int k, int *row) {
    int row_rank = k / count_rows;
    int row_num = k - row_rank * count_rows;

    if (row_rank == rank) {
        copy(&arr[row_num * n], n, row);
    }

    MPI_Bcast(row, n, MPI_INT, row_rank, MPI_COMM_WORLD); 
}

void par_Floyd(int *arr, int n, int count_rows) {
    int *row = calloc(sizeof(int), n);

    for (int k = 0; k < n; k++) {
        row_distr(arr, n, count_rows, k, row);

        for (int i = 0; i < count_rows; i++) {
            for (int j = 0; j < n; j++) {
                if ((arr[i * n + j] != INT_MAX) && (row[j] != INT_MAX)) {
                    arr[i * n + j] = min(arr[i * n + j], arr[i * n + k] + row[j]);
                }
            }
        }
    }

    free(row);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
	
	// int rank;
	int commsize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);


    int n;
    int rem;
    int *arr;
    int *recv_arr;
    int count_rows;

    if (rank == 0) {
        n = (argc > 1) ? atoi(argv[1]) : 0;
        if (n == 0) {
            printf("How to run:\nmpiexec ./main <number of vertices>\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        } else if (n < commsize) {
            printf("Need number of vertices bigger then number of processors\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        rem = n % commsize;

        arr = calloc(sizeof(int), n * n);
        // data_random_init(arr, n);
        dummy_data_init(arr, n);

        print_matrix(arr, n);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    count_rows = n / commsize;

    if (rem > rank) {
        count_rows++;
    }

    recv_arr = calloc(sizeof(int), n * count_rows);
    // printf("rank = %d  n * count_rows = %d\n", rank, n * count_rows);

    int *sendsize = calloc(sizeof(int), commsize);
    int *displs = calloc(sizeof(int), commsize);

    for (int i = 0; i < commsize; i++) {
        sendsize[i] = (n * n) / commsize;

        if (rem > 0) {
            sendsize[i] += n;
            rem--;
            sendsize[i]--;
        }

        displs[i] = (i > 0) ? displs[i - 1] + sendsize[i - 1] : 0;
        if (rank == 0) {
            // printf("sendsize[%d] = %d\n", i, sendsize[i]);
            // printf("displs[%d] = %d\n", i, displs[i]);
        }
    }

    MPI_Scatterv(arr, sendsize, displs, MPI_INT, recv_arr, n * count_rows, MPI_INT, 0, MPI_COMM_WORLD);

    // MPI_Barrier(MPI_COMM_WORLD);
    // par_print_matrix(recv_arr, n, count_rows, rank, commsize);
    printf("\n");
    // if (rank == 0) {
    //     print_matrix(arr, n);
    // }

    par_Floyd(arr, n, count_rows);

    // MPI_Barrier(MPI_COMM_WORLD);
    // par_print_matrix(recv_arr, n, count_rows, rank, commsize);

    int *recv_num = calloc(sizeof(int), commsize);
    int *recv_index = calloc(sizeof(int), commsize);
    int rest_rows = n; 
    count_rows = n / commsize;
    recv_index[0] = 0;
    recv_num[0] = count_rows * n;
    for(int i = 1; i < commsize; i++) {
        rest_rows -= count_rows;
        count_rows = rest_rows / (commsize - i);
        recv_num[i] = count_rows * n;
        recv_index[i] = recv_index[i - 1] + recv_num[i - 1];
    }
    MPI_Gatherv(recv_arr, recv_num[rank], MPI_INT, 
    arr, recv_num, recv_index, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        print_matrix(arr, n);
    }

	MPI_Finalize();

    return 0;
}