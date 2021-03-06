#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <limits.h>
#include <time.h>

#define INF_PERC 30

static int rank;
static int commsize;

static double mpi_total_time;

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
                arr[i * n + j] = INT_MAX;
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

void row_distr(int *arr, int n, int k, int *row) {
    int num[commsize];
    int ind[commsize];
    int rem = n % commsize;
    int count_rows = n / commsize;

    for (int i = 0; i < commsize; i++) {
        num[i] = count_rows;

        if (rem > 0) {
            num[i]++;
            rem--;
        }

        ind[i] = (i > 0) ? ind[i - 1] + num[i - 1] : 0;

        #if 0
        if (rank == 0) {
            printf("k = %d [%d] num[%d] = %d\t", k, rank, i, num[i]);
            printf("ind[%d] = %d\n", i, ind[i]);
        }
        #endif
    }

    int row_rank = -1;

    for (int i = 0; i < commsize; i++) {
        if (k < ind[i] + num[i]) {
            row_rank = i;
            break;
        }
    }

    #if 0
    printf("[%d] k = %d row_rank = %d\n", rank, k, row_rank);
    #endif

    if (row_rank == rank) {
        copy(&arr[(k - ind[rank]) * n], n, row);
    }

    mpi_total_time -= MPI_Wtime();
    MPI_Bcast(row, n, MPI_INT, row_rank, MPI_COMM_WORLD);
    mpi_total_time += MPI_Wtime();
}

void par_Floyd(int *arr, int n, int count_rows) {
    int *row = calloc(n, sizeof(int));
    if (!row) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int k = 0; k < n; k++) {
        row_distr(arr, n, k, row);

        for (int i = 0; i < count_rows; i++) {
            for (int j = 0; j < n; j++) {
                if ((arr[i * n + k] != INT_MAX) && (row[j] != INT_MAX)) {
                    arr[i * n + j] = min(arr[i * n + j], arr[i * n + k] + row[j]);
                }
            }
        }
    }

    free(row);
}

int compare(int *a, int *b, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (a[i * n + j] != b[i * n + j]) {
                return -1;
            }
        }
    }

    return 0;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    mpi_total_time = 0;
    double t = 0;
    t -= MPI_Wtime();
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int n;
    int rem;
    int *arr;
    // int *cp_arr;
    int *recv_arr;
    int count_rows;
    int real_count_rows;

    if (rank == 0) {
        n = (argc > 1) ? atoi(argv[1]) : 0;
        if (n == 0) {
            printf("How to run:\nmpiexec ./main <number of vertices>\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        } else if (n < commsize) {
            printf("Need number of vertices bigger then number of processors\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        arr = calloc(n * n, sizeof(int));
        if (!arr) {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // data_random_init(arr, n);
        dummy_data_init(arr, n);

        #if 0
        printf("arr:\n");
        print_matrix(arr, n);
        #endif

        #if 0
        cp_arr = calloc(n * n, sizeof(int));
        if (!cp_arr) {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        copy(arr, n * n, cp_arr);

        // printf("\nbefore floyd cp_arr:\n");
        // print_matrix(cp_arr, n);

        serial_Floyd(cp_arr, n);

        // printf("\nafter floyd cp_arr:\n");
        // print_matrix(cp_arr, n);
        #endif
    }

    mpi_total_time -= MPI_Wtime();
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    mpi_total_time += MPI_Wtime();

    real_count_rows = n / commsize;
    rem = n % commsize;
    
    if (rem > rank) {
        real_count_rows++;
    }

    count_rows = real_count_rows;

    #if 0
    printf("[%d] count_rows = %d\n", rank, count_rows);
    #endif

    recv_arr = calloc(n * count_rows, sizeof(int));
    if (!recv_arr) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int *send_num = calloc(commsize, sizeof(int)); 
    int *send_ind = calloc(commsize, sizeof(int)); 
    if (!send_num || !send_ind) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int old_rem = rem;
    #if 0
    if (rank == 0) {
        printf("rem = %d\n", rem);
    }
    #endif

    count_rows = n / commsize;

    for (int i = 0; i < commsize; i++) {
        send_num[i] = count_rows * n;

        if (rem > 0) {
            send_num[i] += n;
            rem--;
        }

        send_ind[i] = (i > 0) ? send_ind[i - 1] + send_num[i - 1] : 0;
        #if 0
        if (rank == 0) {
            printf("rem = %d\n", rem);
        }
        #endif
    }

    rem = old_rem;
    count_rows = real_count_rows;
    
    mpi_total_time -= MPI_Wtime();
    MPI_Scatterv(arr, send_num, send_ind, MPI_INT, recv_arr, send_num[rank], MPI_INT, 0, MPI_COMM_WORLD); 
    mpi_total_time += MPI_Wtime();

    #if 0
    if (rank == 0) {
        printf("\nprint before par_Floyd\n");
    }
    par_print_matrix(recv_arr, n, count_rows, rank, commsize);
    printf("\n");
    #endif

    count_rows = real_count_rows;
    
    par_Floyd(recv_arr, n, count_rows);
    
    #if 0
    if (rank == 0) {
        printf("\nprint after par_Floyd\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    par_print_matrix(recv_arr, n, count_rows, rank, commsize);
    #endif

    int *recv_num = calloc(commsize, sizeof(int));
    int *recv_ind = calloc(commsize, sizeof(int));
    if (!recv_num || !recv_ind) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    #if 0
    if (rank == 0) {
        printf("rem = %d\n", rem);
    }
    #endif

    count_rows = n / commsize;

    for (int i = 0; i < commsize; i++) {
        recv_num[i] = count_rows * n;

        if (rem > 0) {
            recv_num[i] += n;
            rem--;
        }

        recv_ind[i] = (i > 0) ? recv_ind[i - 1] + recv_num[i - 1] : 0;
        #if 0
        if (rank == 0) {
            printf("rem = %d\n", rem);
        }
        #endif
    }

    count_rows = real_count_rows;

    #if 0
    if (rank == 0) {
        for (int i = 0; i < commsize; i++) {
            printf("recv_ind[%d] = %d\n", i, recv_ind[i]);
            printf("recv_num[%d] = %d\n", i, recv_num[i]);
        }
    }
    #endif

    mpi_total_time -= MPI_Wtime();
    MPI_Gatherv(recv_arr, real_count_rows, MPI_INT, arr, recv_num, recv_ind, MPI_INT, 0, MPI_COMM_WORLD);
    mpi_total_time += MPI_Wtime();

    #if 0
    if (rank == 0) {
        printf("\nprint ser after gather\n");
        print_matrix(arr, n);
    }
    #endif


    #if 0
    if (rank == 0) {
        if (compare(arr, cp_arr, n) == 0) {
            printf("Compare is bad\n");
        } else {
            printf("Compare is good\n");
            // printf("Elapsed time is %.5f sec\n", t);
        }
    }
    #endif

    t += MPI_Wtime();

    double total_time_max = 0;
    MPI_Reduce(&t, &total_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    double total_wtime_max = 0;
    MPI_Reduce(&mpi_total_time, &total_wtime_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    #if 1
    if (rank == 0) {
        printf("Total time is %.5f sec\n", total_time_max);
        printf("Total mpi time is %.5f sec\n", total_wtime_max);
        printf("Total compute time is %.5f sec\n", total_time_max - total_wtime_max);
        printf("coefficient is %.4f \n", (double)total_wtime_max / (double)(total_time_max - total_wtime_max));
    }
    #endif

	MPI_Finalize();

    return 0;
}