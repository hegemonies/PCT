#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <limits.h>
#include <time.h>

#define INF_PERC 30

static int rank;
static int commsize;

int min(int a, int b) {
    return (a < b) ? a : b;
}

int Min(int A, int B) {
    int Result = (A < B) ? A : B;

    if((A < 0) && (B >= 0)) Result = B;
    if((B < 0) && (A >= 0)) Result = A;
    if((A < 0) && (B < 0)) Result = -1;

    return Result;
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

void SerialFloyd(int *arr, int n) {
    int t1, t2;
    for(int k = 0; k < n; k++)
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                if((arr[i * n + k] != INT_MAX) && (arr[k * n + j] != INT_MAX)) {
                    t1 = arr[i * n + j];
                    t2 = arr[i * n + k] + arr[k * n + j];
                    arr[i * n + j] = min(t1, t2);
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
    int row_rank;
    int rest_rows = n;
    int index = 0;
    int num = n / commsize;

    for (row_rank = 1; row_rank < commsize + 1; row_rank++) {
        if (k < index + num) {
            break;
        }
        rest_rows -= num;
        index += num;
        num = rest_rows / (commsize - row_rank);
    }
    row_rank = row_rank - 1;
    int row_num = k - index;

    if (row_rank == rank) {
        copy(&arr[row_num * n], n, row);
    }
    
    MPI_Bcast(row, n, MPI_INT, row_rank, MPI_COMM_WORLD);
}

void RowDistribution(int *pProcRows, int Size, int RowNum, int k, int *pRow) {
    int ProcRowRank; 
    int ProcRowNum; 

    int RestRows = Size;
    int Ind = 0;
    int Num = Size / commsize;

    for(ProcRowRank = 1; ProcRowRank < commsize + 1; ProcRowRank ++) {
        if(k < Ind + Num ) 
            break;
        RestRows -= Num;
        Ind += Num;
        Num = RestRows / (commsize - ProcRowRank);
    }
    ProcRowRank = ProcRowRank - 1;
    ProcRowNum = k - Ind;
    if(ProcRowRank == rank)
        copy(&pProcRows[ProcRowNum * Size], Size, pRow);
        
    MPI_Bcast(pRow, Size, MPI_INT, ProcRowRank, MPI_COMM_WORLD);
} 

void par_Floyd(int *arr, int n, int count_rows) {
    int *row = calloc(sizeof(int), n);

    for (int k = 0; k < n; k++) {
        // row_distr(arr, n, count_rows, k, row);
        RowDistribution(arr, n, count_rows, k, row);

        #if 0
        for (int p = 0; p < commsize; p++) {
            if (p == rank) {
                printf("rank = %d rr = ", rank);
                for (int i = 0; i < n; i++) {
                    if (row[i] == INT_MAX) {
                        printf("INF ");
                    } else {
                        printf("%d\t\t", row[i]);
                    }
                }
                printf("\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        #endif

        for (int i = 0; i < count_rows; i++) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if ((arr[i * n + k] != INT_MAX) && (row[j] != INT_MAX)) {
                        arr[i * n + j] = min(arr[i * n + j], arr[i * n + k] + row[j]);
                    }
                }
            }
        }
    }

    // free(row);
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
	
	// int rank;
	// int commsize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int n;
    int rem;
    int *arr;
    int *cp_arr;
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

        rem = n % commsize;

        arr = calloc(sizeof(int), n * n);
        // data_random_init(arr, n);
        dummy_data_init(arr, n);

        #if 0
        printf("arr:\n");
        print_matrix(arr, n);
        #endif

        cp_arr = calloc(sizeof(int), n * n);
        copy(arr, n * n, cp_arr);

        // printf("\nbefore floyd cp_arr:\n");
        // print_matrix(cp_arr, n);

        serial_Floyd(cp_arr, n);
        // SerialFloyd(cp_arr, n);

        // printf("\nafter floyd cp_arr:\n");
        // print_matrix(cp_arr, n);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    real_count_rows = n / commsize;
    
    if (rem > rank) {
        real_count_rows++;
    }

    count_rows = real_count_rows;

    #if 1
    printf("[%d] count_rows = %d\n", rank, count_rows);
    #endif

    recv_arr = calloc(sizeof(int), n * count_rows);
    // printf("rank = %d  n * count_rows = %d\n", rank, n * count_rows);

    // int *sendsize = calloc(sizeof(int), commsize);
    // int *displs = calloc(sizeof(int), commsize);

    // int f = 0;
    // if (rem > 0) {
    //     f = -1;
    // }

    // for (int i = 0; i < commsize; i++) {
    //     sendsize[i] = (n * n) / commsize;
    //     sendsize[i] += f;

    //     if (rem > 0) {
    //         sendsize[i] += n;
    //         rem--;
    //         // sendsize[i]--;
    //     }

    //     displs[i] = (i > 0) ? displs[i - 1] + sendsize[i - 1] : 0;
        // if (rank == 0) {
        //     printf("sendsize[%d] = %d\n", i, sendsize[i]);
        //     printf("displs[%d] = %d\n", i, displs[i]);
        // }
    // }

    // MPI_Scatterv(arr, sendsize, displs, MPI_INT, recv_arr, n * count_rows, MPI_INT, 0, MPI_COMM_WORLD);

    //

    int *send_num = calloc(sizeof(int), commsize); 
    int *send_ind = calloc(sizeof(int), commsize);
    int rest_rows = n;

    send_num[0] = count_rows * n; 
    send_ind[0] = 0;

    for (int i = 1; i < commsize; i++) {
        rest_rows -= count_rows;
        count_rows = rest_rows / (commsize - i);
        send_num[i] = count_rows * n;
        send_ind[i] = send_ind[i - 1] + send_num[i - 1];
    }

    #if 1
    if (rank == 0) {
        for (int i = 0; i < commsize; i++) {
            printf("send_num[%d] = %d\n", i, send_num[i]);
            printf("send_ind[%d] = %d\n", i, send_ind[i]);
        }
    }
    #endif
    
    MPI_Scatterv(arr, send_num, send_ind, MPI_INT, recv_arr, send_num[rank], MPI_INT, 0, MPI_COMM_WORLD); 
    //

    MPI_Barrier(MPI_COMM_WORLD);
    #if 0
    if (rank == 0) {
        printf("\nprint before par_Floyd\n");
    }
    par_print_matrix(recv_arr, n, count_rows, rank, commsize);
    printf("\n");
    #endif

    // if (rank == 0) {
    //     print_matrix(arr, n);
    // }

    // count_rows = n / commsize;
    count_rows = real_count_rows;
    MPI_Barrier(MPI_COMM_WORLD);
    par_Floyd(recv_arr, n, count_rows);

    printf("FLOYD OK\n");

    #if 0
    if (rank == 0) {
        printf("\nprint after par_Floyd\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    par_print_matrix(recv_arr, n, count_rows, rank, commsize);
    #endif

    int *recv_num = calloc(sizeof(int), commsize);
    int *recv_index = calloc(sizeof(int), commsize);
    
    rest_rows = n;
    // count_rows = n / commsize;

    recv_index[0] = 0;
    recv_num[0] = count_rows * n;

    for(int i = 1; i < commsize; i++) {
        rest_rows -= count_rows;
        count_rows = rest_rows / (commsize - i);
        recv_num[i] = count_rows * n;
        recv_index[i] = recv_index[i - 1] + recv_num[i - 1];
    }

    #if 1
    if (rank == 0) {
        for (int i = 0; i < commsize; i++) {
            printf("recv_num[%d] = %d\n", i, recv_num[i]);
            printf("recv_ind[%d] = %d\n", i, recv_index[i]);
        }
    }
    #endif

    // printf("[%d] step 0\n", rank);
    if (rank == 0) {
        MPI_Gatherv(MPI_IN_PLACE, recv_num[rank], MPI_INT, arr, recv_num, recv_index, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        // MPI_Gatherv(recv_arr, recv_num[rank], MPI_INT, arr, recv_num, recv_index, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(recv_arr, recv_num[rank], MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
    }

    printf("[%d] GATHER OK\n", rank);

    MPI_Barrier(MPI_COMM_WORLD);

    #if 1
    if (rank == 0) {
        printf("\nprint ser after gather\n");
        print_matrix(arr, n);
    }
    #endif

    if (rank == 0) {
        if (compare(arr, cp_arr, n)) {
            printf("Compare is bad\n");
        } else {
            printf("Compare is good\n");
        }
    }

	MPI_Finalize();

    return 0;
}