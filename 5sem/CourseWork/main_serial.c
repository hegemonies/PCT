#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/time.h>
#include <time.h>
#include <inttypes.h>

#define INF_PERC 30

double w_time()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1E-6;
}

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

void copy(int *start, int n, int *out) {
    if (!start || !out) {
        return;
    }

    for (int i = 0; i < n; i++) {
        out[i] = start[i];
    }
}


int main(int argc, char **argv) {
    double t = 0;
    t -= w_time();
	
    int n;
    int *cp_arr;

    n = (argc > 1) ? atoi(argv[1]) : 0;
    if (n == 0) {
        printf("How to run:\nmpiexec ./main <number of vertices>\n");
        exit(1);
    }

    cp_arr = calloc(n * n, sizeof(int));
    if (!cp_arr) {
        exit(1);
    }

    // data_random_init(arr, n);
    dummy_data_init(cp_arr, n);

    // printf("\nbefore floyd cp_arr:\n");
    // print_matrix(cp_arr, n);

    serial_Floyd(cp_arr, n);

    // printf("\nafter floyd cp_arr:\n");
    // print_matrix(cp_arr, n);

    t += w_time();

    printf("Total time is %.5f sec\n", t);

    return 0;
}