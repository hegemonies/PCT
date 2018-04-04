#include <stdio.h>
#include "hpctimer.h"
#include <omp.h>

enum { 
    N = 4,
    NREPS = 3
};

double A[N * N], B[N * N], C[N * N];

void init_matrix(double *a, double *b, double *c, int n)
{
	int i, j, k;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
                *(a + i * n + j) = 1.0;
                *(b + i * n + j) = 2.0;
                *(c + i * n + j) = 0.0;
			}
		}
	}
}

void dgemm_transpose(double *a, double *b, double *c, int st, int fn)
{
    int i, j, k;
    
    for (i = st; i < fn; i++) {
        for (k = 0; k < fn; k++) {
            for (j = 0; j < fn; j++) {
                *(c + i * fn + j) += *(a + i * fn + k) * *(b + k * fn + j);
            }
        }
    }
}

void printMatrix(double *a, int n)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%.2lf ", *(a + i * n + j));
		}
		printf("\n");
	}
	printf("\n");
}

int main()
{
	init_matrix(A, B, C, N);

	double C1[N * N];
	init_matrix(A, B, C1, N);
	dgemm_transpose(A, B, C1, 0, N);


	double t = hpctimer_getwtime();
	for (int i = 0; i < NREPS; i++) {
		#pragma omp parallel
		{
			int num_threads = omp_get_num_threads();
			int thread_num = omp_get_thread_num();
			int num_rows_for_one_thread = N / num_threads;

			if (thread_num == 0) {
				dgemm_transpose(A, B, C, 0, num_rows_for_one_thread);
			}
			if (thread_num == 1) {
				dgemm_transpose(A, B, C, num_rows_for_one_thread + 1, num_rows_for_one_thread * 2);
			}
			if (thread_num == 2) {
				dgemm_transpose(A, B, C, (num_rows_for_one_thread * 2) + 1, num_rows_for_one_thread * 3);
			}
			if (thread_num == 3) {
				dgemm_transpose(A, B, C, (num_rows_for_one_thread * 3) + 1, N);
			}
			// for (int j = 0; j < num_threads; j++) {
			// 	if (thread_num == j) {
			// 		dgemm_transpose(A, B, C, (num_rows_for_one_thread * j) + ((j == 0) ? 0 : 1), num_rows_for_one_thread * (j + 1));
			// 	}
			// }
		}
		
		// for (int j = 0;)
	}
	t = hpctimer_getwtime() - t;
	t = t / NREPS;

	printMatrix(C, N);
	printMatrix(C1, N);

	printf("Elapsed time: %.6f sec.\n", t);

	int check = 0;
	for (int i = 0; i < N * N; i++) {
		if (C1[i] != C[i]) {
			check++;
		}
	}
	printf("Errors %d\n", check);

	return 0;
}