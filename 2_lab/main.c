#include <stdio.h>
#include "hpctimer.h"
#include <omp.h>

enum { 
	N = 15000,
    NREPS = 1
};

double A[N * N], B[N], C[N], C1[N];

void init_matrix(double *a, double *b, double *c, int n)
{
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
            *(a + i * n + j) = 1.0;
            *(b + i * n + j) = 2.0;
            *(c + i * n + j) = 0.0;
		}
	}
}

void dgemv(double *m, double *v, double *c, int st, int fn, int n)
{
	for (int i = st; i < fn; i++) {
		for (int j = 0; j < n; j++) {
			#pragma omp atomic
			c[j] += *(m + i * n + j) * *(v + j);
		}
	}
}

void dgemm_transpose(double *a, double *b, double *c, int st, int fn, int n)
{
    int i, j, k;
    
    for (i = st; i < fn; i++) {
        for (k = 0; k < n; k++) {
            for (j = 0; j < n; j++) {
            	#pragma omp atomic
                *(c + i * n + j) += *(a + i * n + k) * *(b + k * n + j);
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

void printVector(double *v, int n)
{
	for (int i = 0; i < n; i++) {
		printf("%.3lf ", v[i]);
	}
	printf("\n");
}

// int main()
// {
// 	#pragma omp parallel
// 	{
// 		printf("Hey, thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
// 	}

// 	return 0;
// }


int main()
{
	// init_matrix(A, B, C, N);

	// init_matrix(A, B, C1, N);
	for (int i = 0; i < N; i++) {
		B[i] = 2.0;
		C[i] = 0.0;
		for (int j = 0; j < N; j++) {
			*(A + i * N + j) = 1.0;
		}
	}
/*
	double t = hpctimer_getwtime();
	for (int i = 0; i < NREPS; i++) {
		// dgemm_transpose(A, B, C1, 0, N, N);
		dgemv(A, B, C1, 0, N, N);
	}
	t = hpctimer_getwtime() - t;
	t = t / NREPS;
	printf("Elapsed time: %.6f sec.\n", t);
*/

	// dgemm_transpose(A, B, C1, 0, N);

	printf("\nN = %d\n", N);
	double t1 = omp_get_wtime();
	#pragma omp parallel
	{
		#pragma omp master
		printf("Num Threads = %d\n", omp_get_num_threads());
		int num_threads = omp_get_num_threads();
		int thread_num = omp_get_thread_num();

		int num_rows_for_one_thread = N / num_threads;
		int st = (num_rows_for_one_thread * thread_num);
		int fn = num_rows_for_one_thread * (thread_num + 1);

		// printf("%d :: %d = %d %d %d\n", num_threads, thread_num, num_rows_for_one_thread, st, fn);

		// for (int i = 0; i < NREPS; i++) {
		// 	if (N % num_threads == 0 ) {
		// 		dgemm_transpose(A, B, C, st, fn, N);
		// 	} else if (thread_num == num_threads - 1) {
		// 		dgemm_transpose(A, B, C, st, fn, N);
		// 	}
		// }
		for (int i = 0; i < NREPS; i++) {
			if (N % num_threads == 0) {
				dgemv(A, B, C, st, fn, N);
			} else if (thread_num == num_threads - 1) {
				dgemv(A, B, C, st, fn, N);
			}
		}

		#pragma omp barrier
	}
	t1 = omp_get_wtime() - t1;
	t1 = t1 / NREPS;


	printf("Elapsed time: %.6f sec.\n", t1);
	// printf("Speedup = %.4lf\n", t / t1);


	// printMatrix(C, N);
	// printMatrix(C1, N);

	// printf("Elapsed time: %.6f sec.\n", t);

	// int check = 0;
	// for (int i = 0; i < N * N; i++) {
	// 	if (C1[i] != C[i]) {
	// 		check++;
	// 	}
	// }
	// printf("Errors %d\n", check);

	return 0;
}
