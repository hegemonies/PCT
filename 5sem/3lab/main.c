#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

void get_chunk(int a, int b, int commsize, int rank, int *lb, int *ub)
{
	int n = b - a + 1;
	int q = n / commsize;
	if (n % commsize) {
		q++;
	}
	int r = commsize * q - n;

	int chunk = q;
	if (rank >= commsize - r) {
		chunk = q - 1;
	}

	*lb = a;
	if (rank > 0) {
		if (rank <= commsize - r) {
			*lb += q * rank;
		} else {
			*lb += q * (commsize - r) + (q - 1) * (rank - (commsize - r));
		}
	}
	*ub = *lb + chunk - 1;
}

void dgemv(double *a, double *b, double *c, int m, int n, int lb, int ub, int rank, int max_rank)
{
	int nrows = ub - lb + 1;

	for (int i = 0; i < nrows; i++) {
		c[lb + i] = 0.0;
		for (int j = 0; j < n; j++) {
			c[lb + i] += a[i * n + j] * b[j];
		}
	}

	if (rank == 0) {
		int *displs = malloc(sizeof(*displs) * max_rank);
		int *rcounts = malloc(sizeof(*rcounts) * max_rank);
		for (int i = 0; i < max_rank; i++) {
			int l, u;
			get_chunk(0, m - 1, max_rank, i, &l, &u);
			rcounts[i] = u - l + 1;
			displs[i] = (i > 0) ? displs[i - 1] + rcounts[i - 1] : 0;
		}
		MPI_Gatherv(MPI_IN_PLACE, ub - lb + 1, MPI_DOUBLE, c, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} else {
		MPI_Gatherv(&c[lb], ub - lb + 1, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	int max_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &max_rank);

	int n = 45000;
	int m = n;
	
	double t = 0;
	t -= MPI_Wtime();

	int lb, ub;
	get_chunk(0, m - 1, max_rank, rank, &lb, &ub);

	int nrows = ub - lb + 1;

	double *a = malloc(sizeof(*a) * nrows * n);
	double *b = malloc(sizeof(*b) * n);
	double *c = malloc(sizeof(*c) * m);

	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < n; j++) {
			a[i * n + j] = lb + i + 1;
		}
	}
	for (int j = 0; j < n; j++) {
		b[j] = j + 1;
	}

	dgemv(a, b, c, m, n, lb, ub, rank, max_rank);

	t += MPI_Wtime();
	// printf("[%d] \t%d %d Elapsed time (sec.): %.6f\n", rank, m, n, t);
	
	if (rank == 0) {
		for (int i = 0; i < m; i++) {
			double r = (i + 1) * (n / 2.0 + pow(n, 2) / 2.0);
			if (fabs(c[i] - r) > 1E-6) {
				fprintf(stderr, "Validation failed: elem %d = %f (real value %f)\n", i, c[i], r);
				break;
			}
		}

		printf("DGEMV: matrix-vector product (c[m] = a[m, n] * b[n]; m = %d, n = %d)\n", m, n);
		printf("Memory used: %lu MiB\n", (uint64_t)(((double)m * n + m + n) * sizeof(double)) >> 20);
		double gflop = 2.0 * m * n * 1E-9;
		printf("Elapsed time (%d procs): %.6f sec.\n", max_rank, t);
		printf("Performance: %.2f GFLOPS\n", gflop / t);
	}

	free(a);
	free(b);
	free(c);

	MPI_Finalize();
	return 0;
}
