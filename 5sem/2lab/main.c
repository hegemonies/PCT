#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

double function(double x)
{
	return pow(x, 4) / (0.5 * pow(x, 2) + x + 6);
}

void Runge(int rank, int max_rank)
{
	double a = 0.4;
	double b = 1.5;
	int n0 = 100000000;
	double eps = 1E-6;

	if (rank == 0) {
		printf("Numerical integration: [%f, %f], n0 = %d, EPS = %f\n", a, b, n0, eps);
	}

	double sq[2];

	int n = n0;
	int k;
	double delta = 1;

	for (k = 0; delta > eps; n *= 2, k ^= 1) {
		double sloc = 0.0;
		double h = (b - a) / n;
		
		int points_per_proc = n / max_rank;
		int lb = rank * points_per_proc;
		int ub = (rank == max_rank - 1) ? (n - 1) : (lb + points_per_proc - 1);

		sq[k] = 0.0;
		for (int i = lb; i <= ub; i++) {
			sloc += function(a + h * (i + 0.5));
		}
		
		double res_sum = 0;
		MPI_Allreduce(&sloc, &sq[k], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		sq[k] *= h;
		
		if (n > n0) {
			delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
		}
	}

	if (rank == 0) {
		printf("Result: %.12f; Runge rule: EPS %e, n %d\n", sq[k], eps, n / 2);
	}
}

double getrand(unsigned int *seed)
{
	return (double)rand_r(seed) / RAND_MAX;
}

double func(double x, double y)
{
	return exp(x - y);
}

void Monte_Carlo(int rank, int max_rank, int n)
{
	srand(rank);
	int in = 0;
	double s = 0;

	double s_loc = 0.0;
	int in_loc = 0;
	unsigned int seed = rank;

	int points_per_proc = n / max_rank;
	int lb = rank * points_per_proc;
	int ub = (rank == max_rank - 1) ? (n - 1) : (lb + points_per_proc - 1);

	for (int i = lb; i < ub; i++) {
		double x = getrand(&seed) * (-1);
		double y = getrand(&seed);
		in_loc++;
		s_loc += func(x, y);
	}

	MPI_Reduce(&in_loc, &in, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&s_loc, &s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double v =  (1.0 * in) / n;
		double res = (v * s) / in;
		printf("Result: %.12f, n %d \n", res, n);
	}
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	int max_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &max_rank);
	
	double t = 0;

	if (rank == 0) {
		t -= MPI_Wtime();
	}

	Runge(rank, max_rank);
	
	if (rank == 0) {
		t += MPI_Wtime();
		printf("Elapsed time (sec.): %.6f\n", t);
	}

	for (int i = 0; i < 2; i++) {
		int n = pow(10, 7 + i);
		t = 0;

		if (rank == 0) {
			t -= MPI_Wtime();
		}

		Monte_Carlo(rank, max_rank, n);
		
		if (rank == 0) {
			t += MPI_Wtime();
			printf("Elapsed time (sec.): %.6f\n", t);
		}
	}

	MPI_Finalize();
	return 0;
}
