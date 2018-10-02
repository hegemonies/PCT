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
	// double t;
	// if (rank == 0) {
	// 	t = MPI_Wtime();
	// }
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
		// MPI_Reduce(&sloc, &res_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Allreduce(&sloc, &sq[k], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		sq[k] *= h;
		
		if (n > n0) {
			delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
		}
		/*
		if (rank == 0) {
			sq[k] += res_sum * h;
			
		}
		*/
	}

	if (rank == 0) {
		printf("Result: %.12f; Runge rule: EPS %e, n %d\n", sq[k], eps, n / 2);
		// t += MPI_Wtime();
		// printf("Elapsed time (sec.): %.6f\n", t);
	}
}

double func(double x, double y)
{
	return x / pow(y, 2);
}
/*
void Monte_Carlo(int x)
{
	double t = omp_get_wtime();
	int n = 100000000;

	int in = 0;
	double s = 0;

	{
		printf("num threads = %d\n", omp_get_num_threads());
		double s_loc = 0.0;
		int in_loc = 0;
		unsigned int seed = omp_get_thread_num();

		for (int i = 0; i < n; i++) {
			double x = getrand(&seed);
			double y = getrand(&seed) * 5.0;
			if (y >= 2) {
				in_loc++;
				s_loc += func(x, y);
			}
		}
		s += s_loc;
		in += in_loc;
	}

	double v =  (5.0 * in) / n;
	double res = (v * s) / in;

	t = MPI_Wtime() - t;

	printf("Result: %.12f, n %d \n", res, n);
	printf("Elapsed time = %lf\n", t);
}
*/

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	
	int rank;
	int max_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &max_rank);
	
	double t;
	if (rank == 0) {
		t -= MPI_Wtime();
	}

	Runge(rank, max_rank);
	
	if (rank == 0) {
		t += MPI_Wtime();
		printf("Elapsed time (sec.): %.6f\n", t);
	}

	MPI_Finalize();
	return 0;
}
