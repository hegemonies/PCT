#define _POSIX_C_SOURCE 1
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

double function(double x)
{
	return (1 - exp(0.7 / x)) / (2 + x);
}

void Runge()
{
	double t = omp_get_wtime();
	double a = 1.0;
	double b = 2.0;
	int n0 = 100000000;
	double eps = 1E-5;

	printf("Numerical integration: [%f, %f], n0 = %d, EPS = %f\n", a, b, n0, eps);

	double sq[2];

	#pragma omp parallel
	{
		int n = n0;
		int k;
		double delta = 1;

		for (k = 0; delta > eps; n *= 2, k ^= 1) {
			double sloc = 0.0;
			double h = (b - a) / n;

			sq[k] = 0.0;
			#pragma omp  barrier

			#pragma omp for nowait
			for (int i = 0; i < n; i++) {
				sloc += function(a + h * (i + 0.5));
			}

			#pragma omp atomic
			sq[k] += sloc * h;

			#pragma omp barrier
			if (n > n0) {
				delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
			}
		}

		#pragma omp master
		printf("Result: %.12f; Runge rule: EPS %e, n %d\n", sq[k], eps, n / 2);

	}

	t = omp_get_wtime() - t;
	printf("Elapsed time (sec.): %.6f\n", t);
}

double getrand(unsigned int *seed)
{
	return (double)rand_r(seed) / RAND_MAX;
}

double func(double x, double y)
{
	return x / pow(y, 2);
}

void Monte_Carlo()
{
	double PI = 3.14159265358979323846;
	int n = 10000000;

	int in = 0;
	double s = 0;

	#pragma omp parallel
	{
		double s_loc = 0.0;
		int in_loc = 0;
		unsigned int seed = omp_get_thread_num();

		#pragma omp for nowait
		for (int i = 0; i < n; i++) {
			double x = getrand(&seed);
			double y = getrand(&seed);
			if (y > 2 && y < 5) {
				in++;
				s_loc += func(x, y);
			}
		}
		#pragma omp atomic
		s += s_loc;
		#pragma omp atomic
		in += in_loc;
	}

	double v = PI * in / n;
	double res = v * s / in;

	printf("Result: %.12f, n %d \n", res, n);
}

int main()
{
	// Runge();
	Monte_Carlo();

	return 0;
}