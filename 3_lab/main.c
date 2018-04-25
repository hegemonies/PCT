#include <stdio.h>
#include <math.h>
#include <omp.h>

double function(double x)
{
	return (1 - exp(0.7 / x)) / (2 + x);
}

int main()
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

	return 0;
}