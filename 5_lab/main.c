#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

// #define N 10000000
#define THRESHOLD 1000

void swap(int *a, int *b)
{
	int tmp = *a;
	*a = *b;
	*b = tmp;
}

void partition(int *array, int *i, int *j, int low, int high)
{
	*i = low;
	*j = high;
	int pivot = array[(low + high) / 2];
	do {
		while (array[*i] < pivot) 
			(*i)++;
		while (array[*j] > pivot)
			(*j)--;
		if ((*i) <= (*j)) {
			swap(&array[*i], &array[*j]);
			(*i)++;
			(*j)--;
		}
	} while ((*i) <= (*j));
}

void quicksort(int *array, int low, int high)
{
	int i = 0;
	int j = 0;
	
	partition(array, &i, &j, low, high);

	if (high - low < THRESHOLD || j - low < THRESHOLD || high - i < THRESHOLD) {
		if (low < j) {
			quicksort(array, low, j);
		}
		if (i < high) {
			quicksort(array, i, high);
		}
	} else {
		#pragma omp task untied
		{
			quicksort(array, low, j);
		}
		quicksort(array, i, high);
	}

}

void init_array(int *array, int N)
{
	for (int i = 0; i < N; i++) {
		array[i] = rand() / (RAND_MAX + 1.0) * (100000 - 1) + 1;;
	}
}

int main()
{
	int N[6] = { 1000, 10000, 100000, 1000000, 10000000, 100000000 };
	int nthreads[7] = { 1, 2, 4, 6, 8, 16, 32 };

	for (int j = 0; j < 6; j++) {
		printf("N = %d\n", N[j]);

		int *orig_arr = malloc(sizeof(int) * N[j]);
		int *for_test_arr = malloc(sizeof(int) * N[j]);

		init_array(orig_arr, N[j]);

		for (int i = 0; i < 7; i++) {
			memcpy(for_test_arr, orig_arr, sizeof(int) * N[j]);
			
			double t = omp_get_wtime();
			#pragma omp parallel num_threads(nthreads[i])
			{
				#pragma omp single
				quicksort(for_test_arr, 0, N[j] - 1);
			}
			t = omp_get_wtime() - t;
			printf("Time for %d threads = %.5lf\n", nthreads[i], t);
		}

		printf("\n");
		free(orig_arr);
		free(for_test_arr);
	}

	// printf("orig_arr = ");
	// for (int i = 0; i < N; i++)
	// 	printf("%d ", orig_arr[i]);
	// printf("\n");

	// memcpy(for_test_arr, orig_arr, sizeof(int) * N);

	// printf("for_test_arr = ");
	// for (int i = 0; i < N; i++)
	// 	printf("%d ", for_test_arr[i]);
	// printf("\n");

	// double t = omp_get_wtime();
	// #pragma omp parallel num_threads(nthreads[0])
	// {
	// 	#pragma omp single
	// 	quicksort(for_test_arr, 0, N - 1);
	// }
	// t = omp_get_wtime() - t;

	// printf("After the quicksort = ");
	// for (int i = 0; i < N; i++)
	// 	printf("%d ", for_test_arr[i]);
	// // printf("\n");

	// printf("Time for %d threads = %.5lf\n", nthreads[0], t);

	return 0;
}