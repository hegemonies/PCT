/* 
 * dgemm.c: DGEMM - Double-precision General Matrix Multiply.
 *
 */
 
#include <stdio.h>
#include <stdlib.h>

#include "hpctimer.h"

enum { 
    N = 1024,
    NREPS = 3
};

int BS = 2;

double A[N * N], B[N * N], C[N * N];

void dgemm_def(double *a, double *b, double *c, int n)
{
    int i, j, k;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                *(c + i * n + j) += *(a + i * n + k) * *(b + k * n + j);
			}
		}
	}
}

void dgemm_transpose(double *a, double *b, double *c, int n)
{
    int i, j, k;
    
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            for (j = 0; j < n; j++) {
                *(c + i * n + j) += *(a + i * n + k) * *(b + k * n + j);
            }
        }
    }
}

void dgemm_block(double *a, double *b, double *c, int n)
{
    int i, j, k;
    int i0, j0, k0;
    double *a0, *b0, *c0;

    for (i = 0; i < n; i += BS) {
        for (j = 0; j < n; j += BS) {
            for (k = 0; k < n; k += BS) {
                for (i0 = 0, c0 = (c + i * n + j),
                a0 = (a + i * n + k); i0 < BS;
                ++i0, c0 += n, a0 += n)
                {
                    for (k0 = 0, b0 = (b + k * n + j);
                    k0 < BS; ++k0, b0 += n)
                    {
                        for (j0 = 0; j0 < BS; ++j0) {
                            c0[j0] += a0[k0] * b0[j0];
                        }
                    }
                }
            }
        }
    }
}
void Rec_Mult(double *C, const double *A, const double *B, int n, int rowsize)
{
    if (n == 2)
    {
        const int d11 = 0;
        const int d12 = 1;
        const int d21 = rowsize;
        const int d22 = rowsize + 1;
 
        C[d11] += A[d11] * B[d11] + A[d12] * B[d21];
        C[d12] += A[d11] * B[d12] + A[d12] * B[d22];
        C[d21] += A[d21] * B[d11] + A[d22] * B[d21];
        C[d22] += A[d21] * B[d12] + A[d22] * B[d22];
    }
    else
    {
        const int d11 = 0;
        const int d12 = n / 2;
        const int d21 = (n / 2) * rowsize;
        const int d22 = (n / 2) * (rowsize + 1);
 
        // C11 += A11 * B11
        Rec_Mult(C + d11, A + d11, B + d11, n / 2, rowsize);
        // C11 += A12 * B21
        Rec_Mult(C + d11, A + d12, B + d21, n / 2, rowsize);
 
        // C12 += A11 * B12
        Rec_Mult(C + d12, A + d11, B + d12, n / 2, rowsize);
        // C12 += A12 * B22
        Rec_Mult(C + d12, A + d12, B + d22, n / 2, rowsize);
 
        // C21 += A21 * B11
        Rec_Mult(C + d21, A + d21, B + d11, n / 2, rowsize);
        // C21 += A22 * B21
        Rec_Mult(C + d21, A + d22, B + d21, n / 2, rowsize);
 
        // C22 += A21 * B12
        Rec_Mult(C + d22, A + d21, B + d12, n / 2, rowsize);
        // C22 += A22 * B22
        Rec_Mult(C + d22, A + d22, B + d22, n / 2, rowsize);
    }
}


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

void print_matrix(double *a, int n)
{
	int i, j;
	
	printf("Matrix:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%12.2f", *(a + i * n + j));
		}
		printf("\n");
	}
}

int main(int argc, char **argv)
{
    int i;
    double t;
    		
    init_matrix(A, B, C, N);

    // int size_min_matx = 2;

    // print_matrix(A, N);
    // printf("\n");
    // print_matrix(B, N);

    // t = hpctimer_getwtime();
    // dgemm_def(A, B, C, N);
    // t = hpctimer_getwtime() - t;
    // printf("Elapsed time: %.6f sec.\n", t);
    // print_matrix(C, N);

   // for (int j = 0; j < 7; j++) {
        t = hpctimer_getwtime();
        for (i = 0; i < NREPS; i++) {
            //dgemm_def(A, B, C, N);
            dgemm_transpose(A, B, C, N);
            //dgemm_transpose2(A, B, C, N);
            //dgemm_block(A, B, C, N);
        }
        t = hpctimer_getwtime() - t;
        t = t / NREPS;
        // printf("BS = %d\n", BS);
        printf("Elapsed time: %.6f sec.\n", t);

        //BS *= 2;
    //}
/*
    init_matrix(A, B, C, N);

    t = hpctimer_getwtime();
    for (i = 0; i < NREPS; i++) {
        //dgemm_def(A, B, C, N);
        dgemm_transpose(A, B, C, N);
        // dgemm_transpose2(A, B, C, N);
        //dgemm_block(A, B, C, N);
    }
    t = hpctimer_getwtime() - t;
    t = t / NREPS;
    // printf("BS = %d\n", BS);
    printf("Elapsed time: %.6f sec.\n", t);
*/
    // print_matrix(C, N);    

    // init_matrix(A, B, C, N);

    // printf("\nMy: \n");

    // t = hpctimer_getwtime();
    // for (i = 0; i < NREPS; i++) {
    //     // dgemm_block(A, B, C, N);
    //     Rec_Mult(C, A, B, N, N);
    // }
    // t = hpctimer_getwtime() - t;
    // t = t / NREPS;
    // printf("Elapsed time: %.6f sec.\n", t);

    // print_matrix(C, N);
    
    // printf("Elapsed time: %.6f sec.\n", t);

    return 0;
}

