/*
 * loop.c:
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "hpctimer.h"

enum { n = 64 * 1024 * 1024 * 4 };

int main()
{
    int *v, i, sum;
    double t;
    
    if ( (v = malloc(sizeof(*v) * n)) == NULL) {
        fprintf(stderr, "No enough memory\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < n; i++)
        v[i] = 1;
   
    t = hpctimer_wtime();
    /* TODO: Unroll this loop */
    for (sum = 0, i = 0; i < n; i += 16) {
        sum += v[i];
        sum += v[i+1];
        sum += v[i+2];
        sum += v[i+3];
        sum += v[i+4];
        sum += v[i+5];
        sum += v[i+6];
        sum += v[i+7];
        sum += v[i+8];
        sum += v[i+9];
        sum += v[i+10];
        sum += v[i+11];
        sum += v[i+12];
        sum += v[i+13];
        sum += v[i+14];
        sum += v[i+15];
    }
    t = hpctimer_wtime() - t;

    printf("Sum = %d\n", sum);
    printf("Elapsed time (sec.): %.6f\n", t);

    free(v);
    return 0;
}
