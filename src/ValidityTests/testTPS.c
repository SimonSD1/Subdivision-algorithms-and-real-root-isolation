#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include "../HeaderFiles/taylorShift_implem.h"
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/coeff_truncation.h"


int main()
{
    fmpz_poly_t poly, TrueResult;
    fmpz_poly_init(poly);
    fmpz_poly_init(TrueResult);
    
    readPolyDATA(poly, 0, 60);

    fmpz_poly_t* power_array;
    slong block_len, levels;
    compute_power_array(&power_array, poly, &block_len, &levels);
    
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    iterative_taylor_shift_precompute(TrueResult, poly, power_array, block_len, levels);
    clock_gettime(CLOCK_MONOTONIC, &end);
    double tps = (end.tv_sec - start.tv_sec) * 1e6 + (end.tv_nsec - start.tv_nsec) / 1e3;
    printf("Tps 1 : %f\n", tps);

    clock_gettime(CLOCK_MONOTONIC, &start);
    iterative_taylor_shift_precompute(TrueResult, poly, power_array, block_len, levels);
    clock_gettime(CLOCK_MONOTONIC, &end);
    tps = (end.tv_sec - start.tv_sec) * 1e6 + (end.tv_nsec - start.tv_nsec) / 1e3;
    printf("Tps 2 : %f\n", tps);



    //cleanup the precomputation
    for (int i = 0; i < levels; i++)
        fmpz_poly_clear(power_array[i]);
    free(power_array);

    fmpz_poly_clear(TrueResult);
    return 0;
}