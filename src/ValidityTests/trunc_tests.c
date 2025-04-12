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
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    fmpz_t shift;
    fmpz_init(shift);
    fmpz_set_si(shift, 1);

    fmpz_poly_t TrueResult;
    fmpz_poly_init(TrueResult);

    fmpz_poly_t trunc_poly;
    fmpz_poly_init(trunc_poly);
    fmpz_poly_t result;
    fmpz_poly_init(result);
    
    load_precomputed_polynomials(15);

    for (int i =5; i<=100; i+=10) {
        readPolyDATA(poly, 0, i);

        slong max_bit_before = fmpz_poly_max_bits(poly);
        //printf("max_bit_before = %ld\n", max_bit_before);
        fmpz_poly_taylor_shift_divconquer(TrueResult, poly, shift);
        slong max_bit_after = fmpz_poly_max_bits(TrueResult);
        //printf("max_bit_after = %ld\n\n", max_bit_after);
        slong bit_diff = labs(max_bit_after) - labs(max_bit_before);// labs => absolute value for long int
        printf("length (degree+1) = %ld\n", poly->length);
        printf("bit diff = %ld\n\n", bit_diff);
    }
    
    /*readPolyDATA(poly, 0, 100);
    truncate_coefficients(trunc_poly, poly, poly->length);

    clock_t begin = clock();
    fmpz_poly_taylor_shift_divconquer(TrueResult, poly, shift);
    clock_t end = clock();
    double tps = (double)(end - begin); // / CLOCKS_PER_SEC;
    printf("Tps taylor shift without truncation : %f\n", tps);

    begin = clock();
    fmpz_poly_taylor_shift_divconquer(result, trunc_poly, shift);
    end = clock();
    tps = (double)(end - begin);
    printf("Tps taylor shift with truncation : %f\n", tps);

    

    slong max_bit_before = fmpz_poly_max_bits(poly);
    printf("max_bit_before trunc = %ld\n", max_bit_before);

    slong max_bit_after = fmpz_poly_max_bits(trunc_poly);
    printf("max_bit_after trunc = %ld\n", max_bit_after);

    if (same_signs(trunc_poly, poly))
        printf("yaayy, polys have same signs after truncation\n"); 
    else
        printf("polys don't have same signs after truncation :(\n");


    if (same_signs(result, TrueResult))
        printf("yaayy, polys have same signs after taylor shift\n"); 
    else
        printf("polys don't have same signs after taylor shift :(\n");*/


    fmpz_poly_clear(TrueResult);
    fmpz_poly_clear(poly);
    fmpz_poly_clear(trunc_poly);
    fmpz_poly_clear(result);
    fmpz_clear(shift);
    free_global_precomputed();
    return 0;
}