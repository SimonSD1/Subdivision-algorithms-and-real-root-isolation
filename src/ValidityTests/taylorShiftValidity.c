#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/taylorShift_implem.h"
#include "../HeaderFiles/functionsForTests.h"


int main()
{
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    fmpz_t shift;
    fmpz_init(shift);
    fmpz_set_si(shift, 1);

    fmpz_poly_t TrueResult;
    fmpz_poly_init(TrueResult);

    fmpz_poly_t result;
    fmpz_poly_init(result);

    readPolyDATA(poly, 0, 80);

    fmpz_poly_taylor_shift_divconquer(TrueResult, poly, shift);


    slong threshold = 512;  // Block size threshold
    fmpz_poly_t* power_array = NULL;
    slong levels = compute_power_array(&power_array, poly, threshold);
    flint_printf("Precomputed %ld levels of binomial coefficients\n", levels);

    iterative_taylor_shift_precompute(result, poly, threshold, power_array);

    for (int i = 0; i < levels; i++)
        fmpz_poly_clear(power_array[i]);
    free(power_array);


    if(fmpz_poly_equal(result, TrueResult))
        printf("Implem of iterative DivConq with precomputation table is valid :)\n");
    else
        printf("Implem of iterative DivConq with precomputation table is incorrect :(\n");

    /*naiveShift(result, poly, shift);
    if(fmpz_poly_equal(result, TrueResult))
        printf("Implem of naive Taylor shift is valid :)\n");
    else
        printf("Implem of naive Taylor shift is incorrect :(\n");*/

    fmpz_poly_clear(TrueResult);
    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    fmpz_clear(shift);
    return 0;
}
