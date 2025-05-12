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

    readPolyDATA(poly, 1, 100);

    fmpz_poly_taylor_shift_divconquer(TrueResult, poly, shift);

    fmpz_poly_t* power_array = NULL;
    slong block_len, levels;
    compute_power_array(&power_array, poly, &block_len, &levels);
    flint_printf("Precomputed %ld levels of binomial coefficients\n", levels);

    iterative_taylor_shift_precompute(result, poly, power_array, block_len, levels);

    if(fmpz_poly_equal(result, TrueResult))
        printf("Implem of iterative DivConq with precomputation table is valid :)\n");
    else
        printf("Implem of iterative DivConq with precomputation table is incorrect :(\n");

    /*naiveShift(result, poly, shift);
    if(fmpz_poly_equal(result, TrueResult))
        printf("Implem of naive Taylor shift is valid :)\n");
    else
        printf("Implem of naive Taylor shift is incorrect :(\n");*/

    free_precompute_table(power_array, levels);
    fmpz_poly_clear(TrueResult);
    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    fmpz_clear(shift);
    return 0;
}
