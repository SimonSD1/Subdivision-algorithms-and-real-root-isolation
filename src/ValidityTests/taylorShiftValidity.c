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
    load_precomputed_polynomials(15);


    fmpz_poly_taylor_shift_divconquer(TrueResult, poly, shift);

    //poly_shift_plus_one(result, poly, shift);
    if(fmpz_poly_equal(result, TrueResult))
        printf("Implem of divide and conquer is valid :)\n");
    else
        printf("Implem of divide and conquer is incorrect :(\n");

    
    //poly_shift_plus_one_Precomputed(result, poly);
    if(fmpz_poly_equal(result, TrueResult))
        printf("Implem of DivConq with precomputation table is valid :)\n");
    else
        printf("Implem of DivConq with precomputation table is incorrect :(\n");

    /*naiveShift(result, poly, shift);
    if(fmpz_poly_equal(result, TrueResult))
        printf("Implem of naive Taylor shift is valid :)\n");
    else
        printf("Implem of naive Taylor shift is incorrect :(\n");*/

    fmpz_poly_clear(TrueResult);
    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    fmpz_clear(shift);
    free_global_precomputed();
    return 0;
}
