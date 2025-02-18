#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>


#include <string.h>


void mult_x_plus_1_power(fmpz_poly_t result, fmpz_poly_t poly, slong power){
    // void fmpz_bin_uiui(fmpz_t f, ulong n, ulong k)
    // Sets f to the binomial coefficient (n,k)

    // void _fmpz_poly_shift_left(fmpz *res, const fmpz *poly, slong len, slong n)
    // Sets (res, len + n) to (poly, len) shifted left by coefficients.

    fmpz_poly_init(result);

    fmpz_poly_t temp;
    fmpz_poly_init(temp);

    fmpz_t coef_binom;
    fmpz_init(coef_binom);

    for(slong i=0; i<power+1; i++){
        fmpz_poly_shift_left(temp,poly,i);

        fmpz_bin_uiui(coef_binom,power,i);

        printf("coeff : \n");
        fmpz_print(coef_binom);

        fmpz_poly_scalar_mul_fmpz(temp,temp,coef_binom);

        fmpz_poly_add(result,result,temp);
    }
    


}