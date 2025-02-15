#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/evaluate.h"

void evaluate_0(fmpz_t result, fmpz_poly_t poly)
{
    fmpz_poly_get_coeff_fmpz(result, poly, 0);
}

void evaluate_1(fmpz_t result, fmpz_poly_t poly)
{
    fmpz_init(result);
    fmpz_t temp;
    fmpz_init(temp);

    for (int i = 0; i < poly->length; i++)
    {
        fmpz_poly_get_coeff_fmpz(temp, poly, i);
        fmpz_add(result, result, temp);
    }
    fmpz_clear(temp);
}

void evaluate_half(fmpq_t result, fmpz_poly_t poly)
{
    fmpq_init(result);
    
    fmpz_t numerator;
    fmpz_init(numerator);
    fmpz_zero(numerator); 

    fmpz_t denominator;
    fmpz_init(denominator);

    fmpz_t temp;
    fmpz_init(temp);

    slong d = poly->length - 1; 

    for (slong i = d; i >= 0; i--)
    {
        fmpz_poly_get_coeff_fmpz(temp, poly, i);
        fmpz_mul_2exp(temp, temp, d - i); 
        fmpz_add(numerator, numerator, temp);
    }

    fmpz_set_d_2exp(denominator, 1, d); 

    fmpq_set_fmpz_frac(result, numerator, denominator);

    
    fmpz_clear(numerator);
    fmpz_clear(denominator);
    fmpz_clear(temp);
}