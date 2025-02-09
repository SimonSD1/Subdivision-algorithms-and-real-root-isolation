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

void evaluate_half(fmpz_t result, fmpz_poly_t poly)
{
    fmpz_init(result);
    fmpz_t temp;
    fmpz_init(temp);

    for (slong i = 0; i < poly->length; i++)
    {
        fmpz_poly_get_coeff_fmpz(temp, poly, i);
        fmpz_fdiv_q_2exp(temp,temp, i);
        fmpz_add(result, result, temp);
    }

    fmpz_clear(temp);
}