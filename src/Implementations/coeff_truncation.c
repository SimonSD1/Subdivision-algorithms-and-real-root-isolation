#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include "../HeaderFiles/taylorShift_implem.h"
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/coeff_truncation.h"

int same_signs(const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    slong len_1 = fmpz_poly_length(poly1);
    slong len_2 = fmpz_poly_length(poly2);

    // Vérification rapide : mêmes longueurs ?
    if (len_1 != len_2)
    {
        return 0;
    }

    fmpz_t coeff_1, coeff_2;
    fmpz_init(coeff_1);
    fmpz_init(coeff_2);

    for (slong i = 0; i < len_1; i++)
    {
        fmpz_poly_get_coeff_fmpz(coeff_1, poly1, i);
        fmpz_poly_get_coeff_fmpz(coeff_2, poly2, i);

        int sign_a = fmpz_sgn(coeff_1);
        int sign_b = fmpz_sgn(coeff_2);

        if (sign_a != sign_b)
        {
            fmpz_clear(coeff_1);
            fmpz_clear(coeff_2);
            return 0;
        }
    }

    fmpz_clear(coeff_1);
    fmpz_clear(coeff_2);
    return 1;
}


void truncate_coefficients(fmpz_poly_t result, const fmpz_poly_t poly, slong trunc)
{
    fmpz_poly_fit_length(result,poly->length);
    for(slong i=0; i<poly->length;i++){
        fmpz_tdiv_q_2exp(&(result->coeffs[i]),&(poly->coeffs[i]),trunc);
    }

    result->length=poly->length;
}