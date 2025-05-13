#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include "../HeaderFiles/taylorShift_implem.h"
#include "../HeaderFiles/functionsForTests.h"


#ifndef COEFF_TRUNCATION_H
#define COEFF_TRUNCATION_H

int same_signs(const fmpz_poly_t poly1, const fmpz_poly_t poly2);
void fmpz_trunc(fmpz_t rop, const fmpz_t op, slong keep_bits);
void truncate_coefficients(fmpz_poly_t result, const fmpz_poly_t poly, slong trunc);


#endif // COEFF_TRUNCATION_H