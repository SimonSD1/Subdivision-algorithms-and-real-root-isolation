#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>

#ifndef TAYLOR_SHIFT_IMPLEM_H
#define TAYLOR_SHIFT_IMPLEM_H


void divide_conquer(fmpz_poly_t g, const fmpz_poly_t f, fmpz_poly_t *precomputed, slong k);
void poly_shift(fmpz_poly_t g, fmpz_poly_t poly, fmpz_t a);
void naiveShift(fmpz_poly_t result, fmpz_poly_t poly, fmpz_t a);

#endif