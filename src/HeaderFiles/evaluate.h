#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>

#ifndef EVALUATE_H
#define EVALUATE_H

void evaluate_0(fmpz_t result, fmpz_poly_t poly);
void evaluate_1(fmpz_t result, fmpz_poly_t poly);
void evaluate_half(fmpz_t result, fmpz_poly_t poly);

#endif