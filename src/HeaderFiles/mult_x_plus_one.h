#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>


#ifndef MULT_H
#define MULT_H

void mult_x_plus_1_power(fmpz_poly_t result, fmpz_poly_t poly, slong power);

#endif
