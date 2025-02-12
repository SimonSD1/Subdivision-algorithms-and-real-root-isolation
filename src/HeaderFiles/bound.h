#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>


#ifndef BOUND_H
#define BOUND_H

void Lagrange_bound(fmpz_t bound, fmpz_poly_t poly );
void Cauchy_bound(fmpz_t bound, const fmpz_poly_t poly);
void local_max_bound_implementation(fmpz_t ub, fmpz_poly_t poly);

#endif