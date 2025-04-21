#ifndef ISOLATION_H
#define ISOLATION_H

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include "evaluate.h"
#include "descartes.h"
#include "mult_x_plus_one.h"
#include <stdlib.h>
#include "bound.h"

extern fmpz_t FMPZ_ONE;

typedef struct solution
{
  // solution : (c / 2^k, (c+1) / 2^k)
  // is exact =1 if the solution is exactly c/2^k
  fmpz_t c;
  slong k;
  int is_exact;
} solution;

void div_by_x(fmpz_poly_t pol);

void isolation_recursive(fmpz_poly_t pol, fmpz_t c, slong k, solution *solutions, slong *nb_sol, fmpz_t temp,fmpz_poly_t *power_array, slong threshold );

void isolation(fmpz_poly_t pol, solution **solutions, slong *nb_sol, slong *root_upper_bound);

void compose_mult_2exp(fmpz_poly_t result, fmpz_poly_t pol, slong exp);

#endif
