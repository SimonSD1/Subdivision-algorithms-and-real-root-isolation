#ifndef ISOLATION_H
#define ISOLATION_H

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include "evaluate.h"
#include "descartes.h"

typedef struct solution{
    // solution : (a / 2^p, (a+1) / 2^p)
    fmpz_t c;
    fmpz_t k;
    fmpz_t h;
} solution;

void div_by_x(fmpz_poly_t pol);

void isolation(fmpz_poly_t pol, fmpz_t c, __int32_t k, solution* solutions ,int* nb_sol);

#endif