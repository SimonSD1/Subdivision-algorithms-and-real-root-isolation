#ifndef ISOLATION_H
#define ISOLATION_H

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include "evaluate.h"
#include "descartes.h"
#include "mult_x_plus_one.h"
#include <stdlib.h>

typedef struct solution{
    // solution : (a / 2^p, (a+1) / 2^p)
    int c;
    int k;
    int h;
} solution;

void div_by_x(fmpz_poly_t pol);

void isolation(fmpz_poly_t pol, int c, int k, solution* solutions ,int* nb_sol);

#endif