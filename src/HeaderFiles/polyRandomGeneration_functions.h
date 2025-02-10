#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>

void fmpz_poly_randtest_dense(fmpz_poly_t poly, flint_rand_t state, slong degree, flint_bitcnt_t coeff_size);
// Génère un polynôme dense (pas de coefficient égal à 0) et à coefficients entiers
void writePolyFile(slong degree, flint_bitcnt_t coeff_size, int fixedVariable, int i);
// Ecrit un polynôme généré aléatoirement avec  fmpz_poly_randtest_dense dans un fichier