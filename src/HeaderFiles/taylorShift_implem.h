#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>

#ifndef TAYLOR_SHIFT_IMPLEM_H
#define TAYLOR_SHIFT_IMPLEM_H

extern fmpz **global_precomputed;  // Array of fmpz vectors
extern slong *global_precomputed_len; // Array of vector lengths
extern slong global_precomputed_size;


void poly_shift_plus_one(fmpz_poly_t g, fmpz_poly_t poly, fmpz_t a);
void load_precomputed_polynomials(slong max_m);
void free_global_precomputed();
void poly_shift_plus_one_Precomputed(fmpz_poly_t g, const fmpz_poly_t poly);
void poly_shift_plus_one_Non_Precomputed(fmpz_poly_t result, const fmpz_poly_t poly);
void naiveShift(fmpz_poly_t result, fmpz_poly_t poly, fmpz_t a);

#endif