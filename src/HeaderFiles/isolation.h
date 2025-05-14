#ifndef NEW_ISOLATION_RECURSIVE_H
#define NEW_ISOLATION_RECURSIVE_H

#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

typedef struct
{
    fmpz_t c;
    slong k;
    int is_exact;
} solution;

extern fmpz_t FMPZ_ONE;

void div_by_x(fmpz_poly_t pol);
int fmpz_poly_is_half_root(fmpz_poly_t pol);
void var_change_to_inf(fmpz_poly_t result, fmpz_poly_t pol, fmpz_poly_t *power_array, slong block_len, slong levels);
void split_left_in_place(fmpz_poly_t pol);
void split_right(fmpz_poly_t result, const fmpz_poly_t pol);
slong isolation_recursive(fmpz_poly_t pol, fmpz_t c, slong k, solution *solutions, slong *nb_sol, fmpz_t temp, fmpz_poly_t *power_array, slong block_len, slong levels, fmpz_poly_t poly_temp);
void compose_mult_2exp_in_place(fmpz_poly_t pol, slong exp);
void isolation_pos(fmpz_poly_t pol, solution *solutions, slong *nb_sol, slong *upper_power_of_two, fmpz_poly_t *power_array, slong block_len, slong levels);
void poly_moins_x(fmpz_poly_t res, const fmpz_poly_t poly);
void isolation(fmpz_poly_t pol, solution **solutions, slong *nb_sol, slong *nb_neg_sol, slong *upper_power_of_two_pos, slong *upper_power_of_two_neg);
void isolation_trunc(fmpz_poly_t pol, solution **solutions, slong *nb_sol, slong *nb_neg_sol, slong *upper_power_of_two_pos, slong *upper_power_of_two_neg);
slong sign_changes_trunc(fmpz_poly_t poly);
void isolation_recursive_trunc(fmpz_poly_t pol, fmpz_t c, slong k, solution *solutions, slong *nb_sol, fmpz_t temp, fmpz_poly_t *power_array, slong block_len, slong levels);
void isolation_pos_trunc(fmpz_poly_t pol, solution *solutions, slong *nb_sol, slong *upper_power_of_two, fmpz_poly_t *power_array, slong block_len, slong levels);

#endif