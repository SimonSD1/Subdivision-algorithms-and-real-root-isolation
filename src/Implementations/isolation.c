#include "../HeaderFiles/isolation.h"
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include "../HeaderFiles/taylorShift_implem.h"
#include "../HeaderFiles/coeff_truncation.h"
#include "../HeaderFiles/evaluate.h"
#include "../HeaderFiles/descartes.h"
#include "../HeaderFiles/bound.h"

fmpz_t FMPZ_ONE;

void div_by_x(fmpz_poly_t pol)
{
  fmpz_poly_shift_right(pol, pol, 1);
}

int fmpz_poly_is_half_root(fmpz_poly_t pol)
{
  fmpz_t numerator, temp;
  fmpz_init(numerator);
  fmpz_zero(numerator);
  fmpz_init(temp);

  slong d = pol->length - 1;
  for (slong i = d; i >= 0; i--)
  {
    fmpz_poly_get_coeff_fmpz(temp, pol, i);
    fmpz_mul_2exp(temp, temp, d - i);
    fmpz_add(numerator, numerator, temp);
  }

  int res = fmpz_is_zero(numerator);
  fmpz_clear(numerator);
  fmpz_clear(temp);
  return res;
}

void var_change_to_inf(fmpz_poly_t result, fmpz_poly_t pol, fmpz_poly_t *power_array, slong block_len, slong levels)
{
  int degre = fmpz_poly_degree(pol);
  fmpz_poly_reverse(result, pol, degre + 1);
  iterative_taylor_shift_precompute(result, result, power_array, block_len, levels);
  // fmpz_poly_taylor_shift(result,result,FMPZ_ONE);
}

void compose_divide_2_in_place(fmpz_poly_t pol)
{
  slong degree = pol->length - 1;
  for (slong i = 0; i < degree; i++)
  {
    fmpz *coeff = &(pol->coeffs[i]);
    fmpz_mul_2exp(coeff, coeff, degree - i);
  }
}

void compose_mult_2exp_in_place(fmpz_poly_t pol, slong exp)
{
  slong degree = fmpz_poly_degree(pol);
  for (int i = 1; i <= degree; i++)
  {
    fmpz *coeff = &(pol->coeffs[i]);
    fmpz_mul_2exp(coeff, coeff, exp * i);
  }
}

void compose_div_by_2exp_in_place(fmpz_poly_t pol, slong exp)
{
  slong d = fmpz_poly_degree(pol);
  fmpz *coeff_ptr = fmpz_poly_get_coeff_ptr(pol, 0);
  fmpz_mul_2exp(coeff_ptr, coeff_ptr, d + exp);
  for (slong i = 1; i <= d; i++)
  {
    coeff_ptr = fmpz_poly_get_coeff_ptr(pol, i);
    fmpz_mul_2exp(coeff_ptr, coeff_ptr, d - i);
  }
}

void divide2exp_coeff_in_place(fmpz_poly_t pol, slong exp)
{
  slong d = fmpz_poly_degree(pol);
  for (slong i = 0; i <= d; i++)
  {
    fmpz *coeff_ptr = fmpz_poly_get_coeff_ptr(pol, i);
    fmpz_cdiv_q_2exp(coeff_ptr, coeff_ptr, exp);
  }
}

slong sign_changes_trunc(fmpz_poly_t poly)
{
  slong degree = poly->length - 1;
  slong sign_changes = 0, last_sign = 0;

  fmpz_t precision_limit;
  fmpz_init(precision_limit);
  fmpz_set_si(precision_limit, 1);
  fmpz_mul_2exp(precision_limit, precision_limit, degree + 1);

  for (int i = 0; i <= degree; i++)
  {
    fmpz *coeff_poly = &(poly->coeffs[i]);
    if (fmpz_cmpabs(coeff_poly, precision_limit) > 0)
    {
      if (!fmpz_is_zero(coeff_poly))
      {
        int sign = fmpz_sgn(coeff_poly);
        if (last_sign != 0 && sign != last_sign)
          sign_changes++;
        last_sign = sign;
      }
    }
    if (sign_changes >= 2)
    {
      fmpz_clear(precision_limit);
      return sign_changes;
    }
  }

  fmpz_clear(precision_limit);
  return -1;
}

slong isolation_recursive(fmpz_poly_t pol, fmpz_t c, slong k, solution *solutions, slong *nb_sol, fmpz_t temp, fmpz_poly_t *power_array, slong block_len, slong levels, fmpz_poly_t poly_temp, int trunc_true)
{
  //fmpz_poly_set(poly_temp,pol);
  slong remove =_fmpz_poly_remove_content_2exp(pol->coeffs, pol->length);

  //printf("poly : ");
  //fmpz_poly_print_pretty(pol,"x");
  //printf("\n");
  //printf("poly temp : ");
  //fmpz_poly_print_pretty(poly_temp,"x");
  //printf("\n");
  //printf("\n");

  slong depth_right = 0;

  if (fmpz_is_zero(&(pol->coeffs[0])))
  {
    div_by_x(pol);
    fmpz_set(solutions[*nb_sol].c, c);
    solutions[*nb_sol].k = k;
    solutions[*nb_sol].is_exact = 1;
    (*nb_sol)++;
  }

  slong reliable_sign_changes = 0;
  slong max_bits = fmpz_poly_max_bits(pol);

  // here we truncate hopping that enough coefficient will be big reliable
  if (trunc_true)
  {
    if (max_bits < 0)
      max_bits = -max_bits;

    slong degree = pol->length - 1;
    slong remainding_bit_target = degree * 2;
    slong to_truncate = max_bits - remainding_bit_target;
    slong nb_sign_changes_trunc = -1;

    if (to_truncate > 0)
    {
      to_truncate = 0;
      truncate_coefficients(poly_temp, pol, to_truncate);
      //_fmpz_poly_remove_content_2exp(poly_temp->coeffs, poly_temp->length);

      if (!fmpz_is_zero(&(pol->coeffs[degree])))
      {
        var_change_to_inf(poly_temp, poly_temp, power_array, block_len, levels);
        nb_sign_changes_trunc = sign_changes_trunc(poly_temp);
      }
    }

    if (nb_sign_changes_trunc != -1)
    {
      reliable_sign_changes = nb_sign_changes_trunc;
    }

    else
    {
      var_change_to_inf(poly_temp, pol, power_array, block_len, levels);
      reliable_sign_changes = descartes_rule(poly_temp);
    }
  }
  else
  {
    var_change_to_inf(poly_temp, pol, power_array, block_len, levels);
    reliable_sign_changes = descartes_rule(poly_temp);
  }

  if (reliable_sign_changes == 0)
    return 0;
  else if (reliable_sign_changes == 1)
  {
    fmpz_set(solutions[*nb_sol].c, c);
    solutions[*nb_sol].k = k;
    solutions[*nb_sol].is_exact = 0;
    (*nb_sol)++;
    return 0;
  }

  //fmpz_t c_transformed;
  //fmpz_init(c_transformed);
  fmpz_mul_2exp(c, c, 1);

  compose_divide_2_in_place(pol);
  slong depth_left = isolation_recursive(pol, c, k + 1, solutions, nb_sol, temp, power_array, block_len, levels, poly_temp, trunc_true);

  fmpz_add_ui(c, c, 1);
  iterative_taylor_shift_precompute(pol, pol, power_array, block_len, levels);

  if (depth_left > 0)
  {
    compose_mult_2exp_in_place(pol, depth_left);
    divide2exp_coeff_in_place(pol, (pol->length - 1) * depth_left);
  }

  
  depth_right = isolation_recursive(pol, c, k + 1, solutions, nb_sol, temp, power_array, block_len, levels, poly_temp, trunc_true);
  
  fmpz_add_si(c,c,-1);
  fmpz_tdiv_q_2exp(c,c,1);
  //fmpz_clear(c_transformed);

  fmpz_poly_scalar_mul_2exp(pol,pol,remove);

  return 1 + depth_right;
}

void isolation_pos(fmpz_poly_t pol, solution *solutions, slong *nb_sol, slong *upper_power_of_two, fmpz_poly_t *power_array, slong block_len, slong levels, int trunc_true)
{
  fmpz_t root_upper_bound;
  fmpz_init(root_upper_bound);

  fmpz_poly_t compressed_pol, temp_poly;
  fmpz_poly_init(compressed_pol);
  fmpz_poly_init(temp_poly);

  fmpz_poly_bound_roots(root_upper_bound, pol);
  *upper_power_of_two = fmpz_bits(root_upper_bound);

  fmpz_poly_set(compressed_pol, pol);
  compose_mult_2exp_in_place(compressed_pol, *upper_power_of_two);

  fmpz_t zero, temp;
  fmpz_init_set_ui(zero, 0);
  fmpz_init(temp);
  isolation_recursive(compressed_pol, zero, 0, solutions, nb_sol, temp, power_array, block_len, levels, temp_poly, trunc_true);

  fmpz_poly_clear(compressed_pol);
  fmpz_poly_clear(temp_poly);
  fmpz_clear(temp);
  fmpz_clear(zero);
  fmpz_clear(root_upper_bound);
}

void poly_moins_x(fmpz_poly_t res, const fmpz_poly_t poly)
{
  fmpz_poly_set(res, poly);
  for (slong i = 1; i <= fmpz_poly_degree(res); i += 2)
  {
    fmpz_neg(fmpz_poly_get_coeff_ptr(res, i), fmpz_poly_get_coeff_ptr(res, i));
  }
}

void clear_power_array(fmpz_poly_t **p_array, slong num_elements_in_power_array)
{
  if (*p_array == NULL)
    return;
  for (slong i = 0; i < num_elements_in_power_array; i++)
  {
    fmpz_poly_clear((*p_array)[i]);
  }
  flint_free(*p_array);
  *p_array = NULL;
}

void isolation(fmpz_poly_t pol, solution **solutions, slong *nb_sol, slong *nb_neg_sol, slong *upper_power_of_two_pos, slong *upper_power_of_two_neg, int trunc_true)
{
  slong max_nb_roots = descartes_rule(pol);

  fmpz_poly_t pol_neg;
  fmpz_poly_init(pol_neg);
  poly_moins_x(pol_neg, pol);
  max_nb_roots += descartes_rule(pol_neg);

  *solutions = malloc(max_nb_roots * sizeof(solution));
  *nb_sol = 0;

  for (slong i = 0; i < max_nb_roots; i++)
    fmpz_init((*solutions)[i].c);

  fmpz_poly_t *power_array = NULL;
  slong block_len, levels;
  compute_power_array(&power_array, pol, &block_len, &levels);

  isolation_pos(pol_neg, *solutions, nb_sol, upper_power_of_two_neg, power_array, block_len, levels, trunc_true);
  *nb_neg_sol = *nb_sol;

  isolation_pos(pol, *solutions, nb_sol, upper_power_of_two_pos, power_array, block_len, levels, trunc_true);

  fmpz_poly_clear(pol_neg);
}
