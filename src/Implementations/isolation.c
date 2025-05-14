#include "../HeaderFiles/isolation.h"
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include "../HeaderFiles/taylorShift_implem.h"
#include "../HeaderFiles/coeff_truncation.h"
#include "../HeaderFiles/evaluate.h"
#include "../HeaderFiles/descartes.h"

fmpz_t FMPZ_ONE;

void div_by_x(fmpz_poly_t pol)
{
  fmpz_poly_shift_right(pol, pol, 1);
}

int fmpz_poly_is_half_root(fmpz_poly_t pol)
{
  fmpz_t numerator;
  fmpz_init(numerator);
  fmpz_zero(numerator);

  fmpz_t temp;
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

// on calcul pol(1/(x+1))
void var_change_to_inf(fmpz_poly_t result, fmpz_poly_t pol, fmpz_poly_t *power_array, slong block_len, slong levels)
{
  // fmpz_poly_init(result);

  // taylor shift (result(x)=pol(x+1))
  // on renverse
  // f(1/x) = ( f_0 *x^d + f_1*x^(d-1) ) / x^d
  // pour compter les changement de signe on a pas besoin du denominateur
  int degre = fmpz_poly_degree(pol);

  // reverse et shift

  fmpz_poly_reverse(result, pol, degre + 1);
  iterative_taylor_shift_precompute(result, result, power_array, block_len, levels);
}

// on calcul pol(x/2) * 2^d
void compose_divide_2_in_place(fmpz_poly_t pol)
{
  slong degree = pol->length - 1;
  for (slong i = 0; i < degree; i++)
  {
    fmpz *coeff = &(pol->coeffs[i]);
    fmpz_mul_2exp(coeff, coeff, degree - i);
  }
}

// void mult_2exp_compose(fmpz_poly_t result, const fmpz_poly_t pol, slong exp)
//{
//   // for each coeff a_i, we compute a_i * 2^(d-i)
//   slong d = fmpz_poly_degree(pol);
//
//   fmpz_poly_init(result);
//   fmpz_poly_fit_length(result, d + 1);
//
//   fmpz_t coeff;
//   fmpz_init(coeff);
//
//   for (slong i = 1; i <= d; i++)
//   {
//     fmpz_poly_get_coeff_fmpz(coeff, pol, i);
//     fmpz_mul_2exp(coeff, coeff, exp * i);
//     fmpz_poly_set_coeff_fmpz(result, i, coeff);
//   }
//
//   fmpz_clear(coeff);
// }

/*
  compute pol(x*2^exp)
*/
void compose_mult_2exp_in_place(fmpz_poly_t pol, slong exp)
{

  fmpz *coeff;

  slong degree = fmpz_poly_degree(pol);

  for (int i = 1; i < degree + 1; i++)
  {
    coeff = &(pol->coeffs[i]);
    fmpz_mul_2exp(coeff, coeff, exp * (i));
  }
}

void compose_div_by_2exp_in_place(fmpz_poly_t pol, slong exp)
{
  slong d = fmpz_poly_degree(pol);
  fmpz *coeff_ptr = fmpz_poly_get_coeff_ptr(pol, 0);
  fmpz_mul_2exp(coeff_ptr, coeff_ptr, d + exp);
  for (slong i = 1; i <= d; i++)
  {
    fmpz *coeff_ptr = fmpz_poly_get_coeff_ptr(pol, i);
    fmpz_mul_2exp(coeff_ptr, coeff_ptr, d - i);
  }
}

void divide2exp_coeff_in_place(fmpz_poly_t pol, slong exp){
  slong d = fmpz_poly_degree(pol);
  for (slong i = 0; i <= d; i++)
  {
    fmpz *coeff_ptr = fmpz_poly_get_coeff_ptr(pol, i);
    fmpz_cdiv_q_2exp(coeff_ptr, coeff_ptr, exp);
  }
}

slong isolation_recursive(fmpz_poly_t pol, fmpz_t c, slong k, solution *solutions, slong *nb_sol, fmpz_t temp, fmpz_poly_t *power_array, slong block_len, slong levels, fmpz_poly_t poly_temp)
{


  //printf("dans recursion pol = ");
  //fmpz_poly_print_pretty(pol, "x");
  //printf("\n");

  int exact_root = 0;

  slong depth_left, depth_right;

  if (fmpz_is_zero(&(pol->coeffs[0])))
  {
    //printf("div par x\n");
    // if 0 is solution we can divide by x
    div_by_x(pol);
    // add the solution to the list
    fmpz_set(solutions[*nb_sol].c, c);
    solutions[*nb_sol].k = k;
    solutions[*nb_sol].is_exact = 1;
    (*nb_sol)++;
    exact_root = 1;
  }

  // on peut diviser par (x- 1/2)
  if (fmpz_poly_is_half_root(pol))
  {
   // printf("root /2\n");
    fmpz_set(temp, c);
    fmpz_mul_2exp(temp, temp, 1);
    fmpz_add_si(temp, temp, 1);
    fmpz_set(solutions[*nb_sol].c, temp);
    solutions[*nb_sol].k = k + 1;
    solutions[*nb_sol].is_exact = 1;
    (*nb_sol)++;
    exact_root = 1;
  }

  var_change_to_inf(poly_temp, pol, power_array, block_len, levels);

  int sign_changes = descartes_rule(poly_temp);

  if (sign_changes == 0)
  {
    // nothing in this interval
    //printf("rien ici\n");
    return 0;
  }

  else if (sign_changes == 1)
  {
    // printf("1 sign change\n");
    if (!exact_root)
    {
      // exactly one solution in this interval
      //printf("ajoute interval : c=");
      //fmpz_print(c);
      //printf(" k=%ld\n", k);
      fmpz_set(solutions[*nb_sol].c, c);
      solutions[*nb_sol].k = k;
      solutions[*nb_sol].is_exact = 0;
      (*nb_sol)++;
    }
    return 0;
  }

  else
  {
    // at most 2 or more solution in this interval, we use bisection
    // printf("on bisect");
    fmpz_t c_transformed;
    fmpz_init(c_transformed);
    fmpz_mul_2exp(c_transformed, c, 1);

    compose_divide_2_in_place(pol);

    //printf("left = ");
    //fmpz_poly_print_pretty(pol, "x");
    //printf("\n");
    depth_left = isolation_recursive(pol, c_transformed, k + 1, solutions, nb_sol, temp, power_array, block_len, levels,poly_temp);

    // printf("depth lef = %ld\n", depth_left);
    fmpz_add_ui(c_transformed, c_transformed, 1);

    //(pol, pol, power_array, block_len, levels);

    //printf("avant shift= ");
    //fmpz_poly_print_pretty(pol, "x");
    //printf("\n");
    fmpz_t one;
    fmpz_init_set_si(one, 1);
    fmpz_poly_taylor_shift(pol, pol, one);

    //printf("apres shift= ");
    //fmpz_poly_print_pretty(pol, "x");
    //printf("\n");
    if (depth_left > 0)
    {
      //printf("depth>0, = %ld\n",depth_left);
      //compose_div_by_2exp_in_place(pol, depth_left);
      compose_mult_2exp_in_place(pol,depth_left);
      // on doit encore diviser par 2**depth_left
      divide2exp_coeff_in_place(pol,(pol->length-1)*depth_left);
    }

    //printf("right = ");
    //fmpz_poly_print_pretty(pol, "x");
    //printf("\n");
    depth_right = isolation_recursive(pol, c_transformed, k + 1, solutions, nb_sol, temp, power_array, block_len, levels,poly_temp);

    // printf("depth right = %ld\n", depth_right);

    fmpz_clear(c_transformed);
  }
  
  return 1+depth_right;
}

void isolation_pos(fmpz_poly_t pol, solution *solutions, slong *nb_sol, slong *upper_power_of_two, fmpz_poly_t *power_array, slong block_len, slong levels)
{
  // on fait fait le changement de variable pour avoir les racines sur [0,1]
  fmpz_t root_upper_bound;
  fmpz_init(root_upper_bound);

  fmpz_poly_t compressed_pol;
  fmpz_poly_init(compressed_pol);

  fmpz_poly_t temp_poly;
  fmpz_poly_init(temp_poly);

  printf("bound\n");
  //fmpz_poly_print_pretty(pol, "x");
  printf("\n");

  fmpz_poly_bound_roots(root_upper_bound, pol);

  printf("bound = ");
  fmpz_print(root_upper_bound);
  printf("\n");

  *upper_power_of_two = fmpz_bits(root_upper_bound);

  printf("\nupper power of two = %ld\n", *upper_power_of_two);

  fmpz_poly_set(compressed_pol, pol);
  compose_mult_2exp_in_place(compressed_pol, *upper_power_of_two);

  printf("polynome compresse = ");
  //fmpz_poly_print_pretty(compressed_pol, "x");
  printf("\n");

  fmpz_t zero;
  fmpz_init_set_ui(zero, 0);

  fmpz_t temp;
  fmpz_init(temp);

  isolation_recursive(compressed_pol, zero, 0, solutions, nb_sol, temp, power_array, block_len, levels, temp_poly);

  fmpz_poly_clear(compressed_pol);
  fmpz_clear(temp);
  fmpz_clear(root_upper_bound);
}

void poly_moins_x(fmpz_poly_t res, const fmpz_poly_t poly)
{
  fmpz_poly_set(res, poly); // Copie le polynôme
  for (slong i = 1; i <= fmpz_poly_degree(res); i += 2)
  {
    fmpz_neg(fmpz_poly_get_coeff_ptr(res, i), fmpz_poly_get_coeff_ptr(res, i));
  }
}

void isolation(fmpz_poly_t pol, solution **solutions, slong *nb_sol, slong *nb_neg_sol, slong *upper_power_of_two_pos, slong *upper_power_of_two_neg)
{
  slong max_nb_roots;
  max_nb_roots = descartes_rule(pol);

  printf("max nb roots = %ld\n", max_nb_roots);

  fmpz_poly_t pol_neg;
  fmpz_poly_init(pol_neg);
  poly_moins_x(pol_neg, pol);

  max_nb_roots += descartes_rule(pol_neg);

  *solutions = malloc(max_nb_roots * sizeof(solution));
  *nb_sol = 0;

  fmpz_poly_t *power_array = NULL;
  slong block_len, levels;
  compute_power_array(&power_array, pol, &block_len, &levels);

  isolation_pos(pol_neg, *solutions, nb_sol, upper_power_of_two_neg, power_array, block_len, levels);
  *nb_neg_sol = *nb_sol;

  printf("fini les negatif\n");
  printf("passe a pol = ");
  //fmpz_poly_print_pretty(pol, "x");
  printf("\n");
  isolation_pos(pol, *solutions, nb_sol, upper_power_of_two_pos, power_array, block_len, levels);

  fmpz_poly_clear(pol_neg);
}
