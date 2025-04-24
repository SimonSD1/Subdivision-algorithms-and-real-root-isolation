#include "../HeaderFiles/isolation.h"
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include "../HeaderFiles/taylorShift_implem.h"

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
void var_change_to_inf(fmpz_poly_t result, fmpz_poly_t pol, fmpz_poly_t* power_array, slong threshold)
{
  // fmpz_poly_init(result);

  // taylor shift (result(x)=pol(x+1))
  // on renverse
  // f(1/x) = ( f_0 *x^d + f_1*x^(d-1) ) / x^d
  // pour compter les changement de signe on a pas besoin du denominateur
  int degre = fmpz_poly_degree(pol);

  // reverse et shift

  fmpz_poly_reverse(result, pol, degre + 1);
  iterative_taylor_shift_precompute(result, result, threshold,power_array);
}




// pol(x/2) * 2^d to avoid denominator
void split_left(fmpz_poly_t result, const fmpz_poly_t pol)
{
  // for each coeff a_i, we compute a_i * 2^(d-i)
  slong d = fmpz_poly_degree(pol);

  fmpz_poly_init(result);
  fmpz_poly_fit_length(result, d + 1);

  fmpz_t coeff;
  fmpz_init(coeff);

  for (slong i = 0; i <= d; i++)
  {
    fmpz_poly_get_coeff_fmpz(coeff, pol, i);
    fmpz_mul_2exp(coeff, coeff, d - i);
    fmpz_poly_set_coeff_fmpz(result, i, coeff);
  }

  fmpz_clear(coeff);
}

// pol((1+x)/2) * 2^d to avoid denominator
void split_right(fmpz_poly_t result, const fmpz_poly_t pol)
{
  // for each coeff a_i, we compute a_i * 2^(d-i)
  slong d = fmpz_poly_degree(pol);

  fmpz_poly_init(result);
  fmpz_poly_fit_length(result, d + 1);

  fmpz_t coeff;
  fmpz_init(coeff);

  for (slong i = 0; i <= d; i++)
  {
    fmpz_poly_get_coeff_fmpz(coeff, pol, i);
    fmpz_mul_2exp(coeff, coeff, d - i);
    fmpz_poly_set_coeff_fmpz(result, i, coeff);
  }

  fmpz_poly_taylor_shift_divconquer(result, result, FMPZ_ONE);
  fmpz_clear(FMPZ_ONE);
  fmpz_clear(coeff);
}




void isolation_recursive(fmpz_poly_t pol, fmpz_t c, slong k, solution *solutions, slong *nb_sol, fmpz_t temp,fmpz_poly_t *power_array, slong threshold )
{
  fmpz_t eval;
  fmpz_init(eval);

  fmpz_poly_t var_changed;
  fmpz_poly_init(var_changed);

  evaluate_0(eval, pol);
  

  int exact_root = 0;

  if (fmpz_is_zero(eval))
  {
    printf("div par x\n");
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
    fmpz_set(temp, c);
    fmpz_add_si(temp, temp, 1);
    fmpz_set(solutions[*nb_sol].c, temp);
    solutions[*nb_sol].k = k + 1;
    solutions[*nb_sol].is_exact = 1;
    (*nb_sol)++;
    exact_root = 1;
  }

  var_change_to_inf(var_changed, pol, power_array, threshold);

  int sign_changes = descartes_rule(var_changed);

  if (sign_changes == 0)
  {
    // nothing in this interval
    return;
  }

  else if (sign_changes == 1)
  {
    if(!exact_root)
    {
      // exactly one solution in this interval
      printf("ajoute interval : c=");
      fmpz_print(c);
      printf(" k=%ld\n", k);
      fmpz_set(solutions[*nb_sol].c, c);
      solutions[*nb_sol].k = k;
      solutions[*nb_sol].is_exact = 0;
      (*nb_sol)++;
    }
    return;
  }

  else
  {
    // at most 2 or more solution in this interval, we use bisection
    fmpz_poly_t pol_left;
    fmpz_poly_t pol_right;

    split_left(pol_left, pol);
    split_right(pol_right, pol);

    fmpz_t c_transformed;
    fmpz_mul_2exp(c_transformed, c, 1);

    isolation_recursive(pol_left, c_transformed, k + 1, solutions, nb_sol, temp, power_array,threshold);

    fmpz_add_ui(c_transformed, c_transformed, 1);

    isolation_recursive(pol_right, c_transformed, k + 1, solutions, nb_sol, temp,power_array,threshold);
  }

  fmpz_clear(eval);
  fmpz_poly_clear(var_changed);
}



/*
  compute pol(x*2^exp)
*/
void compose_mult_2exp(fmpz_poly_t result, fmpz_poly_t pol, slong exp)
{

  fmpz_t coeff;
  fmpz_init(coeff);

  slong degree = fmpz_poly_degree(pol);

  fmpz_poly_init2(result, degree + 1);

  fmpz_poly_get_coeff_fmpz(coeff, pol, 0);
  fmpz_poly_set_coeff_fmpz(result, 0, coeff);

  for (int i = 1; i < degree + 1; i++)
  {
    fmpz_poly_get_coeff_fmpz(coeff, pol, i);
    fmpz_mul_2exp(coeff, coeff, exp * (i));
    fmpz_poly_set_coeff_fmpz(result, i, coeff);
  }

  fmpz_clear(coeff);
}




void isolation_pos(fmpz_poly_t pol, solution *solutions, slong *nb_sol, slong *upper_power_of_two)
{
  // on fait fait le changement de variable pour avoir les racines sur [0,1]
  fmpz_t root_upper_bound;
  fmpz_init(root_upper_bound);

  fmpz_poly_t compressed_pol;
  fmpz_poly_init(compressed_pol);

  fmpz_poly_t *power_array;

  slong threshold=2;

  compute_power_array(&power_array, pol, threshold);

  local_max_bound_implementation(root_upper_bound, pol);

  *upper_power_of_two = fmpz_bits(root_upper_bound);

  printf("\nupper power of two = %ld\n", *upper_power_of_two);

  compose_mult_2exp(compressed_pol, pol, *upper_power_of_two);

  fmpz_t zero;
  fmpz_init_set_ui(zero, 0);

  fmpz_t temp;
  fmpz_init(temp);

  
  isolation_recursive(compressed_pol, zero, 0, solutions, nb_sol, temp, power_array, threshold);



  fmpz_poly_clear(compressed_pol);
  fmpz_clear(temp);
  fmpz_clear(root_upper_bound);
}




void poly_moins_x(fmpz_poly_t res, const fmpz_poly_t poly) {
  fmpz_poly_set(res, poly); // Copie le polynôme
  for (slong i = 1; i <= fmpz_poly_degree(res); i += 2) {
      fmpz_neg(fmpz_poly_get_coeff_ptr(res, i), fmpz_poly_get_coeff_ptr(res, i));
  }
}



void isolation(fmpz_poly_t pol, solution **solutions, slong *nb_sol, slong *nb_neg_sol, slong *upper_power_of_two_pos, slong *upper_power_of_two_neg)
{
  slong max_nb_roots;
  max_nb_roots = descartes_rule(pol);

  fmpz_poly_t pol_neg;
  fmpz_poly_init(pol_neg);
  poly_moins_x(pol_neg, pol);

  max_nb_roots += descartes_rule(pol_neg);

  *solutions = malloc(max_nb_roots * sizeof(solution));
  *nb_sol = 0;

  isolation_pos(pol_neg, *solutions, nb_sol, upper_power_of_two_neg);
  for(int i=0; i<*nb_sol; i++) {
    (*solutions)[i].sign = 0;
  }
  *nb_neg_sol = *nb_sol;

  isolation_pos(pol, *solutions, nb_sol, upper_power_of_two_pos);
  for(int i=*nb_neg_sol; i<*nb_sol; i++) {
    (*solutions)[i].sign = 1;
  }


  fmpz_poly_clear(pol_neg);
}