#include "../HeaderFiles/isolation.h"
#include <stdio.h>

int main(void)
{
    // global variable
    fmpz_init_set_ui(FMPZ_ONE, 1);

    fmpz_poly_t pol;
    fmpz_poly_init(pol);

    // Define the polynomial as in SageMath (*840 to avoid division issues)
    fmpz_poly_set_str(pol, "6  840 -26616 26546 -9253 1354 -71");

    fmpz_poly_print_pretty(pol, "x");

    slong nb_sol;
    solution *solutions = NULL;

    // bound used to compressed the polynomial
    // needed to get the solution from the solution on [0,1]
    slong upper_power_two;

    fmpz_t temp;
    fmpz_init(temp);

    isolation(pol, &solutions, &nb_sol, &upper_power_two);

    printf("Isolated solutions: %ld,\n", nb_sol);

    for (int i = 0; i < nb_sol; i++)
    {
        // Compute denominator as 2^k.
        double denom = (double)(1LL << solutions[i].k);

        if (solutions[i].is_exact == 1)
        {
            fmpz_mul_2exp(temp,solutions[i].c,upper_power_two);
            double exact_val = fmpz_get_d(temp) / denom;
            printf("Exact solution: %lf\n", exact_val);
        }
        else
        {
            fmpz_mul_2exp(temp,solutions[i].c,upper_power_two);
            double lower = fmpz_get_d(temp) / denom;

            //double lower = fmpz_get_d(solutions[i].c) / denom;
            fmpz_add(temp,solutions[i].c,FMPZ_ONE);
            fmpz_mul_2exp(temp,temp,upper_power_two);
            double upper = (fmpz_get_d(temp)) / denom;
            printf("Solution in interval: [ %lf , %lf ]\n", lower, upper);
        }

        fmpz_clear(solutions[i].c);
    }

    free(solutions);
    fmpz_poly_clear(pol);

    return 0;
}
