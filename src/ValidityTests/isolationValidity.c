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
    
    slong nb_sol;
    solution *solutions = NULL;

    isolation(pol, &solutions, &nb_sol);

    printf("Isolated solutions: %ld,\n",nb_sol);
    
    for (int i = 0; i < nb_sol; i++)
    {
        // Compute denominator as 2^k.
        double denom = (double)(1LL << solutions[i].k);

        if (solutions[i].is_exact == 1)
        {
            double exact_val = fmpz_get_d(solutions[i].c) / denom;
            printf("Exact solution: %lf\n", exact_val);
        }
        else
        {
            double lower = fmpz_get_d(solutions[i].c) / denom;
            double upper = (fmpz_get_d(solutions[i].c) + 1.0) / denom;
            printf("Solution in interval: [ %lf , %lf ]\n", lower, upper);
        }

        fmpz_clear(solutions[i].c);
    }
    
    free(solutions);
    fmpz_poly_clear(pol);
    
    return 0;
}
