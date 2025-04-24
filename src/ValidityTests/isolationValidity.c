#include "../HeaderFiles/isolation.h"
#include <stdio.h>
#include "../HeaderFiles/functionsForTests.h"

int main(void)
{
    /////////////////// Execute the testIsolation.mw first //////////////////////

    // global variable
    fmpz_init_set_ui(FMPZ_ONE, 1);

    fmpz_poly_t pol;
    fmpz_poly_init(pol);

    FILE* filePoly;
    char pathFile[50];
    sprintf(pathFile, "src/bin/test.out");
    filePoly = fopen(pathFile, "r");
    if(!fmpz_poly_fread(filePoly, pol))
        printf("Could not read the FILE");

    fclose(filePoly);
    //fmpz_poly_set_str(pol, "4  -4 8 -5 1");

    fmpz_poly_print_pretty(pol, "x");

    slong nb_sol, nb_neg_sol;
    solution *solutions = NULL;

    // bound used to compressed the polynomial
    // needed to get the solution from the solution on [0,1]
    slong upper_power_of_two_pos, upper_power_of_two_neg;

    fmpz_t temp;
    fmpz_init(temp);

    isolation(pol, &solutions, &nb_sol, &nb_neg_sol, &upper_power_of_two_pos, &upper_power_of_two_neg);

    printf("Isolated solutions: %ld,\n", nb_sol);



    for (int i=nb_neg_sol-1; i>=0; i--) //we print negative roots intervals first
    {
        double denom = (double)(1LL << solutions[i].k);
        if (solutions[i].is_exact == 1)
        {
            fmpz_mul_2exp(temp,solutions[i].c, upper_power_of_two_neg);
            double exact_val = fmpz_get_d(temp) / denom;
            printf("Exact solution: %lf\n", -exact_val);
        }
        else
        {
            fmpz_mul_2exp(temp,solutions[i].c, upper_power_of_two_neg);
            double lower = fmpz_get_d(temp) / denom;

            //double lower = fmpz_get_d(solutions[i].c) / denom;
            fmpz_add(temp,solutions[i].c,FMPZ_ONE);
            fmpz_mul_2exp(temp,temp, upper_power_of_two_neg);
            double upper = (fmpz_get_d(temp)) / denom;
            printf("Solution in interval: [ %lf , %lf ]\n", -upper, -lower);
        }

        fmpz_clear(solutions[i].c);
    }


    for (int i = nb_neg_sol; i < nb_sol; i++)   //we print positive roots intervals
    {
        // Compute denominator as 2^k.
        double denom = (double)(1LL << solutions[i].k);
        if (solutions[i].is_exact == 1)
        {
            fmpz_mul_2exp(temp,solutions[i].c, upper_power_of_two_pos);
            double exact_val = fmpz_get_d(temp) / denom;
            printf("Exact solution: %lf\n", exact_val);
        }
        else
        {
            fmpz_mul_2exp(temp,solutions[i].c, upper_power_of_two_pos);
            double lower = fmpz_get_d(temp) / denom;

            //double lower = fmpz_get_d(solutions[i].c) / denom;
            fmpz_add(temp,solutions[i].c,FMPZ_ONE);
            fmpz_mul_2exp(temp,temp,upper_power_of_two_pos);
            double upper = (fmpz_get_d(temp)) / denom;
            printf("Solution in interval: [ %lf , %lf ]\n", lower, upper);
        }

        fmpz_clear(solutions[i].c);
    }

    free(solutions);
    fmpz_poly_clear(pol);
    return 0;
}
