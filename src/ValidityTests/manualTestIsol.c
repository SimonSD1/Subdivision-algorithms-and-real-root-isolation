#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/bound.h"
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/descartes.h"
#include "../HeaderFiles/isolation.h"

int main(int argc, char const *argv[])
{
    fmpz_poly_t poly;

    fmpz_init_set_ui(FMPZ_ONE, 1);

    fmpz_poly_init(poly);

    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, time(NULL), clock()); // Seed with current time & clock ticks

    //fmpz_poly_randtest(poly, state, 100, 10000);
    //fmpz_poly_set_coeff_fmpz(poly,0,FMPZ_ONE);

    //fmpz_poly_set_str(poly,"6  840 -26616 26546 -9253 1354 -71");


    fmpz_poly_set_str(poly,"3  3 -4 1");

    //-71*x^5 + 1354*x^4 - 9253*x^3 + 26546*x^2 - 26616*x + 840

    fmpz_poly_print_pretty(poly, "x");

    solution *solutions = NULL;

    slong nb_sol, nb_neg_sol, upper_power_of_two, upper_power_of_two_neg;

    isolation(poly, &solutions, &nb_sol, &nb_neg_sol, &upper_power_of_two, &upper_power_of_two_neg);

    printf("Isolated solutions: %ld //////// ", nb_sol);

    fmpz_t temp;
    fmpz_init(temp);

    for (int i = nb_neg_sol - 1; i >= 0; i--) // we print negative roots intervals first
    {
        double denom = (double)(1LL << solutions[i].k);
        if (solutions[i].is_exact == 1)
        {
            fmpz_mul_2exp(temp, solutions[i].c, upper_power_of_two_neg);
            double exact_val = fmpz_get_d(temp) / denom;
            printf("Exact solution: %lf\n", -exact_val);
        }
        else
        {
            fmpz_mul_2exp(temp, solutions[i].c, upper_power_of_two_neg);
            double lower = fmpz_get_d(temp) / denom;

            // double lower = fmpz_get_d(solutions[i].c) / denom;
            fmpz_add(temp, solutions[i].c, FMPZ_ONE);
            fmpz_mul_2exp(temp, temp, upper_power_of_two_neg);
            double upper = (fmpz_get_d(temp)) / denom;
            printf("Solution in interval: [ %lf , %lf ]", -upper, -lower);
        }

        fmpz_clear(solutions[i].c);
    }

    for (int i = nb_neg_sol; i < nb_sol; i++) // we print positive roots intervals
    {
        // Compute denominator as 2^k.
        double denom = (double)(1LL << solutions[i].k);
        if (solutions[i].is_exact == 1)
        {
            fmpz_mul_2exp(temp, solutions[i].c, upper_power_of_two);
            double exact_val = fmpz_get_d(temp) / denom;
            printf("Exact solution: %lf\n", exact_val);
        }
        else
        {
            fmpz_mul_2exp(temp, solutions[i].c, upper_power_of_two);
            double lower = fmpz_get_d(temp) / denom;

            // double lower = fmpz_get_d(solutions[i].c) / denom;
            fmpz_add(temp, solutions[i].c, FMPZ_ONE);
            fmpz_mul_2exp(temp, temp, upper_power_of_two);
            double upper = (fmpz_get_d(temp)) / denom;
            printf("Solution in interval: [ %lf , %lf ]", lower, upper);
        }

        fmpz_clear(solutions[i].c);
    }

    free(solutions);
    fmpz_poly_clear(poly);
    return 0;
}
