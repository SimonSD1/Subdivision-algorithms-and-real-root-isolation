#include "../HeaderFiles/isolation.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "../HeaderFiles/functionsForTests.h"

#define MAX_ROOTS 10000


double* read_maple_roots(const char* filename, int* count) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        printf("Failed to open %s\n", filename);
        exit(1);
    }
    double* roots = malloc(sizeof(double) * MAX_ROOTS);
    *count = 0;
    while (fscanf(f, "%lf", &roots[*count]) == 1) {
        (*count)++;
        if (*count >= MAX_ROOTS) break;
    }
    fclose(f);
    return roots;
}


int run_maple_script() {
    // Essaie sur Linux natif
    int ret = system("\"/mnt/c/Program Files/Maple 2024/bin.X86_64_WINDOWS/cmaple.exe\"  src/ValidityTests/isolationValidityScript.mpl");
    if (ret == 0) return 0;
    return ret;
}




int main(void)
{
    /////////////////// Execute the testIsolation.mw first //////////////////////

    // global variable
    fmpz_init_set_ui(FMPZ_ONE, 1);

    fmpz_poly_t pol;
    fmpz_poly_init(pol);

    if (run_maple_script() != 0) {
        printf("Erreur lors de l'exécution du script Maple.\n");
        return 1;
    }


    FILE* filePoly;
    char pathFile[50];
    sprintf(pathFile, "src/bin/test.out");
    filePoly = fopen(pathFile, "r");
    if(!fmpz_poly_fread(filePoly, pol))
        printf("Could not read the FILE");

    fclose(filePoly);
    //fmpz_poly_set_str(pol, "4  -4 8 -5 1");

    //fmpz_poly_print_pretty(pol, "x");

    slong nb_sol, nb_neg_sol;
    solution *solutions = NULL;

    // bound used to compressed the polynomial
    // needed to get the solution from the solution on [0,1]
    slong upper_power_of_two_pos, upper_power_of_two_neg;

    fmpz_t temp;
    fmpz_init(temp);

    isolation(pol, &solutions, &nb_sol, &nb_neg_sol, &upper_power_of_two_pos, &upper_power_of_two_neg);
    printf("Isolated solutions: %ld //////// ", nb_sol);

    int nb_maple_roots = 0;
    double* maple_roots = read_maple_roots("src/bin/maple_roots.txt", &nb_maple_roots);
    printf("Maple found %d roots\n", nb_maple_roots);

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
            printf("Solution in interval: [ %lf , %lf ]", -upper, -lower);
            double root = maple_roots[nb_neg_sol-1 - i];
            if (root >= -upper && root <= -lower)
                printf(" -> Maple root %lf is inside this interval yayyy\n", root);
            else
                printf(" -> No Maple root found inside this interval :(\n");
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
            printf("Solution in interval: [ %lf , %lf ]", lower, upper);
            double root = maple_roots[i];
            if (root >= lower && root <= upper)
                printf(" -> Maple root %lf is inside this interval yayyy\n", root);
            else
                printf(" -> No Maple root found inside this interval :(\n");
        }

        fmpz_clear(solutions[i].c);
    }

    free(solutions);
    fmpz_poly_clear(pol);
    return 0;
}
