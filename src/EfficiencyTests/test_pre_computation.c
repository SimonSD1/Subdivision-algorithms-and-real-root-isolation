#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/mult_x_plus_one.h"

// 0 compar with changing degree
void compar_x_plus_1(FILE *fileResults)
{
    

    clock_t start, end;
    double *tabImplem = calloc(101, sizeof(double)); // car il y a 101 polynômes dans DATA
    double *tabFlint = calloc(101, sizeof(double));

    for (slong iter = 0; iter < 15; iter++)
    {

        fmpz_poly_t *precomputed = flint_malloc(iter * sizeof(fmpz_poly_t));

        for (slong i = 0; i < iter; i++)
        {
            fmpz_poly_init(precomputed[i]);
        }

        fmpz_poly_t *precomputed_flint = flint_malloc(iter * sizeof(fmpz_poly_t));

        for (slong i = 0; i < iter; i++)
        {
            fmpz_poly_init(precomputed_flint[i]);
        }

        start = clock();
        fmpz_poly_set_coeff_si(precomputed_flint[0], 1, 1);
        fmpz_poly_set_coeff_si(precomputed_flint[0], 0, 1);
        // square (x+a) m times by using the previous one
        for (slong i = 1; i < iter; i++)
        {
            // square the previous one in the i-th
            fmpz_poly_sqr(precomputed_flint[i], precomputed_flint[i - 1]);
        }

        end = clock();
        tabFlint[iter] = (double)(end - start);

        start = clock();
        fmpz_poly_set_coeff_si(precomputed[0], 1, 1);
        fmpz_poly_set_coeff_si(precomputed[0], 0, 1);
        // square (x+a) m times by using the previous one
        for (slong i = 1; i < iter; i++)
        {
            // square the previous one in the i-th
            mult_x_plus_1_power(precomputed[i], precomputed[i - 1], 1 << (i - 1));
        }
        end = clock();
        tabImplem[iter] = (double)(end - start);

        for (slong i = 0; i < iter; i++)
        {

            if (!fmpz_poly_equal(precomputed[i], precomputed_flint[i]))
            {
                printf("incorrect\n");
                break;
            }
        }
    }

    fprintf(fileResults, "Pre computation test\n"); // Title of the plot
    fprintf(fileResults, "Implem\n");               // labels of the plot
    fprintTab(tabImplem, 101, fileResults);
    fprintf(fileResults, "Flint\n");
    fprintTab(tabFlint, 101, fileResults);

    free(tabImplem);
    free(tabFlint);

}

int main()
{
    printf("running");

    FILE *fileResultsCoeffsSize = fopen("EfficiencyTests/Results/precompute.txt", "w");
    if (fileResultsCoeffsSize == NULL)
    {
        printf("Failed to open.\n");
        exit(0);
    }
    compar_x_plus_1(fileResultsCoeffsSize);
    fclose(fileResultsCoeffsSize);
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/precompute.txt 'time' 0");

    return 0;
}