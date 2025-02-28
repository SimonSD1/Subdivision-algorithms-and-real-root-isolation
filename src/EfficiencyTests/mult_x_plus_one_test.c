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
void compar_x_plus_1(FILE *fileResults, int fixedVariable)
{
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    fmpz_poly_t x_plus_one;
    fmpz_poly_init(x_plus_one);

    fmpz_poly_t result;
    fmpz_poly_init(result);

    fmpz_poly_t result_flint;
    fmpz_poly_init(result_flint);

    fmpz_poly_set_str(x_plus_one, "2  1 1");

    fmpz_poly_set_str(poly, "1  1");

    printf("poly\n");
    fmpz_poly_print(poly);
    printf("\n");

    printf("+1\n");
    fmpz_poly_print(x_plus_one);
    printf("\n");

    slong power=100;

    fmpz_poly_pow_multinomial(x_plus_one,x_plus_one,power);

    clock_t start, end;
    double *tabImplem = calloc(101, sizeof(double)); // car il y a 101 polynômes dans DATA
    double *tabFlint = calloc(101, sizeof(double));

    for (slong i = 0; i <= 100; i++)
    {
        readPolyDATA(poly, fixedVariable, i);

        start = clock();
        mult_x_plus_1_power(result, poly, power);
        end = clock();
        tabImplem[i] = (double)(end - start);
        
        start = clock();
        fmpz_poly_mul(result_flint, poly, x_plus_one);
        end = clock();
        tabFlint[i] = (double)(end - start);
        


        if (!fmpz_poly_equal(result, result_flint))
        {
            printf("incorrect\n");
            break;
        }
    }

    fprintf(fileResults, "Mult test\n"); // Title of the plot
    fprintf(fileResults, "Implem\n");                         // labels of the plot
    fprintTab(tabImplem, 101, fileResults);
    fprintf(fileResults, "Flint\n");
    fprintTab(tabFlint, 101, fileResults);

    free(tabImplem);
    free(tabFlint);

    fmpz_poly_clear(poly);
}


int main()
{
    printf("running");

    FILE *fileResultsCoeffsSize = fopen("EfficiencyTests/Results/x_plus_1.txt", "w");
    if (fileResultsCoeffsSize == NULL)
    {
        printf("Failed to open.\n");
        exit(0);
    }
    compar_x_plus_1(fileResultsCoeffsSize, 1);
    fclose(fileResultsCoeffsSize);
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/x_plus_1.txt 'time' 0");

    return 0;
}