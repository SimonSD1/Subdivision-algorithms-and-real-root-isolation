#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/bound.h"
#include "../HeaderFiles/taylorShift_implem.h"

void comparShift(FILE *fileResults, int fixedVariable)
{
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    clock_t start, end;
    double *tabTimesImplem = calloc(101, sizeof(double)); // car il y a 101 polynômes dans DATA
    double *tabTimesFlint = calloc(101, sizeof(double));
    double *tabTimesNaive = calloc(101, sizeof(double));
    double *tabTimesImplem50 = calloc(101, sizeof(double));

    fmpz_poly_t result_div;
    fmpz_poly_init(result_div);

    fmpz_poly_t result_flint;
    fmpz_poly_init(result_flint);

    fmpz_poly_t result_naive;
    fmpz_poly_init(result_naive);

    fmpz_t a_fmpz;
    fmpz_init_set_si(a_fmpz, 50);
    printf("debut");
    for (slong i = 0; i <= 20; i++)
    {
        readPolyDATA(poly, fixedVariable, i, 0);

        printf("debut");
        start = clock();
        poly_shift_plus_one(result_div, poly, a_fmpz,0);
        end = clock();
        tabTimesImplem[i] = (double)(end - start);
        printf("poly shift");

        printf("debut");
        start = clock();
        poly_shift_plus_one(result_div, poly, a_fmpz,50);
        end = clock();
        tabTimesImplem50[i] = (double)(end - start);
        printf("poly shift");

        start = clock();
        fmpz_poly_taylor_shift_divconquer(result_flint, poly, a_fmpz);
        end = clock();
        tabTimesFlint[i] = (double)(end - start);
        printf("flint");
        start = clock();

        start = clock();
        naiveShift(result_naive, poly, a_fmpz);
        end = clock();
        tabTimesNaive[i] = (double)(end - start);
        printf("naive");
        start = clock();

        if (!fmpz_poly_equal(result_div, result_flint) || !fmpz_poly_equal(result_naive, result_flint))
        {
            printf("error, bad result");
            return;
        }
    }

    if (fixedVariable == 0)
    {
        fprintf(fileResults, "Taylor shift efficiency test, changing degree\n"); // Title of the plot
    }
    else
    {
        fprintf(fileResults, "Taylor shift efficiency test, changing coeffs size\n"); // Title of the plot
    }
    fprintf(fileResults, "Implem divide conquer\n"); // labels of the plot
    fprintTab(tabTimesImplem, 101, fileResults);
    fprintf(fileResults, "Implem divide conquer 50\n"); // labels of the plot
    fprintTab(tabTimesImplem50, 101, fileResults);
    fprintf(fileResults, "Flint div conquer\n");
    fprintTab(tabTimesFlint, 101, fileResults);
    fprintf(fileResults, "Naive\n");
    fprintTab(tabTimesNaive, 101, fileResults);
    fclose(fileResults);

    free(tabTimesFlint);
    free(tabTimesImplem);
    free(tabTimesNaive);
    free(tabTimesImplem50);

    fmpz_poly_clear(poly);
    fmpz_poly_clear(result_div);
    fmpz_poly_clear(result_naive);
    fmpz_poly_clear(result_flint);
}

int main()
{
    // to run or re-run the tests, pass argument : -runTests
    // if not it will only make the graphs as png files in the results folder

    printf("running");
    FILE *fileResultsDegree;
    fileResultsDegree = fopen("EfficiencyTests/Results/shiftChangingDegree.txt", "w");
    if (fileResultsDegree == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }

    FILE *fileResultsCoeffsSize;
    fileResultsCoeffsSize = fopen("EfficiencyTests/Results/shiftChangingCoeffsSize.txt", "w");
    if (fileResultsCoeffsSize == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }

    comparShift(fileResultsDegree, 0);
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/shiftChangingDegree.txt 'time' 0");

    comparShift(fileResultsCoeffsSize, 1);
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/shiftChangingCoeffsSize.txt 'time' 1");

    return 0;
}