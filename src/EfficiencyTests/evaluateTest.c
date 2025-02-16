#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/evaluate.h"

// 0 compar with changing degree
void comparEvaluate(FILE *fileResults, int fixedVariable)
{
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    fmpq_t half;
    fmpq_init(half);
    fmpq_set_si(half, 1, 2);

    clock_t start, end;
    double *tabImplem = calloc(101, sizeof(double)); // car il y a 101 polynômes dans DATA
    double *tabFlint = calloc(101, sizeof(double));

    fmpq_t resultFlint;
    fmpq_init(resultFlint);

    fmpq_t resultImplem;
    fmpq_init(resultImplem);

    for (slong i = 0; i <= 100; i++)
    {
        readPolyDATA(poly, fixedVariable, i);

        start = clock();
        evaluate_half(resultImplem, poly);
        end = clock();
        tabImplem[i] = (double)(end - start);

        start = clock();
        fmpz_poly_evaluate_fmpq(resultFlint, poly, half);
        end = clock();
        tabFlint[i] = (double)(end - start);

        if (!fmpq_equal(resultFlint, resultImplem))
        {
            
            printf("error, bad result");

            fmpq_print(resultFlint);
            printf("\n");
            fmpq_print(resultImplem);
            return;
        }
    }

    fprintf(fileResults, "Evaluation test\n"); // Title of the plot
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
    
    FILE *fileResultsDegree = fopen("EfficiencyTests/Results/EvaluateDegree.txt", "w");
    if (fileResultsDegree == NULL)
    {
        printf("Failed to open EvaluateDegree.txt.\n");
        exit(0);
    }
    comparEvaluate(fileResultsDegree, 0);
    fclose(fileResultsDegree); 
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/EvaluateDegree.txt 'time' 0");

    FILE *fileResultsCoeffsSize = fopen("EfficiencyTests/Results/EvaluateCoeffs.txt", "w");
    if (fileResultsCoeffsSize == NULL)
    {
        printf("Failed to open EvaluateCoeffs.txt.\n");
        exit(0);
    }
    comparEvaluate(fileResultsCoeffsSize, 1);
    fclose(fileResultsCoeffsSize); 
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/EvaluateCoeffs.txt 'time' 1");

    return 0;
}