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
void comparEvaluate(FILE *fileResultsHalf,FILE * fileResults1, FILE * fileResults0, int fixedVariable)
{
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    fmpq_t half;
    fmpq_init(half);
    fmpq_set_si(half, 1, 2);

    fmpz_t zero;
    fmpz_init_set_si(zero,0);

    fmpz_t one;
    fmpz_init_set_si(one,1);

    clock_t start, end;
    double *tabImplemHalf = calloc(101, sizeof(double)); // car il y a 101 polynômes dans DATA
    double *tabFlintHalf = calloc(101, sizeof(double));

    double *tabImplem0 = calloc(101, sizeof(double)); // car il y a 101 polynômes dans DATA
    double *tabFlint0 = calloc(101, sizeof(double));

    double *tabImplem1 = calloc(101, sizeof(double)); // car il y a 101 polynômes dans DATA
    double *tabFlint1 = calloc(101, sizeof(double));

    fmpq_t resultFlint;
    fmpq_init(resultFlint);

    fmpq_t resultImplem;
    fmpq_init(resultImplem);

    fmpz_t resultImplemZ;
    fmpz_init(resultImplemZ);

    fmpz_t resultFlintZ;
    fmpz_init(resultFlintZ);

    for (slong i = 0; i <= 100; i++)
    {
        readPolyDATA(poly, fixedVariable, i, 0);


        // half
        start = clock();
        evaluate_half(resultImplem, poly);
        end = clock();
        tabImplemHalf[i] = (double)(end - start);

        start = clock();
        fmpz_poly_evaluate_fmpq(resultFlint, poly, half);
        end = clock();
        tabFlintHalf[i] = (double)(end - start);

        if (!fmpq_equal(resultFlint, resultImplem))
        {

            printf("error, bad result");

            fmpq_print(resultFlint);
            printf("\n");
            fmpq_print(resultImplem);
            return;
        }


        // zero
        start = clock();
        evaluate_0(resultImplemZ, poly);
        end = clock();
        tabImplem0[i] = (double)(end - start);

        start = clock();
        fmpz_poly_evaluate_fmpz(resultFlintZ, poly, zero);
        end = clock();
        tabFlint0[i] = (double)(end - start);

        if (!fmpz_equal(resultFlintZ, resultImplemZ))
        {
            printf("error, bad result");

            fmpz_print(resultFlintZ);
            printf("\n");
            fmpz_print(resultImplemZ);
            return;
        }


        // one
        start = clock();
        evaluate_1(resultImplemZ, poly);
        end = clock();
        tabImplem1[i] = (double)(end - start);

        start = clock();
        fmpz_poly_evaluate_fmpz(resultFlintZ, poly, one);
        end = clock();
        tabFlint1[i] = (double)(end - start);

        if (!fmpz_equal(resultFlintZ, resultImplemZ))
        {
            printf("error, bad result");

            fmpz_print(resultFlintZ);
            printf("\n");
            fmpz_print(resultImplemZ);
            return;
        }
    }

    fprintf(fileResultsHalf, "Evaluation test half\n"); // Title of the plot
    fprintf(fileResultsHalf, "Implem\n");          // labels of the plot
    fprintTab(tabImplemHalf, 101, fileResultsHalf);
    fprintf(fileResultsHalf, "Flint\n");
    fprintTab(tabFlintHalf, 101, fileResultsHalf);

    fprintf(fileResults1, "Evaluation test 1\n"); // Title of the plot
    fprintf(fileResults1, "Implem\n");          // labels of the plot
    fprintTab(tabImplem1, 101, fileResults1);
    fprintf(fileResults1, "Flint\n");
    fprintTab(tabFlint1, 101, fileResults1);

    fprintf(fileResults0, "Evaluation test 0\n"); // Title of the plot
    fprintf(fileResults0, "Implem\n");          // labels of the plot
    fprintTab(tabImplem0, 101, fileResults0);
    fprintf(fileResults0, "Flint\n");
    fprintTab(tabFlint0, 101, fileResults0);

    free(tabImplemHalf);
    free(tabFlintHalf);

    free(tabImplem0);
    free(tabFlint0);

    free(tabImplem1);
    free(tabFlint1);

    fmpz_poly_clear(poly);
}

int main()
{
    printf("running");

    FILE *fileResultsDegreeHalf = fopen("EfficiencyTests/Results/EvaluateHalfDegree.txt", "w");
    if (fileResultsDegreeHalf == NULL)
    {
        printf("Failed to open EvaluateHalfDegree.txt.\n");
        exit(0);
    }
    FILE *fileResultsDegree1 = fopen("EfficiencyTests/Results/Evaluate1Degree.txt", "w");
    if (fileResultsDegree1 == NULL)
    {
        printf("Failed to open Evaluate1Degree.txt.\n");
        exit(0);
    }
    FILE *fileResultsDegree0 = fopen("EfficiencyTests/Results/Evaluate0Degree.txt", "w");
    if (fileResultsDegree0 == NULL)
    {
        printf("Failed to open Evaluate0Degree.txt.\n");
        exit(0);
    }


    comparEvaluate(fileResultsDegreeHalf,fileResultsDegree0,fileResultsDegree1, 0);
    fclose(fileResultsDegreeHalf);
    fclose(fileResultsDegree0);
    fclose(fileResultsDegree1);
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/EvaluateHalfDegree.txt 'time' 0");
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/Evaluate0Degree.txt 'time' 0");
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/Evaluate1Degree.txt 'time' 0");


    
    FILE *fileResultsCoeffsHalf = fopen("EfficiencyTests/Results/EvaluateHalfCoeffs.txt", "w");
    if (fileResultsCoeffsHalf == NULL)
    {
        printf("Failed to open EvaluateHalfCoeffs.txt.\n");
        exit(0);
    }
    FILE *fileResultsCoeffs1 = fopen("EfficiencyTests/Results/Evaluate1Coeffs.txt", "w");
    if (fileResultsCoeffs1 == NULL)
    {
        printf("Failed to open Evaluate1Coeffs.txt.\n");
        exit(0);
    }
    FILE *fileResultsCoeffs0 = fopen("EfficiencyTests/Results/Evaluate0Coeffs.txt", "w");
    if (fileResultsCoeffs0 == NULL)
    {
        printf("Failed to open Evaluate0Coeffs.txt.\n");
        exit(0);
    }
    
    comparEvaluate(fileResultsCoeffsHalf,fileResultsCoeffs0,fileResultsCoeffs1,1);
    fclose(fileResultsCoeffsHalf);
    fclose(fileResultsCoeffs0);
    fclose(fileResultsCoeffs1);
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/EvaluateHalfCoeffs.txt 'time' 1");
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/Evaluate1Coeffs.txt 'time' 1");
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/Evaluate0Coeffs.txt 'time' 1");

    return 0;
}