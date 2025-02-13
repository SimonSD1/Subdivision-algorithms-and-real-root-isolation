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

void comparShift(FILE *fileResults)
{
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    clock_t start, end;
    double *tabTimesDegreedivide_conquer_implem = malloc(sizeof(double) * (101)); // car il y a 101 polynômes dans DATA
    double *tabTimesDegreedivide_conquer_flint = malloc(sizeof(double) * (101));
    double *tabTimesDegreeNaive = malloc(sizeof(double) * (101));

    fmpz_poly_t result;
    fmpz_poly_init(result);

    // test for changing degree
    int fixedVariable = 0;
    fmpz_t a_fmpz;
    fmpz_init_set_si(a_fmpz, 50);
    printf("debut");
    for (slong i = 0; i <= 25; i++)
    {
        readPolyDATA(poly, fixedVariable, i);

        printf("debut");
        start = clock();
        poly_shift(result, poly, a_fmpz);
        end = clock();
        tabTimesDegreedivide_conquer_implem[i] = (double)(end - start);
        printf("poly shift");

        start = clock();
        fmpz_poly_taylor_shift_divconquer(result, poly, a_fmpz);
        end = clock();
        tabTimesDegreedivide_conquer_flint[i] = (double)(end - start);
        printf("flint");
        start = clock();

        start = clock();
        naiveShift(result, poly, a_fmpz);
        end = clock();
        tabTimesDegreeNaive[i] = (double)(end - start);
        printf("naive");
        start = clock();
    }

    fprintf(fileResults, "Taylor shift efficiency test\n"); // Title of the plot
    fprintf(fileResults, "Implem divide consuer\n");                                          // labels of the plot
    fprintTab(tabTimesDegreedivide_conquer_implem, 101, fileResults);
    fprintf(fileResults, "Flint div conquer\n");
    fprintTab(tabTimesDegreedivide_conquer_flint, 101, fileResults);
    fprintf(fileResults, "Naive\n");
    fprintTab(tabTimesDegreeNaive, 101, fileResults);
    fclose(fileResults);

    free(tabTimesDegreedivide_conquer_flint);
    free(tabTimesDegreedivide_conquer_implem);
    free(tabTimesDegreeNaive);

    fmpz_poly_clear(poly);
    fmpz_poly_clear(result);
}

int main(int argc, char *argv[])
{
    // to run or re-run the tests, pass argument : -runTests
    // if not it will only make the graphs as png files in the results folder
    printf("debut");
    if (argc > 1 && strcmp(argv[1], "-runTests") == 0)
    {
        printf("running");
        FILE *fileResults;
        fileResults = fopen("EfficiencyTests/Results/shiftChangingDegree.txt", "w");
        if (fileResults == NULL)
        {
            printf("The file is not opened. The program will "
                   "now exit.\n");
            exit(0);
        }

        comparShift(fileResults);
        system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/shiftChangingDegree.txt 'time' 0");

    }

    return 0;
}