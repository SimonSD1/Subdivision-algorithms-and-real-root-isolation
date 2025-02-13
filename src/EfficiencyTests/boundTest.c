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

void comparBounds(FILE *fileResultsSpeed,FILE *fileResultsValues)
{
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    clock_t start, end;
    double *tabTimesDegreeCauchy = malloc(sizeof(double) * (101)); // car il y a 101 polynômes dans DATA
    double *tabTimesDegreeLagrange = malloc(sizeof(double) * (101));
    double *tabTimesDegreeModifiedCauchy = malloc(sizeof(double) * (101));
    double *tabTimesDegreeFlint = malloc(sizeof(double) * (101));

    fmpz *tabBoundDegreeCauchy = _fmpz_vec_init(101);
    fmpz *tabBoundDegreeLagrange = _fmpz_vec_init(101);
    fmpz *tabBoundDegreeModifiedCauchy = _fmpz_vec_init(101);
    fmpz *tabBoundDegreeFlint = _fmpz_vec_init(101);
    fmpz_t bound;

   

    // test for changing degree
    int fixedVariable = 0;
    for (slong i = 0; i <= 100; i++)
    {
        readPolyDATA(poly, fixedVariable, i);


        start = clock();
        Cauchy_bound(bound, poly);
        end = clock();
        tabTimesDegreeCauchy[i] = (double)(end - start);
        fmpz_set(&tabBoundDegreeCauchy[i], bound);

        start = clock();
        Lagrange_bound(bound, poly);
        end = clock();
        tabTimesDegreeLagrange[i] = (double)(end - start);
        fmpz_set(&tabBoundDegreeLagrange[i], bound);

        start = clock();
        local_max_bound_implementation(bound, poly);
        end = clock();
        tabTimesDegreeCauchy[i] = (double)(end - start);
        fmpz_set(&tabBoundDegreeModifiedCauchy[i], bound);

        start = clock();
        //  Fujiwara’s bound
        fmpz_poly_bound_roots(bound,poly);
        end = clock();
        tabTimesDegreeCauchy[i] = (double)(end - start);
        fmpz_set(&tabBoundDegreeFlint[i], bound);
    }

    fprintf(fileResultsValues, "Bound values test\n"); // Title of the plot
    fprintf(fileResultsValues, "Chauchy\n");  
    fprintFmpzTab(tabBoundDegreeCauchy, 101, fileResultsValues);
    fprintf(fileResultsValues, "Lagrange\n");
    fprintFmpzTab(tabBoundDegreeLagrange, 101, fileResultsValues);
    fprintf(fileResultsValues, "Modified Cauchy\n");
    fprintFmpzTab(tabBoundDegreeModifiedCauchy, 101, fileResultsValues);
    fprintf(fileResultsValues, "Flint\n");
    fprintFmpzTab(tabBoundDegreeFlint, 101, fileResultsValues);

    fprintf(fileResultsSpeed, "Bound compputing speed test\n"); // Title of the plot
    fprintf(fileResultsSpeed, "Chauchy\n");                                          // labels of the plot
    fprintTab(tabTimesDegreeCauchy, 101, fileResultsSpeed);
    fprintf(fileResultsSpeed, "Lagrange\n");
    fprintTab(tabTimesDegreeLagrange, 101, fileResultsSpeed);
    fprintf(fileResultsSpeed, "Modified Cauchy\n");
    fprintTab(tabTimesDegreeModifiedCauchy, 101, fileResultsSpeed);
    fprintf(fileResultsSpeed, "Flint\n");
    fprintTab(tabTimesDegreeFlint, 101, fileResultsSpeed);


    fclose(fileResultsSpeed);
    fclose(fileResultsValues);


    fmpz_poly_clear(poly);
}

int main(int argc, char *argv[])
{
    // to run or re-run the tests, pass argument : ./flintMultTest -runTests
    // if not it will only display the graphs
    if (argc > 1 && strcmp(argv[1], "-runTests") == 0)
    {
        printf("running");
        FILE *fileResultsSpeed;
        FILE *fileResultsValues;
        fileResultsSpeed = fopen("EfficiencyTests/Results/boundSpeedChangingDegree.txt", "w");
        fileResultsValues = fopen("EfficiencyTests/Results/boundValuesChangingDegree.txt", "w");

        if (fileResultsSpeed == NULL || fileResultsValues == NULL)
        {
            printf("The file is not opened. The program will "
                   "now exit.\n");
            exit(0);
        }

        comparBounds(fileResultsSpeed,fileResultsValues);
        system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/boundSpeedChangingDegree.txt 'time' 0");
        system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/boundValuesChangingDegree.txt 'time' 0");

    }

    return 0;
}