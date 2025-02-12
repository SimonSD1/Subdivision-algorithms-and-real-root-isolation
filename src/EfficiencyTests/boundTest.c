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

void comparBounds(FILE *fileResults)
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

    fprintf(fileResults, "Bound_efficiency_test\n"); // Title of the plot
    fprintf(fileResults, "Chauchy\n");  
    fprintFmpzTab(tabBoundDegreeCauchy, 101, fileResults);
    fprintf(fileResults, "Lagrange\n");
    fprintFmpzTab(tabBoundDegreeLagrange, 101, fileResults);
    fprintf(fileResults, "Modified Cauchy\n");
    fprintFmpzTab(tabBoundDegreeModifiedCauchy, 101, fileResults);
    fprintf(fileResults, "Flint\n");
    fprintFmpzTab(tabBoundDegreeFlint, 101, fileResults);

    //fprintf(fileResults, "Bound_efficiency_test\n"); // Title of the plot
    //fprintf(fileResults, "Chauchy\n");                                          // labels of the plot
    //fprintTab(tabTimesDegreeCauchy, 101, fileResults);
    //fprintf(fileResults, "Lagrange\n");
    //fprintTab(tabTimesDegreeLagrange, 101, fileResults);
    //fprintf(fileResults, "Modified Cauchy\n");
    //fprintTab(tabTimesDegreeModifiedCauchy, 101, fileResults);
    //fprintf(fileResults, "Flint\n");
    //fprintTab(tabTimesDegreeFlint, 101, fileResults);


    fclose(fileResults);

    ///fprintf(fileResults, "Taylor shift efficiency test\n"); // Title of the plot
    ///fprintf(fileResults, "Implem divide consuer\n");                                          // labels of the plot
    ///fprintTab(tabTimesDegreedivide_conquer_implem, 101, fileResults);
    ///fprintf(fileResults, "Flint div conquer\n");
    ///fprintTab(tabTimesDegreedivide_conquer_flint, 101, fileResults);
    ///fprintf(fileResults, "Naive\n");
    ///fprintTab(tabTimesDegreeNaive, 101, fileResults);
    ///fclose(fileResults);

    fmpz_poly_clear(poly);
}

int main(int argc, char *argv[])
{
    // to run or re-run the tests, pass argument : ./flintMultTest -runTests
    // if not it will only display the graphs
    if (argc > 1 && strcmp(argv[1], "-runTests") == 0)
    {
        printf("running");
        FILE *fileResults;
        fileResults = fopen("EfficiencyTests/Results/boundChangingDegree.txt", "w");
        if (fileResults == NULL)
        {
            printf("The file is not opened. The program will "
                   "now exit.\n");
            exit(0);
        }

        comparBounds(fileResults);
        system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/boundChangingDegree.txt 'time' 0");

    }

    return 0;
}