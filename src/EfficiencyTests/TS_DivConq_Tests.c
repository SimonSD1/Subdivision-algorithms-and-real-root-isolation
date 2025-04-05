#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/thread_pool.h>
#include <flint/thread_support.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/taylorShift_implem.h"

////////// clock() will sum the execution time of all threads ! instead, use gettimeofday()

void benchmark_taylor_shiftMultiThread(fmpz_t shift, slong maxLen, int fixedVariable, FILE *fileResults, int numThreads)
{
    fmpz_poly_t poly, result;
    fmpz_poly_init(result);
    fmpz_poly_init(poly);
    struct timespec start, finish;
    double *tabTps = malloc(sizeof(double) * (maxLen));

    for (slong i = 0; i < maxLen; i++)
    {
        readPolyDATA(poly, fixedVariable, i);
        // Multi-threaded
        flint_set_num_threads(numThreads);
        clock_gettime(CLOCK_MONOTONIC, &start);
        fmpz_poly_taylor_shift_divconquer(result, poly, shift);
        clock_gettime(CLOCK_MONOTONIC, &finish);
        tabTps[i] = (finish.tv_sec - start.tv_sec);
        tabTps[i] += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    }

    fprintTab(tabTps, maxLen, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}

void benchmark_DivConq_Flint(fmpz_t shift, slong maxLen, int fixedVariable, FILE *fileResults)
{
    fmpz_poly_t poly, result;
    fmpz_poly_init(result);
    fmpz_poly_init(poly);
    double *tabTps = malloc(sizeof(double) * (maxLen));

    for (slong i = 0; i < maxLen; i++)
    {
        readPolyDATA(poly, fixedVariable, i);

        flint_set_num_threads(1);
        clock_t begin = clock();
        fmpz_poly_taylor_shift_divconquer(result, poly, shift);
        clock_t end = clock();
        tabTps[i] = (double)(end - begin); // / CLOCKS_PER_SEC;
    }

    fprintTab(tabTps, maxLen, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}

void benchmark_DivConq_Implem(slong maxLen, int fixedVariable, FILE *fileResults)
{
    printf("no table//////////////\n");
    fmpz_poly_t poly, result;
    fmpz_poly_init(result);
    fmpz_poly_init(poly);
    double *tabTps = malloc(sizeof(double) * (maxLen));

    for (slong i = 0; i < maxLen; i++)
    {
        readPolyDATA(poly, fixedVariable, i);
        printf("bonjou");
        clock_t begin = clock();
        // poly_shift_plus_one(result, poly, shift);
        poly_shift_plus_one_Non_Precomputed(result, poly);
        clock_t end = clock();
        printf("temps total = %lf",(double)(end - begin));
        tabTps[i] = (double)(end - begin); // / CLOCKS_PER_SEC;
    }

    fprintTab(tabTps, maxLen, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}

void benchmark_DivConq_Implem_Table(slong maxLen, int fixedVariable, FILE *fileResults)
{
    printf("table//////////////\n");
    fmpz_poly_t poly, result;
    fmpz_poly_init(poly);
    fmpz_poly_init(result);
    double *tabTps = malloc(sizeof(double) * (maxLen));

    for (slong i = 0; i < maxLen; i++)
    {
        readPolyDATA(poly, fixedVariable, i);

        clock_t begin = clock();
        poly_shift_plus_one_Precomputed(result, poly);
        clock_t end = clock();
        tabTps[i] = (double)(end - begin); // / CLOCKS_PER_SEC;
    }

    fprintTab(tabTps, maxLen, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}

int main(int argc, char *argv[])
{
    slong maxLen = 101;
    fmpz_t shift;
    fmpz_init_set_si(shift, 1);

    load_precomputed_polynomials(15);

    mkdir("src/EfficiencyTests/Results/TS_DivConq", 0777);

    // to run or re-run the tests, pass argument : ./flintMultTest -r
    // if not it will only make the graphs as png files in the results folder

    // note : on effectue les tests qu'une fois par taille, on ne fait pas de moyenne sur plusieurs tests pour une même taille
    FILE *fileResults;
    fileResults = fopen("src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree_Threading.txt", "w");
    if (fileResults == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }

    /*fprintf(fileResults, "DivConq Taylor Shift time efficiencies (coeffSize = 1500)\n"); // Title of the plot
    fprintf(fileResults, "Flint Single-Threaded\n");                                     // labels of the plot
    benchmark_taylor_shiftMultiThread(shift, maxLen, 0, fileResults, 1);                 // datas
    fprintf(fileResults, "Flint Multi-Threaded (6 threads)\n");
    benchmark_taylor_shiftMultiThread(shift, maxLen, 0, fileResults, 6);
    fclose(fileResults);
    */

   /*
    fileResults = fopen("src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize_Threading.txt", "w");
    if (fileResults == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }
    fprintf(fileResults, "DivConq Taylor Shift time efficiencies (degree = 1500)\n"); // Title of the plot
    fprintf(fileResults, "Flint Single-Threaded\n");                                  // labels of the plot
    benchmark_taylor_shiftMultiThread(shift, maxLen, 1, fileResults, 1);              // datas
    fprintf(fileResults, "Flint Multi-Threaded (6 threads)\n");
    benchmark_taylor_shiftMultiThread(shift, maxLen, 1, fileResults, 6);
    fclose(fileResults);
    */

    fileResults = fopen("src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree.txt", "w");
    if (fileResults == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }
    fprintf(fileResults, "DivConq Taylor Shift time efficiencies (coeffSize = 10000)\n"); // Title of the plot
    fprintf(fileResults, "Flint\n");
    //benchmark_DivConq_Flint(shift, maxLen, 0, fileResults);
    fprintf(fileResults, "Implem without table\n");
    benchmark_DivConq_Implem(maxLen, 0, fileResults);
    fprintf(fileResults, "Implem with table\n");
    benchmark_DivConq_Implem_Table(maxLen, 0, fileResults);
    fclose(fileResults);

    /*
    fileResults = fopen("src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize.txt", "w");
    if (fileResults == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }
    fprintf(fileResults, "DivConq Taylor Shift time efficiencies (degree = 1500)\n"); // Title of the plot
    fprintf(fileResults, "Flint\n");
    benchmark_DivConq_Flint(shift, maxLen, 1, fileResults);
    fprintf(fileResults, "Implem without table\n");
    benchmark_DivConq_Implem(maxLen, 1, fileResults);
    fprintf(fileResults, "Implem with table\n");
    benchmark_DivConq_Implem_Table(maxLen, 1, fileResults);
    fclose(fileResults);
    */

    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree_Threading.txt 'time' 0");
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize_Threading.txt 'time' 1");
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree.txt 'time' 0");
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize.txt 'time' 1");

    fmpz_clear(shift);
    free_global_precomputed();

    return 0;
}
