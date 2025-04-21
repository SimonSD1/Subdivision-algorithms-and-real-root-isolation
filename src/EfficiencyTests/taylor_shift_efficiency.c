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
#include <math.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/taylorShift_implem.h"

void benchmark_DivConq_Flint(fmpz_t shift, slong maxLen, int fixedVariable, FILE *fileResults)
{
    fmpz_poly_t poly, result;
    fmpz_poly_init(result);
    fmpz_poly_init(poly);
    struct timespec start, finish;
    double *tabTps = malloc(sizeof(double) * (maxLen));

    for (slong i = 0; i < maxLen; i++)
    {
        readPolyDATA(poly, fixedVariable, i);

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

void benchmark_DivConq_Implem_Table(slong maxLen, int fixedVariable, FILE *fileResults)
{
    printf("table//////////////\n");
    fmpz_poly_t poly, result;
    fmpz_poly_init(poly);
    fmpz_poly_init(result);
    struct timespec start, finish;
    double *tabTps = malloc(sizeof(double) * (maxLen));

    for (slong i = 0; i < maxLen; i++)
    {
        readPolyDATA(poly, fixedVariable, i);

        clock_gettime(CLOCK_MONOTONIC, &start);
        poly_shift_plus_one_Precomputed(result, poly);
        clock_gettime(CLOCK_MONOTONIC, &finish);
        tabTps[i] = (finish.tv_sec - start.tv_sec);
        tabTps[i] += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    }

    fprintTab(tabTps, maxLen, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}

void benchmark_DecoupePrecompute(fmpz_t shift, slong maxLen, int fixedVariable, FILE *fileResults)
{
    printf("Decoupe with precompute//////////////\n");
    fmpz_poly_t poly, result;
    fmpz_poly_init(result);
    fmpz_poly_init(poly);
    struct timespec start_compute, finish_compute, start_shift, finish_shift;
    double *tabTps = malloc(sizeof(double) * (maxLen));

    

    for (slong i = 10; i < maxLen; i++)
    {
        

        readPolyDATA(poly, fixedVariable, i);

        slong threshold = 256 * 2; // Example threshold


        //clock_gettime(CLOCK_MONOTONIC, &start_compute);
        fmpz_poly_t *power_array;
        slong l = compute_power_array(&power_array, poly, threshold);
        //clock_gettime(CLOCK_MONOTONIC, &finish_compute);


        //double time_compute = (finish_compute.tv_sec - start_compute.tv_sec) + (finish_compute.tv_nsec - start_compute.tv_nsec) / 1000000000.0;
        //printf("Time to compute power array (length %ld): %lf seconds\n", poly->length, time_compute);

        clock_gettime(CLOCK_MONOTONIC, &start_shift);
        iterative_taylor_shift_precompute(result, poly, threshold, power_array);
        clock_gettime(CLOCK_MONOTONIC, &finish_shift);

        tabTps[i] = (finish_shift.tv_sec - start_shift.tv_sec);
        //tabTps[i] += (finish_shift.tv_nsec - start_shift.tv_nsec) / 1000000000.0;

        for (int i = 0; i < l; i++)
            fmpz_poly_clear(power_array[i]);
        free(power_array);

        fmpz_poly_clear(result); // Re-initialize for the next iteration
        fmpz_poly_init(result);
    }

    fprintTab(tabTps, maxLen, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}

int main(/*int argc, char *argv[]*/)
{
    slong maxLen = 101;
    fmpz_t shift;
    fmpz_init_set_si(shift, 1);

    load_precomputed_polynomials(5001);

    mkdir("src/EfficiencyTests/Results/TS_DivConq", 0777);

    FILE *fileResults;

    fileResults = fopen("src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree_Comparison.txt", "w");
    if (fileResults == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }
    fprintf(fileResults, "DivConq Taylor Shift time efficiencies (coeffSize = 10000)\n"); // Title of the plot
    fprintf(fileResults, "Flint\n");
    benchmark_DivConq_Flint(shift, maxLen, 0, fileResults);
    printf("\ndivconquer\n");
    fprintf(fileResults, "Implem with table\n");
    benchmark_DivConq_Implem_Table(maxLen, 0, fileResults);
    printf("\nimplem table\n");
    fprintf(fileResults, "Decoupe with Precompute (shift +1)\n");
    benchmark_DecoupePrecompute(shift, maxLen, 0, fileResults);
    printf("\ndecoupage\n");
    fclose(fileResults);

    fileResults = fopen("src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize_Comparison.txt", "w");
    if (fileResults == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }
    fprintf(fileResults, "DivConq Taylor Shift time efficiencies (degree = 1500)\n"); // Title of the plot
    fprintf(fileResults, "Flint\n");
    benchmark_DivConq_Flint(shift, maxLen, 1, fileResults);
    fprintf(fileResults, "Implem with table\n");
    benchmark_DivConq_Implem_Table(maxLen, 1, fileResults);
    fprintf(fileResults, "Decoupe with Precompute (shift +1)\n");
    benchmark_DecoupePrecompute(shift, maxLen, 1, fileResults);
    fclose(fileResults);

    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree_Comparison.txt 'time' 0");
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize_Comparison.txt 'time' 1");

    fmpz_clear(shift);
    free_global_precomputed();

    return 0;
}