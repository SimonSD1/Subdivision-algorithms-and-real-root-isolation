#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/thread_pool.h>
#include <flint/thread_support.h>
#include <flint/fmpq.h>
#include <flint/fmpz_vec.h>
#include <flint/ulong_extras.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/taylorShift_implem.h"

#define NB_RUNS 5


void benchmark_TS_DivConq(slong maxLen, int fixedVariable, FILE *fileResults)
{
    fmpz_t shift;
    fmpz_init_set_si(shift, 1);

    fmpz_poly_t poly, result;
    fmpz_poly_init(poly);
    fmpz_poly_init(result);
    double *tabTpsFlint = malloc(sizeof(double) * (maxLen));
    double *tabTpsImplem = malloc(sizeof(double) * (maxLen));
    double *tabTpsHorner = malloc(sizeof(double) * (maxLen));
    struct timespec start, end;
    fmpz_poly_t* power_array = NULL;
    slong block_len, levels;

    for (slong i = 0; i < maxLen; i++)
    {
        printf("poly iterative %ld\n", i);
        flint_cleanup();
        readPolyDATA(poly, fixedVariable, i);
        compute_power_array(&power_array, poly, &block_len, &levels);

        tabTpsImplem[i] = 0;
        //On effectue une moyenne sur NB_RUNS mesures
        for(int k=0; k<NB_RUNS; k++) {

            clock_gettime(CLOCK_MONOTONIC, &start);
            iterative_taylor_shift_precompute(result, poly, power_array, block_len, levels);
            clock_gettime(CLOCK_MONOTONIC, &end);
            tabTpsImplem[i] += (end.tv_sec - start.tv_sec) * 1e6 + (end.tv_nsec - start.tv_nsec) / 1e3;
        }
        tabTpsImplem[i] /= NB_RUNS;
        free_precompute_table(power_array, levels);
    }
    for (slong i = 0; i < maxLen; i++)
    {
        printf("poly flint %ld\n", i);
        flint_cleanup();
        readPolyDATA(poly, fixedVariable, i);

        tabTpsFlint[i] = 0;
        //On effectue une moyenne sur NB_RUNS mesures
        for(int k=0; k<NB_RUNS; k++) {
            flint_set_num_threads(1);

            clock_gettime(CLOCK_MONOTONIC, &start);
            fmpz_poly_taylor_shift_divconquer(result, poly, shift);
            clock_gettime(CLOCK_MONOTONIC, &end);
            tabTpsFlint[i] += (end.tv_sec - start.tv_sec) * 1e6 + (end.tv_nsec - start.tv_nsec) / 1e3;         
        }
        tabTpsFlint[i] /= NB_RUNS;
    }
    /*for (slong i = 0; i < maxLen; i++)
    {
        printf("poly Horner %ld\n", i);
        flint_cleanup();
        readPolyDATA(poly, fixedVariable, i);

        tabTpsFlint[i] = 0;
        //On effectue une moyenne sur NB_RUNS mesures
        for(int k=0; k<NB_RUNS; k++) {
            flint_set_num_threads(1);

            clock_gettime(CLOCK_MONOTONIC, &start);
            fmpz_poly_taylor_shift_horner(result, poly, shift);
            clock_gettime(CLOCK_MONOTONIC, &end);
            tabTpsHorner[i] += (end.tv_sec - start.tv_sec) * 1e6 + (end.tv_nsec - start.tv_nsec) / 1e3;         
        }
        tabTpsHorner[i] /= NB_RUNS;
    }*/


    // Title of the plot
    if(!fixedVariable)
        fprintf(fileResults, "DivConq Taylor Shift time efficiencies (coeffSize = 10000)\n");
    else
        fprintf(fileResults, "DivConq Taylor Shift time efficiencies (degree = 1500)\n");


    fprintf(fileResults, "Flint\n");
    fprintTab(tabTpsFlint, maxLen, fileResults);
    fprintf(fileResults, "Implem Iterative\n");
    fprintTab(tabTpsImplem, maxLen, fileResults);
    //fprintf(fileResults, "Horner (flint)\n");
    //fprintTab(tabTpsHorner, maxLen, fileResults);


    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    fmpz_clear(shift);
    free(tabTpsFlint);
    free(tabTpsImplem);
    free(tabTpsHorner);
}




int main(int argc, char *argv[])
{
    // to run or re-run the tests, pass argument : ./flintMultTest -r
    // if not it will only make the graphs as png files in the results folder

    if(argc > 1 && strcmp(argv[1], "-r") == 0) {            
        slong maxLen = 101;
        mkdir("src/EfficiencyTests/Results/TaylorShift", 0777);

        FILE *fileResults;
        fileResults = fopen("src/EfficiencyTests/Results/TaylorShift/TaylorShift_ChangingCoeffSize3.txt", "w");
        if (fileResults == NULL)
        {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        benchmark_TS_DivConq(maxLen, 1, fileResults);
        fclose(fileResults);


        
        fileResults = fopen("src/EfficiencyTests/Results/TaylorShift/TaylorShift_ChangingDegree3.txt", "w");
        if (fileResults == NULL)
        {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        benchmark_TS_DivConq(maxLen, 0, fileResults);
        fclose(fileResults);        
    }

    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/TaylorShift/TaylorShift_ChangingDegree3.txt 'time' 0");
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/TaylorShift/TaylorShift_ChangingCoeffSize3.txt 'time' 1");



    return 0;
}
