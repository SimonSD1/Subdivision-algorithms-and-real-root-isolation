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


void benchmark_taylor_shiftMultiThread(fmpz_t shift, slong maxPow, int fixedVariable, FILE* fileResults, int numThreads)
{
    fmpz_poly_t poly, result;
    fmpz_poly_init(result);
    fmpz_poly_init(poly);
    struct timespec start, finish;
    double* tabTps = malloc(sizeof(double)*(15));

    for(slong i=0; i<maxPow; i++) {
        readPolyDATA(poly, fixedVariable, i+1, 1);
        // Multi-threaded
        flint_set_num_threads(numThreads);
        clock_gettime(CLOCK_MONOTONIC, &start);
        fmpz_poly_taylor_shift_divconquer(result, poly, shift);
        clock_gettime(CLOCK_MONOTONIC, &finish);
        tabTps[i] = (finish.tv_sec - start.tv_sec);
        tabTps[i] += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    }
    
    fprintTab(tabTps, 15, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}


void benchmark_DivConq_Flint(fmpz_t shift, slong maxPow, int fixedVariable, FILE* fileResults) {
    fmpz_poly_t poly, result;
    fmpz_poly_init(result);
    fmpz_poly_init(poly);
    double* tabTps = malloc(sizeof(double)*(15));

    for(slong i=0; i<maxPow; i++) {
        readPolyDATA(poly, fixedVariable, i+1, 1);

        flint_set_num_threads(1);
        clock_t begin = clock();
        fmpz_poly_taylor_shift_divconquer(result, poly, shift);
        clock_t end = clock();
        tabTps[i] = (double)(end - begin);// / CLOCKS_PER_SEC;

    }
    
    fprintTab(tabTps, 15, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}



void benchmark_DivConq_Implem(fmpz_t shift, slong maxPow, int fixedVariable, FILE* fileResults) {
    fmpz_poly_t poly, result;
    fmpz_poly_init(result);
    fmpz_poly_init(poly);
    double* tabTps = malloc(sizeof(double)*(15));

    for(slong i=0; i<maxPow; i++) {
        readPolyDATA(poly, fixedVariable, i+1, 1);

        clock_t begin = clock();
        poly_shift_plus_one(result, poly, shift, 129);
        clock_t end = clock();
        tabTps[i] = (double)(end - begin);// / CLOCKS_PER_SEC;

    }
    
    fprintTab(tabTps, 15, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}


void benchmark_DivConq_Implem_Table(slong maxPow, int fixedVariable, FILE* fileResults) {
    fmpz_poly_t poly, result;
    fmpz_poly_init(poly);
    fmpz_poly_init(result);
    double* tabTps = malloc(sizeof(double)*(15));

    for(slong i=0; i<maxPow; i++) {
        readPolyDATA(poly, fixedVariable, i+1, 1);
        
        clock_t begin = clock();
        poly_shift_plus_one_Precomputed2(result, poly, 129);
        clock_t end = clock();
        tabTps[i] = (double)(end - begin);// / CLOCKS_PER_SEC;

    }
    
    fprintTab(tabTps, 15, fileResults);

    fmpz_poly_clear(result);
    fmpz_poly_clear(poly);
    free(tabTps);
}


int main(int argc, char* argv[])
{
    slong maxPow = 15;
    fmpz_t shift;
    fmpz_init_set_si(shift, 1);

    load_precomputed_polynomials(15);

    mkdir("EfficiencyTests/Results/TS_DivConq", 0777);

    // to run or re-run the tests, pass argument : ./flintMultTest -r
    // if not it will only make the graphs as png files in the results folder
    if(argc > 1 && strcmp(argv[1], "-r") == 0) {
        //note : on effectue les tests qu'une fois par taille, on ne fait pas de moyenne sur plusieurs tests pour une même taille
        FILE* fileResults;
        fileResults = fopen("EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree_Threading.txt", "w");
        if (fileResults == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        fprintf(fileResults, "DivConq Taylor Shift time efficiencies (coeffSize = 511)\n");    //Title of the plot
        fprintf(fileResults, "Flint Single-Threaded\n");    //labels of the plot
        benchmark_taylor_shiftMultiThread(shift, maxPow, 0, fileResults, 1);      //datas
        fprintf(fileResults, "Flint Multi-Threaded (6 threads)\n");
        benchmark_taylor_shiftMultiThread(shift, maxPow, 0, fileResults, 6);
        fclose(fileResults);


        fileResults = fopen("EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize_Threading.txt", "w");
        if (fileResults == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        fprintf(fileResults, "DivConq Taylor Shift time efficiencies (degree = 512)\n");    //Title of the plot
        fprintf(fileResults, "Flint Single-Threaded\n");    //labels of the plot
        benchmark_taylor_shiftMultiThread(shift, maxPow, 1, fileResults, 1);      //datas
        fprintf(fileResults, "Flint Multi-Threaded (6 threads)\n");
        benchmark_taylor_shiftMultiThread(shift, maxPow, 1, fileResults, 6);
        fclose(fileResults);



        fileResults = fopen("EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize.txt", "w");
        if (fileResults == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        fprintf(fileResults, "DivConq Taylor Shift time efficiencies (degree = 512)\n");    //Title of the plot
        fprintf(fileResults, "Flint\n");
        benchmark_DivConq_Flint(shift, maxPow, 1, fileResults);
        fprintf(fileResults, "Implem without table\n");
        benchmark_DivConq_Implem(shift, maxPow, 1, fileResults);
        fprintf(fileResults, "Implem with table\n");
        benchmark_DivConq_Implem_Table(maxPow, 1, fileResults);
        fclose(fileResults);

        fileResults = fopen("EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree.txt", "w");
        if (fileResults == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        fprintf(fileResults, "DivConq Taylor Shift time efficiencies (coeffSize = 511)\n");    //Title of the plot
        fprintf(fileResults, "Flint\n");
        benchmark_DivConq_Flint(shift, maxPow, 0, fileResults);
        fprintf(fileResults, "Implem without table\n");
        benchmark_DivConq_Implem(shift, maxPow, 0, fileResults);
        fprintf(fileResults, "Implem with table\n");
        benchmark_DivConq_Implem_Table(maxPow, 0, fileResults);
        fclose(fileResults);

    }

    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree_Threading.txt 'time' 0 'exp'");
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize_Threading.txt 'time' 1 'exp'");
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingDegree.txt 'time' 0 'exp'");
    system("python3 EfficiencyTests/plotGenerator.py EfficiencyTests/Results/TS_DivConq/TS_DivConq_ChangingCoeffSize.txt 'time' 1 'exp'");

    fmpz_clear(shift);
    for (slong i = 0; i < global_precomputed_size; i++) {
        fmpz_poly_clear(global_precomputed[i]);
    }
    flint_free(global_precomputed);

    return 0;
}
