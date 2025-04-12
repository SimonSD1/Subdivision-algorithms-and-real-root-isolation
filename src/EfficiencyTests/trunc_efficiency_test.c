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
#include "../HeaderFiles/coeff_truncation.h"

void benchmark_degree() {
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    fmpz_t shift;
    fmpz_init(shift);
    fmpz_set_si(shift, 1);

    fmpz_poly_t TrueResult;
    fmpz_poly_init(TrueResult);

    fmpz_poly_t trunc_poly;
    fmpz_poly_init(trunc_poly);
    fmpz_poly_t result;
    fmpz_poly_init(result);
    
    double *tabTps = malloc(sizeof(double) * 101);
    double *tabTps2 = malloc(sizeof(double) * 101);

    for (int i=0; i<=100; i++) {
        readPolyDATA(poly, 0, i);
        truncate_coefficients(trunc_poly, poly, poly->length);

        clock_t begin = clock();
        fmpz_poly_taylor_shift_divconquer(TrueResult, poly, shift);
        clock_t end = clock();
        tabTps[i] = (double)(end - begin); // / CLOCKS_PER_SEC;

        begin = clock();
        fmpz_poly_taylor_shift_divconquer(result, trunc_poly, shift);
        end = clock();
        tabTps2[i] = (double)(end - begin); // / CLOCKS_PER_SEC;

        if (same_signs(result, TrueResult))
            printf("yaayy, polys have same signs after taylor shift : %d\n", i); 
        else
            printf("polys don't have same signs after taylor shift :( : %d\n", i);
    }


    mkdir("src/EfficiencyTests/Results/Coeff_Truncation", 0777);
    FILE *fileResults;
    fileResults = fopen("src/EfficiencyTests/Results/Coeff_Truncation/Coeff_Truncation_ChangingDegree.txt", "w");
    if (fileResults == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }


    fprintf(fileResults, "Taylor Shift with truncated coefficient execution time (coeffSize = 10000)\n"); // Title of the plot
    fprintf(fileResults, "No truncation\n");
    fprintTab(tabTps, 101, fileResults);
    fprintf(fileResults, "With truncation\n");
    fprintTab(tabTps2, 101, fileResults);

    
    fclose(fileResults);
    fmpz_poly_clear(TrueResult);
    fmpz_poly_clear(poly);
    fmpz_poly_clear(trunc_poly);
    fmpz_poly_clear(result);
    fmpz_clear(shift);
}




void benchmark_coeffSize() {
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    fmpz_t shift;
    fmpz_init(shift);
    fmpz_set_si(shift, 1);

    fmpz_poly_t TrueResult;
    fmpz_poly_init(TrueResult);

    fmpz_poly_t trunc_poly;
    fmpz_poly_init(trunc_poly);
    fmpz_poly_t result;
    fmpz_poly_init(result);
    
    double *tabTps = malloc(sizeof(double) * 101);
    double *tabTps2 = malloc(sizeof(double) * 101);

    for (int i=0; i<=100; i++) {
        readPolyDATA(poly, 1, i);
        truncate_coefficients(trunc_poly, poly, poly->length);

        clock_t begin = clock();
        fmpz_poly_taylor_shift_divconquer(TrueResult, poly, shift);
        clock_t end = clock();
        tabTps[i] = (double)(end - begin); // / CLOCKS_PER_SEC;

        begin = clock();
        fmpz_poly_taylor_shift_divconquer(result, trunc_poly, shift);
        end = clock();
        tabTps2[i] = (double)(end - begin); // / CLOCKS_PER_SEC;

        if (same_signs(result, TrueResult))
            printf("yaayy, polys have same signs after taylor shift : %d\n", i); 
        else
            printf("polys don't have same signs after taylor shift :( : %d\n", i);
    }


    mkdir("src/EfficiencyTests/Results/Coeff_Truncation", 0777);
    FILE *fileResults;
    fileResults = fopen("src/EfficiencyTests/Results/Coeff_Truncation/Coeff_Truncation_ChangingCoeffSize.txt", "w");
    if (fileResults == NULL)
    {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }


    fprintf(fileResults, "Taylor Shift with truncated coefficient execution time (degree = 1500)\n"); // Title of the plot
    fprintf(fileResults, "No truncation\n");
    fprintTab(tabTps, 101, fileResults);
    fprintf(fileResults, "With truncation\n");
    fprintTab(tabTps2, 101, fileResults);

    
    fclose(fileResults);
    fmpz_poly_clear(TrueResult);
    fmpz_poly_clear(poly);
    fmpz_poly_clear(trunc_poly);
    fmpz_poly_clear(result);
    fmpz_clear(shift);
}





int main() {
    load_precomputed_polynomials(15);

    benchmark_degree();
    benchmark_coeffSize();
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/Coeff_Truncation/Coeff_Truncation_ChangingDegree.txt 'time' 0");
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/Coeff_Truncation/Coeff_Truncation_ChangingCoeffSize.txt 'time' 1");
    
    free_global_precomputed();
    return 0;
}