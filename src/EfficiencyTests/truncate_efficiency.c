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
#include "../HeaderFiles/coeff_truncation.h"


void benchmark_TS_DivConq(slong maxLen, int fixedVariable, FILE *fileResults)
{
    fmpz_poly_t poly, trunc_poly, result, TrueResult;
    fmpz_poly_init(poly);
    fmpz_poly_init(result);
    fmpz_poly_init(TrueResult);
    fmpz_poly_init(trunc_poly);
    double *tabTps = malloc(sizeof(double) * (maxLen));
    double *tabTpsTrunc = malloc(sizeof(double) * (maxLen));
    fmpz_poly_t* power_array;

    for (slong i = 0; i < maxLen; i++)
    {
        flint_cleanup();
        power_array = NULL;
        readPolyDATA(poly, fixedVariable, i);
        truncate_coefficients(trunc_poly, poly);

        slong block_len, levels;
        compute_power_array(&power_array, poly, &block_len, &levels);

        
        clock_t begin = clock();
        iterative_taylor_shift_precompute(result, trunc_poly, power_array, block_len, levels);
        clock_t end = clock();
        tabTpsTrunc[i] = (double)(end - begin); // / CLOCKS_PER_SEC;

        begin = clock();
        iterative_taylor_shift_precompute(TrueResult, poly, power_array, block_len, levels);
        end = clock();
        tabTps[i] = (double)(end - begin);


        if (!same_signs(result, TrueResult)) {
            printf("polys don't have same signs after taylor shift :(\n");
            printf("iteration %ld, fixedVariable %d\n", i, fixedVariable);
        }
        
        free_precompute_table(power_array, levels);
    }


    // Title of the plot
    if(!fixedVariable)
        fprintf(fileResults, "Truncation before Taylor Shift time efficiencies (coeffSize = 10000)\n");
    else
        fprintf(fileResults, "Truncation before Taylor Shift time efficiencies (degree = 1500)\n");


    fprintf(fileResults, "No truncation\n");
    fprintTab(tabTps, maxLen, fileResults);
    fprintf(fileResults, "With truncation\n");
    fprintTab(tabTpsTrunc, maxLen, fileResults);


    fmpz_poly_clear(TrueResult);
    fmpz_poly_clear(poly);
    fmpz_poly_clear(trunc_poly);
    fmpz_poly_clear(result);
    free(tabTps);
    free(tabTpsTrunc);
}




int main(int argc, char *argv[])
{
    // to run or re-run the tests, pass argument : ./flintMultTest -r
    // if not it will only make the graphs as png files in the results folder

    if(argc > 1 && strcmp(argv[1], "-r") == 0) {            
        slong maxLen = 101;
        mkdir("src/EfficiencyTests/Results/Truncation", 0777);

        FILE *fileResults;
        fileResults = fopen("src/EfficiencyTests/Results/Truncation/Truncation_ChangingDegree.txt", "w");
        if (fileResults == NULL)
        {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        benchmark_TS_DivConq(maxLen, 0, fileResults);
        fclose(fileResults);


        fileResults = fopen("src/EfficiencyTests/Results/Truncation/Truncation_ChangingCoeffSize.txt", "w");
        if (fileResults == NULL)
        {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        benchmark_TS_DivConq(maxLen, 1, fileResults);
        fclose(fileResults);
    }

    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/Truncation/Truncation_ChangingDegree.txt 'time' 0");
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/Truncation/Truncation_ChangingCoeffSize.txt 'time' 1");



    return 0;
}
