#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/descartes.h"

char x = 'x';


void benchmarkDescartes(int fixedVariable, FILE* fileResults) {
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    clock_t start, end;
    double* tabTps = malloc(sizeof(double)*(101));  // car il y a 101 polynômes dans DATA

    for(slong i=0; i <= 100; i++) {
        readPolyDATA(poly, fixedVariable, i);

        start = clock();
        descartes_rule(poly);
        end = clock();
        tabTps[i] = (double)(end - start);
    }

    fprintTab(tabTps, 101, fileResults);

    fmpz_poly_clear(poly);
    free(tabTps);
}



int main()
{   
    mkdir("src/EfficiencyTests/Results/descartesTest", 0777);
    FILE* fileResults;

    fileResults = fopen("src/EfficiencyTests/Results/descartesTest/descartesChangingDegree.txt", "w");
    if (fileResults == NULL) {
        printf("The file is not opened. The program will "
            "now exit.\n");
        exit(0);
    }
    fprintf(fileResults, "Sign changes time efficiency (coeffSize=500)\n");    //Title of the plot
    fprintf(fileResults, "Implementation\n");    //labels of the plot
    benchmarkDescartes(0, fileResults);  //datas
    fclose(fileResults);


    fileResults = fopen("src/EfficiencyTests/Results/descartesTest/descartesChangingCoeffSize.txt", "w");
    if (fileResults == NULL) {
        printf("The file is not opened. The program will "
            "now exit.\n");
        exit(0);
    }
    fprintf(fileResults, "Sign changes time efficiency (degree=500)\n");    //Title of the plot
    fprintf(fileResults, "Implementation\n");    //labels of the plot
    benchmarkDescartes(1, fileResults);  //datas
    fclose(fileResults);

    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/descartesTest/descartesChangingDegree.txt 'time' 0");
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/descartesTest/descartesChangingCoeffSize.txt 'time' 1");

    return 0;
}