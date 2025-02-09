#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/polyRandomGeneration_functions.h"


void fprintTab(double* tab, int size, FILE* file) {
    fprintf(file, "[");
    for(int i=0; i<size; i++) {
        if(i==(size-1))
            fprintf(file, "%.2f]\n", tab[i]);
        else
            fprintf(file, "%.2f, ", tab[i]);
    }
}



void readPolyDATA(fmpz_poly_t poly, int fixedVariable, int i) {
    FILE* filePoly;
    char pathFile[50];
    // fixedVariable decides which variable is fixed and which one is varying for tests
    if(!fixedVariable) {
        sprintf(pathFile, "../DATA/FixedCoeffSize_ChangingDegree/poly%d.txt", i);
        filePoly = fopen(pathFile, "r");
    }
    else {
        sprintf(pathFile, "../DATA/FixedDegree_ChangingCoeffSize/poly%d.txt", i);
        filePoly = fopen(pathFile, "r");
    }

    if (filePoly == NULL) {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }

    if(!fmpz_poly_fread(filePoly, poly))
        printf("Could not read the FILE");

    fclose(filePoly);
}



void timeEfficiencyMultiplication(void (*func)(fmpz_poly_t,const fmpz_poly_t,const fmpz_poly_t), int fixedVariable, FILE* fileResults) {
    fmpz_poly_t poly1, poly2, result;
    fmpz_poly_init(poly1);
    fmpz_poly_init(poly2);
    fmpz_poly_init(result);        

    clock_t start, end;
    double* tabTps = malloc(sizeof(double)*(101));  // car il y a 101 polynômes dans DATA

    for(slong i=0; i <= 100; i++) {
        readPolyDATA(poly1, fixedVariable, i);
        readPolyDATA(poly2, fixedVariable, i);  // je ne sais pas si c'est une bonne idée de tester la multiplaction sur 2 mêmes polynômes (carré du polynôme)

        start = clock();
        func(result, poly1, poly2);
        end = clock();
        tabTps[i] = (double)(end - start);
    }

    fprintTab(tabTps, 101, fileResults);

    fmpz_poly_clear(poly1);
    fmpz_poly_clear(poly2);
    fmpz_poly_clear(result);
}