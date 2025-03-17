#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include "../HeaderFiles/functionsForTests.h"


void fprintTab(double* tab, int size, FILE* file) {
    fprintf(file, "[");
    for(int i=0; i<size; i++) {
        if(i==(size-1))
            fprintf(file, "%.2f]\n", tab[i]);
        else
            fprintf(file, "%.2f, ", tab[i]);
    }
}


void fprintFmpzTab(fmpz* tab, int size, FILE* file) {
    fprintf(file, "[");
    for(int i = 0; i < size; i++) {
        char *str = fmpz_get_str(NULL, 10, &tab[i]);  
        
        if(i == (size - 1)) {
            fprintf(file, "%s", str);
        } else {
            fprintf(file, "%s, ", str);
        }
        
        free(str);
    }
    fprintf(file, "]\n");
}



void readPolyDATA(fmpz_poly_t poly, int fixedVariable, int i, int LinearOrExponent) {
    FILE* filePoly;
    char pathFile[80];

    //LinearOrExponent decides whether we generate polynomials with sizes growing linearly or expenentially
    if(!LinearOrExponent) {
        // fixedVariable decides which variable is fixed and which one is varying for tests
        if(!fixedVariable) {
            sprintf(pathFile, "../DATA/LinearSizes_ChangingDegree/poly%d.txt", i);
            filePoly = fopen(pathFile, "r");
        }
        else {
            sprintf(pathFile, "../DATA/LinearSizes_ChangingCoeffSize/poly%d.txt", i);
            filePoly = fopen(pathFile, "r");
        }

        if (filePoly == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
    }
    else {
        if(!fixedVariable) {
            sprintf(pathFile, "../DATA/ExponentialSizes_ChangingDegree/poly%d.txt", i);
            filePoly = fopen(pathFile, "r");
        }
        else {
            sprintf(pathFile, "../DATA/ExponentialSizes_ChangingCoeffSize/poly%d.txt", i);
            filePoly = fopen(pathFile, "r");
        }

        if (filePoly == NULL) {
            printf("The file is not opened. The program will "
                "now exit. i=%d\n",i);
            exit(0);
        }
    }
    

    if(!fmpz_poly_fread(filePoly, poly))
        printf("Could not read the FILE");

    fclose(filePoly);
}