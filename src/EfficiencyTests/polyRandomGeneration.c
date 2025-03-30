#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <sys/stat.h>
#include <time.h>


void fmpz_poly_randtest_dense(fmpz_poly_t poly, flint_rand_t state, slong degree, flint_bitcnt_t coeff_size) {
    fmpz_t coeff;
    fmpz_init(coeff);

    for(int i=0; i<=degree; i++) {
        fmpz_randtest_not_zero(coeff, state, coeff_size);
        fmpz_poly_set_coeff_fmpz(poly, i, coeff);
    }

    fmpz_clear(coeff);
    return;
}



void writePolyFile(slong degree, flint_bitcnt_t coeff_size, int fixedVariable, int i) {
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, time(NULL)+i, 0);

    fmpz_poly_randtest_dense(poly, state, degree, coeff_size);

    FILE* filePoly;
    char pathFile[50];
    // fixedVariable decides which variable is fixed and which one is varying for tests
    if(!fixedVariable) {
        sprintf(pathFile, "DATA/Poly_ChangingDegree/poly%d.txt", i);
        filePoly = fopen(pathFile, "w");
    }
    else {
        sprintf(pathFile, "DATA/Poly_ChangingCoeffSize/poly%d.txt", i);
        filePoly = fopen(pathFile, "w");
    }

    if (filePoly == NULL) {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }

    fmpz_poly_fprint(filePoly, poly);
    fclose(filePoly);

    fmpz_poly_clear(poly);
    flint_randclear(state);
}



int main() {
    mkdir("DATA/Poly_ChangingDegree", 0777);
    mkdir("DATA/Poly_ChangingCoeffSize", 0777);

    slong maxDegree = 10000;
    flint_bitcnt_t maxCoeffSize = 100000;
    slong fixedDegree = 1500;
    flint_bitcnt_t fixedCoeffSize = 10000;
    
    int j=0;
    // generating polynomials with fixed coeff size and changing degree
    for(slong i=0; i<=maxDegree; i=i+100) {
        writePolyFile(i, fixedCoeffSize, 0, j);
        j++;
    }

    j=0;
    // generating polynomials with fixed degree and changing coeff size
    for(flint_bitcnt_t i=1; i<=(maxCoeffSize+1); i=i+1000) {         //cannot generate coeffs of size 0, so start at 1 
        writePolyFile(fixedDegree, i, 1, j);
        j++;
    }

    return 0;
}