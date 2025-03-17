#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <sys/stat.h>

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



void writePolyFile(slong degree, flint_bitcnt_t coeff_size, int fixedVariable, int i, int LinearOrExponent) {
    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    flint_rand_t state;
    flint_randinit(state);
    flint_randseed(state, time(NULL)+i, 0);

    fmpz_poly_randtest_dense(poly, state, degree, coeff_size);

    // Inside writePolyFile, before opening the file:
    mkdir("../DATA", 0777);  // Create DATA if not exists
    if (!fixedVariable && !LinearOrExponent) {
        mkdir("../DATA/LinearSizes_ChangingDegree", 0777);
    } else if (fixedVariable && !LinearOrExponent) {
        mkdir("../DATA/LinearSizes_ChangingCoeffSize", 0777);
    } else if (!fixedVariable && LinearOrExponent) {
        mkdir("../DATA/ExponentialSizes_ChangingDegree", 0777);
    } else {
        mkdir("../DATA/ExponentialSizes_ChangingCoeffSize", 0777);
    }

    FILE* filePoly;
    char pathFile[80];

    //LinearOrExponent decides whether we generate polynomials with sizes growing linearly or expenentially
    if(!LinearOrExponent) {

        // fixedVariable decides which variable is fixed and which one is varying for tests
        if(!fixedVariable) {
            sprintf(pathFile, "../DATA/LinearSizes_ChangingDegree/poly%d.txt", i);
            filePoly = fopen(pathFile, "w");
        }
        else {
            sprintf(pathFile, "../DATA/LinearSizes_ChangingCoeffSize/poly%d.txt", i);
            filePoly = fopen(pathFile, "w");
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
            filePoly = fopen(pathFile, "w");
        }
        else {
            sprintf(pathFile, "../DATA/ExponentialSizes_ChangingCoeffSize/poly%d.txt", i);
            filePoly = fopen(pathFile, "w");
        }

        if (filePoly == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }

    }
        
    fmpz_poly_fprint(filePoly, poly);
    fclose(filePoly);

    fmpz_poly_clear(poly);
    flint_randclear(state);
}


void generateLinear(slong maxDegree, flint_bitcnt_t maxCoeffSize, slong fixedDegree, flint_bitcnt_t fixedCoeffSize){
    int j=0;
    // generating polynomials with fixed coeff size and changing degree
    for(slong i=0; i<=maxDegree; i=i+200) {
        writePolyFile(i, fixedCoeffSize, 0, j, 0);
        j++;
    }

    j=0;
    // generating polynomials with fixed degree and changing coeff size
    for(flint_bitcnt_t i=1; i<=(maxCoeffSize+1); i=i+200) {         //cannot generate coeffs of size 0, so start at 1 
        writePolyFile(fixedDegree, i, 1, j, 0);
        j++;
    }
}

void generateExponent(int maxPower, int fixed){
    // generating polynomials with fixed coeff size and changing degree
    for(int i=0; i<=maxPower; i++) {
        slong degree = (slong)(1 << i)-1;
        writePolyFile(degree, (flint_bitcnt_t)fixed, 0, i, 1);  // sizes from 2^0 to 2^maxPower    
    }

    // generating polynomials with fixed degree and changing coeff size
    for(int i=0; i<=maxPower; i++) {
        flint_bitcnt_t coeff_size = (flint_bitcnt_t)1 << i;
        writePolyFile((slong)fixed, coeff_size, 1, i, 1);
    }
}

int main() {
    //generateLinear(20000, 20000, 500, 500);
    generateExponent(15, 2047);

    return 0;
}