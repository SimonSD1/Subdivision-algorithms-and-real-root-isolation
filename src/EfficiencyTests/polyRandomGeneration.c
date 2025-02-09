#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include "../HeaderFiles/polyRandomGeneration_functions.h"



int main() {
    slong maxDegree = 20000;
    flint_bitcnt_t maxCoeffSize = 20000;
    slong fixedDegree = 500;        // 17min pour multiplication
    flint_bitcnt_t fixedCoeffSize = 500;    // 25min pour multiplication

    int j=0;
    // generating polynomials with fixed coeff size and changing degree
    for(slong i=0; i<=maxDegree; i=i+200) {
        writePolyFile(i, fixedCoeffSize, 0, j);
        j++;
    }

    j=0;
    // generating polynomials with fixed degree and changing coeff size
    for(flint_bitcnt_t i=1; i<=(maxCoeffSize+1); i=i+200) {         //cannot generate coeffs of size 0, so start at 1 
        writePolyFile(fixedDegree, i, 1, j);
        j++;
    }

    return 0;
}