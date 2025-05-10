#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/bound.h"

// Enumeration for bound types
typedef enum {
    CAUCHY,
    LAGRANGE,
    MODIFIED_CAUCHY,
    FLINT
} bound_type_t;

int main() {
    slong i;
    fmpz_poly_t poly;
    fmpz_t fixedVariable;
    fmpz_t bound;

    
    fmpz_t tabBoundDegreeCauchy[101];

    fmpz_t tabBoundDegreeLagrange[101];

    fmpz_t tabBoundDegreeModifiedCauchy[101];

    fmpz_t tabBoundDegreeFlint[101];

    // To track the smallest bound and its type (using enum)
    fmpz_t smallestBound[101];
    bound_type_t smallestBoundType[101];

    fmpz_poly_init(poly);
    fmpz_init(fixedVariable);
    fmpz_init(bound);

    int lowest[4] = {0,0,0,0};

    for (i = 0; i <= 100; i++) {
        readPolyDATA(poly, 1, i);

        // Cauchy Bound
        Cauchy_bound(bound, poly);
        fmpz_set(tabBoundDegreeCauchy[i], bound);

        // Lagrange Bound
        Lagrange_bound(bound, poly);
        fmpz_set(tabBoundDegreeLagrange[i], bound);

        local_max_bound_implementation(bound, poly);
        fmpz_set(tabBoundDegreeModifiedCauchy[i], bound);

        
        fmpz_poly_bound_roots(bound, poly);
        fmpz_set(tabBoundDegreeFlint[i], bound);

        fmpz_init(smallestBound[i]);
        fmpz_set(smallestBound[i], tabBoundDegreeCauchy[i]);
        
        int min_bound=CAUCHY; // 0 cauchy, 1 lagrange, 2 modified, 3 flint

        if (fmpz_cmp(tabBoundDegreeLagrange[i], smallestBound[i]) < 0) {
            fmpz_set(smallestBound[i], tabBoundDegreeLagrange[i]);
            min_bound = LAGRANGE;
        }

        if (fmpz_cmp(tabBoundDegreeModifiedCauchy[i], smallestBound[i]) < 0) {
            fmpz_set(smallestBound[i], tabBoundDegreeModifiedCauchy[i]);
            min_bound = MODIFIED_CAUCHY;
        }

        if (fmpz_cmp(tabBoundDegreeFlint[i], smallestBound[i]) < 0) {
            fmpz_set(smallestBound[i], tabBoundDegreeFlint[i]);
            min_bound = FLINT;
        }

        lowest[min_bound]++;
    }

    for(int i=0; i<4; i++){
        printf("min %d, %d\n",i, lowest[i]);
    }


    fmpz_poly_clear(poly);
    fmpz_clear(fixedVariable);
    fmpz_clear(bound);
    for (i = 0; i <= 100; i++) {
        fmpz_clear(tabBoundDegreeCauchy[i]);
        fmpz_clear(tabBoundDegreeLagrange[i]);
        fmpz_clear(tabBoundDegreeModifiedCauchy[i]);
        fmpz_clear(tabBoundDegreeFlint[i]);
        fmpz_clear(smallestBound[i]);
    }

    return 0;
}