#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>

//https://flintlib.org/doc/fmpz_poly.html



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

void printTab(double* tab, int size) {
    printf("[");
    for(int i=0; i<size; i++) {
        if(i==(size-1))
            printf("%.2f]\n", tab[i]);
        else
            printf("%.2f, ", tab[i]);
    }
}

void testMult_fixedCoeff(void (*func)(fmpz_poly_t,const fmpz_poly_t,const fmpz_poly_t), slong maxDegree, flint_bitcnt_t coeff_size, int nbTests) {
    fmpz_poly_t poly1, poly2, result;
    fmpz_poly_init(poly1);
    fmpz_poly_init(poly2);
    fmpz_poly_init(result);

    flint_rand_t state;
    flint_randinit(state);

    flint_randseed(state, time(NULL), 0);

    clock_t start, end;
    double tempsTotal = 0.0;
    double memTotal = 0.0;

    double* tabTps = malloc(sizeof(double)*(maxDegree + 1));
    double* tabMem = malloc(sizeof(double)*(maxDegree + 1));

    for(slong i=0; i <= maxDegree; i++) {
        for(int j=0; j<nbTests; j++) {
            fmpz_poly_randtest_dense(poly1, state, 1<<i, coeff_size);   //degree of 2**i
            fmpz_poly_randtest_dense(poly2, state, 1<<i, coeff_size);

            start = clock();
            struct mallinfo2 mem_before = mallinfo2();
            func(result, poly1, poly2);
            struct mallinfo2 mem_after = mallinfo2();
            end = clock();

            memTotal += (double)(mem_after.uordblks - mem_before.uordblks);
            tempsTotal += (double)(end - start);
        }
        tabTps[i] = tempsTotal/nbTests;
        tabMem[i] = memTotal/nbTests;
    }
    
    printf("Temps avec fixed size :\n");
    printTab(tabTps, maxDegree+1);
    printf("Memoire avec fixed size :\n");
    printTab(tabMem, maxDegree+1);
    printf("\n");

    fmpz_poly_clear(poly1);
    fmpz_poly_clear(poly2);
    fmpz_poly_clear(result);
    flint_randclear(state);
    return;
}


void testMult_fixedDeg(void (*func)(fmpz_poly_t,const fmpz_poly_t,const fmpz_poly_t), slong degree, flint_bitcnt_t maxCoeff_size, int nbTests) {
    fmpz_poly_t poly1, poly2, result;
    fmpz_poly_init(poly1);
    fmpz_poly_init(poly2);
    fmpz_poly_init(result);

    flint_rand_t state;
    flint_randinit(state);

    flint_randseed(state, time(NULL), 0);

    clock_t start, end;
    double tempsTotal = 0.0;
    double memTotal = 0.0;

    double* tabTps = malloc(sizeof(double)*(maxCoeff_size + 1));
    double* tabMem = malloc(sizeof(double)*(maxCoeff_size + 1));

    for(flint_bitcnt_t i=0; i <= maxCoeff_size; i++) {
        for(int j=0; j<nbTests; j++) {
            fmpz_poly_randtest_dense(poly1, state, degree, 1<<i);   //coeff size of 2**i bits
            fmpz_poly_randtest_dense(poly2, state, degree, 1<<i);

            start = clock();
            struct mallinfo2 mem_before = mallinfo2();
            func(result, poly1, poly2);
            struct mallinfo2 mem_after = mallinfo2();
            end = clock();

            memTotal += (double)(mem_after.uordblks - mem_before.uordblks);
            tempsTotal += (double)(end - start);
        }
        tabTps[i] = tempsTotal/nbTests;
        tabMem[i] = memTotal/nbTests;
    }

    printf("Temps avec fixed degree :\n");
    printTab(tabTps, maxCoeff_size+1);
    printf("Memoire avec fixed degree :\n");
    printTab(tabMem, maxCoeff_size+1);
    printf("\n");

    fmpz_poly_clear(poly1);
    fmpz_poly_clear(poly2);
    fmpz_poly_clear(result);
    flint_randclear(state);
    return;
}




int main() {
    int nbTest = 20; //moyenne sur nbTests tests
    flint_bitcnt_t fixedSize = 16;
    slong maxDeg = 16; //2^maxDeg

    flint_bitcnt_t maxSize = 16; //2^maxSiza
    slong deg = 16;

    printf("Classical :\n\n");
    testMult_fixedCoeff(fmpz_poly_mul_classical, maxDeg, fixedSize, nbTest);
    testMult_fixedDeg(fmpz_poly_mul_classical, deg, maxSize, nbTest);

    printf("\n\nKaratsuba :\n\n");
    testMult_fixedCoeff(fmpz_poly_mul_karatsuba, maxDeg, fixedSize, nbTest);
    testMult_fixedDeg(fmpz_poly_mul_karatsuba, deg, maxSize, nbTest);

    printf("\n\nAdapted algorithm :\n\n");
    testMult_fixedCoeff(fmpz_poly_mul, maxDeg, fixedSize, nbTest);
    testMult_fixedDeg(fmpz_poly_mul, deg, maxSize, nbTest);

    return 0;
}