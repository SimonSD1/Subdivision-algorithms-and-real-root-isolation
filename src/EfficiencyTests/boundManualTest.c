#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/bound.h"

char x = 'x';

void testBound()
{
    fmpz_poly_t polynomials[4];
    for (int i = 0; i < 4; i++) {
        fmpz_poly_init(polynomials[i]);
    }

    fmpz_poly_set_str(polynomials[0], "4  2 -3 0 1");  
    fmpz_poly_set_str(polynomials[1], "3  1 10 -1");    
    fmpz_poly_set_str(polynomials[2], "3  1 10 -5");    
    fmpz_poly_set_str(polynomials[3], "3  1 -3 2"); 

    fmpz_t ub;
    fmpz_init(ub);

    for (int i = 0; i < 4; i++) {
        printf("Polynomial: ");
        fmpz_poly_print(polynomials[i]);
        

        local_max_bound_implementation(ub,polynomials[i]);

        printf("\nlocal max implementation : ");
        fmpz_print(ub);
        printf(" ");

        Cauchy_bound(ub,polynomials[i]);

        printf("cauchy: ");
        fmpz_print(ub);
        printf(" ");

        fmpz_poly_bound_roots(ub,polynomials[i]);

        printf("flint : ");
        fmpz_print(ub);
        printf(" \n");

    }
}

int main()
{
    testBound();
}