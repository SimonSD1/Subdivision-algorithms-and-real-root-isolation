#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/descartes.h"

char x = 'x';

void testDescartes()
{
    fmpz_poly_t polynomials[4];
    for (int i = 0; i < 4; i++) {
        fmpz_poly_init(polynomials[i]);
    }

    fmpz_poly_set_str(polynomials[0], "4  2 -3 0 1");  
    fmpz_poly_set_str(polynomials[1], "3  1 -2 1");    
    fmpz_poly_set_str(polynomials[2], "3  1 0 -4");    
    fmpz_poly_set_str(polynomials[3], "3  1 -3 2"); 

    for (int i = 0; i < 4; i++) {
        printf("Polynomial: ");
        fmpz_poly_print(polynomials[i]);
        printf("\nNumber of sign changes: %ld\n", descartes_rule(polynomials[i]));
        fmpz_poly_clear(polynomials[i]);
    }
}

int main()
{
    testDescartes();
}