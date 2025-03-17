#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>


#ifndef TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H

void fprintTab(double* tab, int size, FILE* file);
void fprintFmpzTab(fmpz* tab, int size, FILE* file);
void readPolyDATA(fmpz_poly_t poly, int fixedVariable, int i, int LinearOrExponent);
#endif