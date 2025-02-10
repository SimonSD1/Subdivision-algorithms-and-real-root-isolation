#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include "../HeaderFiles/polyRandomGeneration_functions.h"


#ifndef TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H

void fprintTab(double* tab, int size, FILE* file);
void readPolyDATA(fmpz_poly_t poly, int fixedVariable, int i);
void timeEfficiencyMultiplication(void (*func)(fmpz_poly_t,const fmpz_poly_t,const fmpz_poly_t), int fixedVariable, FILE* fileResults);

#endif