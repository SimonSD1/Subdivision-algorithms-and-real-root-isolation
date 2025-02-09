#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include "../HeaderFiles/functionsForTests.h"

//https://flintlib.org/doc/fmpz_poly.html

int main() {
    //note : on effectue les tests qu'une fois par taille, on ne fait pas de moyenne sur plusieurs tests pour une même taille
    FILE* fileResults;
    fileResults = fopen("EfficiencyTests/Results/multChangingDegree.txt", "w");
    if (fileResults == NULL) {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }
    fprintf(fileResults, "Flint multiplications time efficiency (coeffSize=500)\n");    //Title of the plot
    fprintf(fileResults, "Classical\n");    //labels of the plot
    timeEfficiencyMultiplication(fmpz_poly_mul_classical, 0, fileResults);  //datas
    fprintf(fileResults, "Karatsuba\n");
    timeEfficiencyMultiplication(fmpz_poly_mul_karatsuba, 0, fileResults);
    fprintf(fileResults, "Kronecker Substitution\n");
    timeEfficiencyMultiplication(fmpz_poly_mul_KS, 0, fileResults);
    fprintf(fileResults, "Schonhage Strassen\n");
    timeEfficiencyMultiplication(fmpz_poly_mul_SS, 0, fileResults);
    fprintf(fileResults, "Adapted algorithm\n");
    timeEfficiencyMultiplication(fmpz_poly_mul, 0, fileResults);

    fclose(fileResults);
    system("python.exe EfficiencyTests/plotGenerator.py EfficiencyTests/Results/multChangingDegree.txt 'time' 0");


    
    fileResults = fopen("EfficiencyTests/Results/multChangingCoeffSize.txt", "w");
    if (fileResults == NULL) {
        printf("The file is not opened. The program will "
               "now exit.\n");
        exit(0);
    }
    fprintf(fileResults, "Flint multiplications time efficiency (degree=500)\n");    //Title of the plot
    fprintf(fileResults, "Classical\n");    //labels of the plot
    timeEfficiencyMultiplication(fmpz_poly_mul_classical, 0, fileResults);  //datas
    fprintf(fileResults, "Karatsuba\n");
    timeEfficiencyMultiplication(fmpz_poly_mul_karatsuba, 0, fileResults);
    fprintf(fileResults, "Kronecker Substitution\n");
    timeEfficiencyMultiplication(fmpz_poly_mul_KS, 0, fileResults);
    fprintf(fileResults, "Schonhage Strassen\n");
    timeEfficiencyMultiplication(fmpz_poly_mul_SS, 0, fileResults);
    fprintf(fileResults, "Adapted algorithm\n");
    timeEfficiencyMultiplication(fmpz_poly_mul, 0, fileResults);

    fclose(fileResults);
    system("python.exe EfficiencyTests/plotGenerator.py EfficiencyTests/Results/multChangingCoeffSize.txt 'time' 1");

    return 0;
};