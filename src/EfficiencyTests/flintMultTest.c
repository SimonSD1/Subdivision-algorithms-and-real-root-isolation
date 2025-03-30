#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include "../HeaderFiles/functionsForTests.h"


void timeEfficiencyMultiplication(void (*func)(fmpz_poly_t,const fmpz_poly_t,const fmpz_poly_t), int fixedVariable, FILE* fileResults) {
    fmpz_poly_t poly1, poly2, result;
    fmpz_poly_init(poly1);
    fmpz_poly_init(poly2);
    fmpz_poly_init(result);        

    clock_t start, end;
    double* tabTps = malloc(sizeof(double)*(101));  // car il y a 101 polynômes dans DATA

    for(slong i=0; i <= 100; i++) {
        readPolyDATA(poly1, fixedVariable, i);
        readPolyDATA(poly2, fixedVariable, i);  // je ne sais pas si c'est une bonne idée de tester la multiplaction sur 2 mêmes polynômes (carré du polynôme)

        start = clock();
        func(result, poly1, poly2);
        end = clock();
        tabTps[i] = (double)(end - start);
    }

    fprintTab(tabTps, 101, fileResults);

    fmpz_poly_clear(poly1);
    fmpz_poly_clear(poly2);
    fmpz_poly_clear(result);
    free(tabTps);
}





int main(int argc, char* argv[]) {
    // to run or re-run the tests, pass argument : ./flintMultTest -r
    // if not it will only make the graphs as png files in the results folder
    
    mkdir("src/EfficiencyTests/Results/flintMultTest", 0777);
    if(argc > 1 && strcmp(argv[1], "-r") == 0) {
        //note : on effectue les tests qu'une fois par taille, on ne fait pas de moyenne sur plusieurs tests pour une même taille

        FILE* fileResults;
        fileResults = fopen("src/EfficiencyTests/Results/flintMultTest/multChangingDegree.txt", "w");
        if (fileResults == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        fprintf(fileResults, "Flint multiplications time src/Efficiency (coeffSize=500)\n");    //Title of the plot
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


        fileResults = fopen("src/EfficiencyTests/Results/flintMultTest/multChangingCoeffSize.txt", "w");
        if (fileResults == NULL) {
            printf("The file is not opened. The program will "
                "now exit.\n");
            exit(0);
        }
        fprintf(fileResults, "Flint multiplications time src/Efficiency (degree=500)\n");    //Title of the plot
        fprintf(fileResults, "Classical\n");    //labels of the plot
        timeEfficiencyMultiplication(fmpz_poly_mul_classical, 1, fileResults);  //datas
        fprintf(fileResults, "Karatsuba\n");
        timeEfficiencyMultiplication(fmpz_poly_mul_karatsuba, 1, fileResults);
        fprintf(fileResults, "Kronecker Substitution\n");
        timeEfficiencyMultiplication(fmpz_poly_mul_KS, 1, fileResults);
        fprintf(fileResults, "Schonhage Strassen\n");
        timeEfficiencyMultiplication(fmpz_poly_mul_SS, 1, fileResults);
        fprintf(fileResults, "Adapted algorithm\n");
        timeEfficiencyMultiplication(fmpz_poly_mul, 1, fileResults);
        fclose(fileResults);
    }


    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/flintMultTest/multChangingDegree.txt 'time' 0");
    system("python3 src/EfficiencyTests/plotGenerator.py src/EfficiencyTests/Results/flintMultTest/multChangingCoeffSize.txt 'time' 1");

    return 0;
}