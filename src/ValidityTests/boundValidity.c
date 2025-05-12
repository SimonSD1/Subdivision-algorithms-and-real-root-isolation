#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/bound.h"
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/descartes.h"

#define MAX_ROOTS 1000

void verifyBounds()
{
    fmpz_poly_t poly;
    fmpz_poly_init(poly);
    fmpz_t bound;

    // test for changing degree
    int fixedVariable = 0;
    for (slong i = 0; i <= 100; i++)
    {
        readPolyDATA(poly, fixedVariable, i);

        // local_max_bound_implementation(bound, poly);
        fmpz_poly_bound_roots(bound, poly);

        fmpz_add_si(bound, bound, 1);

        // shift of bound should result in no positive real roots
        fmpz_poly_taylor_shift_divconquer(poly, poly, bound);

        slong max_roots = descartes_rule(poly);

        // the number of roots is less than
        if (max_roots > 0)
        {
            printf("maybe bad bound\n");
        }
        else
        {
            printf("correct bound\n");
        }
    }

    fmpz_poly_clear(poly);
}


void testBounds(fmpz_t bound, fmpz_poly_t poly) {
    printf("Flint Implementation //////////////\n");
    printf("Testing polynomial: ");
    //fmpz_poly_print_pretty(poly, "x");

    Lagrange_bound(bound, poly);
    printf("\nLagrange Bound: ");
    fmpz_print(bound);

    Cauchy_bound(bound, poly);
    printf("\nCauchy Bound: ");
    fmpz_print(bound);

    local_max_bound_implementation(bound, poly);
    printf("\nLocal Max Bound: ");
    fmpz_print(bound);
    printf("\n");

    FILE* tmpFile;
    char pathFile[100];
    sprintf(pathFile, "src/bin/tmp_bound.txt");
    tmpFile = fopen(pathFile, "w");
    if (tmpFile == NULL) {
        perror("Erreur lors de l'ouverture du fichier");
        return;
    }
    fmpz_poly_fprint_pretty(tmpFile, poly, "x");
    fclose(tmpFile);
    int ret = system("\"/mnt/c/Program Files/Maple 2024/bin.X86_64_WINDOWS/cmaple.exe\"  src/ValidityTests/bound.mpl");
    if (ret != 0) {
        printf("Erreur lors de l'exécution du script Maple.\n");
        return;
    }
    
    sprintf(pathFile, "src/bin/maple_roots_bound.txt");
    tmpFile = fopen(pathFile, "r");
    if (!tmpFile) {
        perror("Failed to open Maple output file");
        exit(1);
    }

    double roots[MAX_ROOTS];
    int count = 0;
    while (fscanf(tmpFile, "%lf", &roots[count]) == 1) {
        count++;
        if (count >= MAX_ROOTS) break;
    }
    fclose(tmpFile);

    // Calcul du maximum
    double max_root = -1;
    for (int i = 0; i < count; i++) {
        if (roots[i] > max_root)
            max_root = roots[i];
    }

    // Affichage et comparaison
    printf("\nMax real positive root from Maple: %.15f\n", max_root);

    if(max_root > (slong)bound)
        printf("Bound implem invalid :(\n");
    else
        printf("Bound is valid yayyy\n");
}




int main()
{
    //printf("running");
    //verifyBounds();

    fmpz_poly_t poly;
    fmpz_poly_init(poly);
    fmpz_t bound;
    fmpz_init(bound);
    readPolyDATA(poly, 0, 2);

    testBounds(bound, poly);
    
    fmpz_poly_clear(poly);
    fmpz_clear(bound);

    return 0;
}
