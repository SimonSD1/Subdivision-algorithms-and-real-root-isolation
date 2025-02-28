#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/bound.h"
#include "../HeaderFiles/functionsForTests.h"
#include "../HeaderFiles/descartes.h"

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
    printf("Testing polynomial:\n");
    fmpz_poly_print_pretty(poly, "x");

    Lagrange_bound(bound, poly);
    printf("\nLagrange Bound: ");
    fmpz_print(bound);

    Cauchy_bound(bound, poly);
    printf("\nCauchy Bound: ");
    fmpz_print(bound);

    local_max_bound_implementation(bound, poly);
    printf("\nLocal Max Bound: ");
    fmpz_print(bound);
}

int main()
{
    //printf("running");
    //verifyBounds();

    fmpz_poly_t poly;
    fmpz_poly_init(poly);
    fmpz_t bound;
    fmpz_init(bound);

    // Define polynomial P(x) = x^4 - 16
    fmpz_poly_set_str(poly, "5  -16 0 0 0 1");

    testBounds(bound, poly);
    printf("\n\nSageMath Implementation //////////////\n");
    system("sage ValidityTests/bound.sage 'x^4 - 16'");

    fmpz_poly_clear(poly);
    fmpz_clear(bound);

    return 0;
}
