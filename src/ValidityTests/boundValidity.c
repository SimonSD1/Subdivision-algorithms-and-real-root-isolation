#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/bound.h"
#include "../HeaderFiles/functionsForTests.h"

#include <string.h>

slong descartes_rule(fmpz_poly_t poly)
{
    slong length = poly->length;

    slong sign_changes = 0;

    for (slong i = 1; i < length; i++)
    {
        if ((poly->coeffs[i] < 0 && poly->coeffs[i - 1] > 0) || (poly->coeffs[i] > 0 && poly->coeffs[i - 1] < 0))
        {
            sign_changes++;
        }
    }

    return sign_changes;
}

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

        //local_max_bound_implementation(bound, poly);
        fmpz_poly_bound_roots(bound,poly);

        fmpz_add_si(bound,bound,1);

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

int main()
{
    // to run or re-run the tests, pass argument : ./flintMultTest -runTests
    // if not it will only display the graphs

    printf("running");
    verifyBounds();

    return 0;
}
