#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/bound.h"


int main()
{

    char x= 'x';

    // polynomials
    fmpz_poly_t P;

    // used for the randomness
    flint_rand_t state;

    // number of coefficients of the polynomial
    slong len;

    // number of bits per coefficient
    flint_bitcnt_t bits;

    len = 5;
    bits = 10;
    flint_randinit(state);
    // bug if time(NULL) used twice
    flint_randseed(state, time(NULL), 0);

    fmpz_poly_randtest(P, state, len, bits);

    printf("\npoly length=%ld\n",P->length);

    printf("\nPolynôme aléatoire :\n");
    fmpz_poly_print(P);
    printf("\n");
    fmpz_poly_print_pretty(P,&x);
    printf("\n\n");

    fmpz_t bound_flint;
    fmpz_init(bound_flint);

    fmpz_poly_bound_roots(bound_flint, P);

    printf("bound_flint= : ");
    fmpz_print(bound_flint);
    printf("\n");

    // fmpq_t bound;
    fmpz_t bound;
    fmpz_init(bound);
    Lagrange_bound(bound, P);

    printf("bound lagrange = : ");
    fmpz_print(bound);
    printf("\n");

    fmpz_t bound_pos;
    local_max_bound_implementation(bound_pos,P);
    printf("bound_pos= : ");
    fmpz_print(bound_pos);
    printf("\n");

    return 0;
}
