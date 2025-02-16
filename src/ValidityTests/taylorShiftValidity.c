#include <flint/flint.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <time.h>
#include <math.h>
#include "../HeaderFiles/taylorShift_implem.h"


int main()
{
    char x = 'x';
    // polynomials
    fmpz_poly_t P;
    fmpz_poly_init(P);

    // used for the randomness
    flint_rand_t state;

    // number of coefficients of the polynomial
    slong len;

    // number of bits per coefficient
    flint_bitcnt_t bits;

    len = 11;
    bits = 10;
    flint_randinit(state);
    // bug if time(NULL) used twice
    flint_randseed(state, time(NULL), 0);

    fmpz_poly_randtest(P, state, len, bits);
    while (P->length != len)
    {
        fmpz_poly_randtest(P, state, len, bits);
    }

    printf("\npoly length=%ld\n", P->length);

    printf("\nPolynôme aléatoire :\n");
    fmpz_poly_print(P);
    printf("\n");
    fmpz_poly_print_pretty(P, &x);
    printf("\n\n");

    slong a = 50;

    fmpz_t a_fmpz;
    fmpz_init_set_si(a_fmpz,a);

    fmpz_t a_flint;
    fmpz_init_set_ui(a_flint, a);

    fmpz_poly_t result;
    fmpz_poly_init(result);

    fmpz_poly_t result_flint;
    fmpz_poly_init(result_flint);

    poly_shift_plus_one(result, P, a_fmpz,0);

    printf("\nResult implem :\n");
    fmpz_poly_print(result);
    printf("\n");
    fmpz_poly_print_pretty(result, &x);
    printf("\n\n");

    fmpz_poly_taylor_shift_horner(result_flint, P, a_flint);
    printf("\nResult flint :\n");
    fmpz_poly_print(result_flint);
    printf("\n");
    fmpz_poly_print_pretty(result_flint, &x);
    printf("\n\n");

    fmpz_poly_clear(result);
    fmpz_poly_clear(result_flint);
    fmpz_poly_clear(P);
    return 0;
}
